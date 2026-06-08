%% benchmark6D.m
%  Head-to-head 3D-vs-6D tracking benchmark for software6D.
%
%  For each (density × seed), generates ONE ground-truth trajectory set,
%  then runs the u-track3D linker twice on identical detections:
%      (a) wOrient = 0    -> spatial-only baseline ("3D")
%      (b) wOrient = 1    -> orientation-augmented   ("6D")
%  and computes Jaqaman-style validation metrics. Uses orientation Kalman
%  Level 2 (adaptive Brownian) — the most principled level with no
%  multi-hypothesis bias (see costMat6DSMOLMLink.m header).
%
%  STRESS REGIME: heavy blinking creates frame gaps where multiple
%  candidate detections compete for the gap-closing linker. This is where
%  orientation information should help most. (Pure Brownian-no-blinking
%  is too easy for spatial linking at any realistic density.)
%
%  OUTPUT:
%      tests/benchmark_output/benchmark6D_results.csv
%      tests/benchmark_output/benchmark6D_comparison.png
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM Benchmark: 3D vs 6D\n');
fprintf('============================================\n\n');

%% ========================================================================
%  SETUP
%  ========================================================================
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));

utrack3DPath = 'C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D';
if exist(utrack3DPath, 'dir')
    addpath(genpath(fullfile(utrack3DPath, 'software')));
end
if exist('trackCloseGapsKalmanSparse', 'file') ~= 2
    error('u-track3D not on path; required for benchmark.');
end

outDir = fullfile(thisDir, 'benchmark_output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ========================================================================
%  SWEEP DESIGN
%  ========================================================================
densities    = [5, 15, 30];     % particles in a 3000x3000 nm field
% With normalizeSpatial=true (default in cost function), spatial cost
% is median-normalized to ~1 so wOrient lives on a unit-ish scale.
% Meaningful values: 0 (3D), ~0.3-1 (typical 6D), ~3+ (over-emphasize).
wOrientVals  = [0, 0.3, 1, 3];   % 0 = 3D baseline; others span 6D regime
seeds        = [101, 202, 303]; % replicates for error bars
nFrames      = 200;
fovHalfWidth = 1500;            % nm  -> 3000x3000 nm field (denser than default 6000)

% u-track3D parameters held FIXED across the sweep (so any benefit is
% attributable to orientation, not to retuning the spatial linker).
ut.minSearchRadius = 50;
ut.maxSearchRadius = 500;
ut.brownStdMult    = 3;
ut.timeWindow      = 10;   % gap closing horizon (frames)
ut.minTrackLen     = 2;

% Orientation Kalman params — Level 2 chosen on purpose (no min-bias).
ok.orientKalmanLevel = 2;
ok.D_rot_prior       = 0.05;
ok.orientCRLB        = 5 * pi/180;
ok.frameTime         = 0.05;

% Sim parameters held FIXED across density / wOrient sweep.
simBase = {
    'NumFrames',     nFrames, ...
    'FrameTime',     ok.frameTime, ...
    'FieldOfView',   [-fovHalfWidth fovHalfWidth -fovHalfWidth fovHalfWidth], ...
    'ZRange',        [-500 500], ...
    'Dtrans',        500, ...
    'Drot',          0.05, ...
    'SpatialModel',  'brownian', ...
    'RotationalModel','brownian', ...
    'GammaMean',     0.7, ...
    'GammaStd',      0.1, ...
    'Signal',        1500, ...
    'SignalStd',     300, ...
    'BlinkOffRate',  0.30, ...   % HEAVY blinking — stress gap closing
    'BlinkOnRate',   0.40, ...
    'BleachRate',    0.0005, ...
    'AddNoise',      true, ...
    'NoiseXY',       20, ...
    'NoiseZ',        50, ...
    'NoiseTheta',    5, ...
    'NoisePhi',      5
};

fprintf('Sweep grid:\n');
fprintf('  densities    : %s\n', mat2str(densities));
fprintf('  wOrient      : %s  (0 = 3D baseline, 1 = 6D)\n', mat2str(wOrientVals));
fprintf('  seeds        : %s\n', mat2str(seeds));
fprintf('  nFrames      : %d\n', nFrames);
fprintf('  FoV          : %d x %d nm\n', 2*fovHalfWidth, 2*fovHalfWidth);
fprintf('  Blinking     : OffRate=%.2f / OnRate=%.2f (heavy)\n', 0.30, 0.40);
fprintf('  Kalman level : %d (adaptive Brownian)\n', ok.orientKalmanLevel);
fprintf('  Search radius: min=%d max=%d nm, brownStdMult=%d\n', ...
    ut.minSearchRadius, ut.maxSearchRadius, ut.brownStdMult);
fprintf('  Gap window   : %d frames\n\n', ut.timeWindow);

nRuns = length(densities) * length(wOrientVals) * length(seeds);
fprintf('Total runs: %d\n\n', nRuns);

%% ========================================================================
%  RESULTS TABLE
%  ========================================================================
cols = {'density','seed','wOrient','nGT','nTracked','nRecovered', ...
        'nSplit','nMerged','linkingAccuracy','meanCompleteness', ...
        'meanPurity','nLinksTracked','nLinksCorrect','runtime_s'};
results = nan(nRuns, length(cols));
rowIdx = 0;

%% ========================================================================
%  MAIN SWEEP
%  ========================================================================
runT0 = tic;
for iDens = 1:length(densities)
    nPart = densities(iDens);

    for iSeed = 1:length(seeds)
        seed = seeds(iSeed);

        %% --- Generate sim ONCE per (density, seed) ---
        fprintf('\n=== density=%d, seed=%d ===\n', nPart, seed);
        tStage = tic;
        [data, gtRaw] = generate6DTrajectories(simBase{:}, ...
            'NumParticles', nPart, 'Seed', seed);
        movieInfo = convert6DSMOLMtoMovieInfo(data, ...
            'AngleUnits', 'degrees', 'WobbleInput', 'gamma', 'Verbose', false);
        fprintf('  sim+convert: %.1fs (%d detections across %d frames)\n', ...
            toc(tStage), size(data,1), length(movieInfo));

        % Build ground truth in validateTracking6D format:
        % [frame, trackID, x, y, z, theta_rad, phi_rad, omega, signal]
        allGT = vertcat(gtRaw.trajectories{:});
        gtMat = zeros(size(allGT,1), 9);
        gtMat(:,1) = allGT(:,1);                 % frame
        gtMat(:,2) = allGT(:,9);                 % trackID
        gtMat(:,3) = allGT(:,3);                 % x
        gtMat(:,4) = allGT(:,4);                 % y
        gtMat(:,5) = allGT(:,5);                 % z
        gtMat(:,6) = allGT(:,6) * (pi/180);      % theta deg -> rad
        gtMat(:,7) = allGT(:,7) * (pi/180);      % phi deg -> rad
        % column 8 of allGT is gamma; convert to omega via the same
        % quadratic inverse used in convert6DSMOLMtoMovieInfo.
        g = max(0, min(1, allGT(:,8)));
        gtMat(:,8) = pi * (3 - sqrt(1 + 8*g));   % omega
        gtMat(:,9) = allGT(:,2);                 % signal
        % drop interpolated GT rows (col 10 = isInterp)
        if size(allGT,2) >= 10
            gtMat = gtMat(allGT(:,10) == 0, :);
        end

        nGT = length(unique(gtMat(:,2)));

        for iW = 1:length(wOrientVals)
            wO = wOrientVals(iW);
            if wO == 0
                modeStr = '3D';
            else
                modeStr = sprintf('6D w=%g', wO);
            end
            fprintf('  [%s] wOrient=%.2f ... ', modeStr, wO);
            tTrack = tic;

            %% --- Linking cost matrix ---
            costMatrices = struct([]);
            costMatrices(1).funcName = 'costMat6DSMOLMLink';
            costMatrices(1).parameters = struct( ...
                'linearMotion',       0, ...
                'minSearchRadius',    ut.minSearchRadius, ...
                'maxSearchRadius',    ut.maxSearchRadius, ...
                'brownStdMult',       ut.brownStdMult, ...
                'useLocalDensity',    1, ...
                'nnWindow',           4, ...
                'diagnostics',        [], ...
                'wSpatial',           1.0, ...
                'wOrient',            wO, ...
                'wOmega',             0.0, ...
                'maxAngularDist',     pi/2, ...
                'useOrientation',     wO > 0, ...
                'orientKalmanLevel',  ok.orientKalmanLevel, ...
                'D_rot_prior',        ok.D_rot_prior, ...
                'orientCRLB',         ok.orientCRLB, ...
                'frameTime',          ok.frameTime, ...
                'useSignalWeighting', false, ...
                'saveCostMatrix',     false);

            %% --- Gap-closing cost matrix ---
            costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
            costMatrices(2).parameters = struct( ...
                'linearMotion',       0, ...
                'minSearchRadius',    ut.minSearchRadius, ...
                'maxSearchRadius',    ut.maxSearchRadius, ...
                'brownStdMult',       ut.brownStdMult * ones(ut.timeWindow,1), ...
                'linStdMult',         ut.brownStdMult * ones(ut.timeWindow,1), ...
                'timeReachConfB',     3, ...
                'timeReachConfL',     3, ...
                'brownScaling',       [0.5 0.01], ...
                'linScaling',         [0.5 0.01], ...
                'useLocalDensity',    1, ...
                'nnWindow',           4, ...
                'lenForClassify',     5, ...
                'maxAngleVV',         45, ...
                'gapPenalty',         1.5, ...
                'resLimit',           10, ...
                'wSpatial',           1.0, ...
                'wOrient',            wO, ...
                'wOmega',             0.0, ...
                'maxAngularDist',     pi/2, ...
                'useOrientation',     wO > 0, ...
                'useBoScaling',       true, ...
                'saveCostMatrix',     false);

            gapCloseParam.timeWindow  = ut.timeWindow;
            gapCloseParam.mergeSplit  = 0;
            gapCloseParam.minTrackLen = ut.minTrackLen;

            kalmanFunctions.reserveMem  = 'kalmanResMemLM';
            kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
            kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
            kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

            %% --- Run tracker (suppress its console chatter) ---
            evalc(sprintf('%s', ...
                '[tracksFinal, ~, ~] = trackCloseGapsKalmanSparse(movieInfo, costMatrices, gapCloseParam, kalmanFunctions, 3, 0, 1);'));

            %% --- Extract + validate ---
            trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
            if isempty(trajRaw)
                fprintf('NO TRACKS\n');
                continue;
            end

            % Quiet validate
            [metrics, ~] = evalc_validate(trajRaw, gtMat);

            runtime = toc(tTrack);
            rowIdx = rowIdx + 1;
            results(rowIdx, 1)  = nPart;
            results(rowIdx, 2)  = seed;
            results(rowIdx, 3)  = wO;
            results(rowIdx, 4)  = metrics.nGTParticles;
            results(rowIdx, 5)  = metrics.nTrackedTraj;
            results(rowIdx, 6)  = metrics.nRecovered;
            results(rowIdx, 7)  = metrics.nSplit;
            results(rowIdx, 8)  = metrics.nMerged;
            results(rowIdx, 9)  = metrics.linkingAccuracy;
            results(rowIdx, 10) = metrics.meanCompleteness;
            results(rowIdx, 11) = metrics.meanPurity;
            results(rowIdx, 12) = metrics.nLinksTracked;
            results(rowIdx, 13) = metrics.nLinksCorrect;
            results(rowIdx, 14) = runtime;

            fprintf('acc=%.2f%% comp=%.2f%% pur=%.2f%% split=%d merge=%d (%.1fs)\n', ...
                100*metrics.linkingAccuracy, 100*metrics.meanCompleteness, ...
                100*metrics.meanPurity, metrics.nSplit, metrics.nMerged, runtime);
        end
    end
end

results = results(1:rowIdx, :);
fprintf('\n=== Sweep complete in %.1fs ===\n', toc(runT0));

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================
T = array2table(results, 'VariableNames', cols);
csvPath = fullfile(outDir, 'benchmark6D_results.csv');
writetable(T, csvPath);
fprintf('Results saved: %s\n', csvPath);

matPath = fullfile(outDir, 'benchmark6D_results.mat');
save(matPath, 'T', 'results', 'cols', 'densities', 'wOrientVals', 'seeds', 'ut', 'ok', 'simBase');
fprintf('Workspace saved: %s\n', matPath);

%% ========================================================================
%  COMPARISON PLOT
%  ========================================================================
fig = figure('Position', [50 50 1300 380], 'Color', 'w', 'Visible', 'off');
metricCols = {'linkingAccuracy','meanCompleteness','meanPurity','nSplit'};
metricLabel= {'Linking accuracy','Mean completeness','Mean purity','# splits'};
metricFmt  = {'%','%','%','count'};

for iM = 1:4
    subplot(1, 4, iM); hold on; grid on; box on;
    colName = metricCols{iM};
    isPct = ~strcmp(metricFmt{iM}, 'count');
    scale = 100 * isPct + 1 * (~isPct);

    colors = lines(length(wOrientVals));
    markers = {'o','s','^','d','v','p'};
    for iW = 1:length(wOrientVals)
        wO = wOrientVals(iW);
        meansW = zeros(size(densities)); sesW = zeros(size(densities));
        for iD = 1:length(densities)
            v = T.(colName)(T.density==densities(iD) & T.wOrient==wO);
            meansW(iD) = mean(v) * scale;
            sesW(iD)   = std(v)/sqrt(length(v)) * scale;
        end
        mk = markers{1 + mod(iW-1, length(markers))};
        if wO == 0
            dispName = '3D (wOrient=0)';
        else
            dispName = sprintf('6D wOrient=%g', wO);
        end
        errorbar(densities, meansW, sesW, ['-' mk], 'Color', colors(iW,:), ...
            'LineWidth', 2, 'MarkerFaceColor', colors(iW,:), 'MarkerSize', 7, ...
            'DisplayName', dispName);
    end

    xlabel('# particles', 'FontWeight','bold', 'FontSize', 12);
    if isPct
        ylabel(sprintf('%s (%%)', metricLabel{iM}), 'FontWeight','bold','FontSize', 12);
        ylim([0 105]);
    else
        ylabel(metricLabel{iM}, 'FontWeight','bold','FontSize', 12);
    end
    title(metricLabel{iM}, 'FontSize', 13, 'FontWeight','bold');
    set(gca, 'XTick', densities, 'FontSize', 11, 'LineWidth', 1.5, 'TickDir','in');
    if iM == 1
        legend('Location','southwest','FontSize',10);
    end
end

sgtitle(sprintf('3D vs 6D tracking — heavy blinking, Level 2 Kalman (%d seeds/cond)', ...
    length(seeds)), 'FontSize', 14, 'FontWeight', 'bold');

pngPath = fullfile(outDir, 'benchmark6D_comparison.png');
exportgraphics(fig, pngPath, 'Resolution', 300);
fprintf('Plot saved: %s\n', pngPath);
close(fig);

%% ========================================================================
%  SUMMARY TABLE (3D vs 6D, averaged across seeds)
%  ========================================================================
fprintf('\n=== Summary: 3D vs 6D (mean over %d seeds) ===\n', length(seeds));
fprintf('  %-8s | %-12s | %-10s %-10s %-10s %-8s %-8s\n', ...
    'density','mode','linkAcc%','completen%','purity%','#split','#merge');
fprintf('  %s\n', repmat('-', 1, 82));
for iD = 1:length(densities)
    for iW = 1:length(wOrientVals)
        wO = wOrientVals(iW);
        if wO == 0
            modeLabel = '3D';
        else
            modeLabel = sprintf('6D w=%g', wO);
        end
        mask = T.density==densities(iD) & T.wOrient==wO;
        fprintf('  %-8d | %-12s | %-10.2f %-10.2f %-10.2f %-8.1f %-8.1f\n', ...
            densities(iD), modeLabel, ...
            100*mean(T.linkingAccuracy(mask)), ...
            100*mean(T.meanCompleteness(mask)), ...
            100*mean(T.meanPurity(mask)), ...
            mean(T.nSplit(mask)), mean(T.nMerged(mask)));
    end
end
fprintf('\nDone.\n');

%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================
function s = ternaryStr(cond, a, b)
    if cond, s = a; else, s = b; end
end

function [metrics, figH] = evalc_validate(trajRaw, gtMat)
    % validateTracking6D prints a lot — capture and discard.
    evalc('[metrics, figH] = validateTracking6D(trajRaw, gtMat, ''PlotResults'', false, ''Verbose'', false);');
end
