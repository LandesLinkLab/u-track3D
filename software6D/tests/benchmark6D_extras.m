%% benchmark6D_extras.m
%  Two additional sweeps on top of benchmark6D_sweeps.m:
%    (A) Blinking density — 1D sweep over BlinkOffRate, 3D vs 6D
%    (B) wSpatial x wOrient heatmap — 2D grid showing the full cost-weight
%        landscape at the baseline regime
%
%  Output (in tests/benchmark_output/sweeps/):
%    extras_blinking.csv / .png    — blinking density sweep
%    extras_weights.csv / *.png    — weight grid heatmaps
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM Extras: blinking + weight grid\n');
fprintf('============================================\n\n');

%% Setup
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));
addpath(genpath('C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D\software'));
if exist('trackCloseGapsKalmanSparse', 'file') ~= 2
    error('u-track3D not on path.');
end

outDir = fullfile(thisDir, 'benchmark_output', 'sweeps');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Baseline (matches benchmark6D_sweeps.m)
B.density        = 30;
B.nFrames        = 200;
B.fovHalfWidth   = 1500;
B.dtrans         = 500;
B.drot           = 0.05;
B.signal         = 1500;
B.noiseXY        = 20;
B.noiseTheta     = 5;
B.blinkOff       = 0.30;
B.blinkOn        = 0.40;
B.minSearchRadius = 50;
B.maxSearchRadius = 500;
B.brownStdMult   = 3;
B.timeWindow     = 10;
B.minTrackLen    = 2;
B.orientKalmanLevel = 2;
B.D_rot_prior    = 0.05;
B.orientCRLB     = 5 * pi/180;
B.frameTime      = 0.05;
B.chenouardGate  = 100;

seeds = [101, 202, 303];

%% =====================================================================
%  (A) BLINKING DENSITY SWEEP
%  =====================================================================
fprintf('\n=== (A) Blinking-density sweep ===\n');
blinkOffVals = [0, 0.05, 0.15, 0.30, 0.50, 0.80];  % BlinkOnRate fixed at 0.40
% Steady-state on-fraction = onRate / (onRate + offRate)
ssOnFrac     = B.blinkOn ./ (B.blinkOn + blinkOffVals);

rowsBlink = [];
for vIdx = 1:length(blinkOffVals)
    bOff = blinkOffVals(vIdx);
    fprintf(' -- BlinkOffRate=%.2f (~%.0f%% on) --\n', bOff, 100*ssOnFrac(vIdx));
    for seedIdx = 1:length(seeds)
        seed = seeds(seedIdx);
        cfg = B; cfg.blinkOff = bOff;
        [data, gtRaw] = makeSim(cfg, seed);
        movieInfo = convert6DSMOLMtoMovieInfo(data, ...
            'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);
        gtMat = gtFromRaw(gtRaw);

        for wO = [0, 0.5]
            tRun = tic;
            tracksFinal = runTracker(movieInfo, cfg, 1.0, wO);
            runtime = toc(tRun);
            trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
            if isempty(trajRaw), continue; end

            jq = struct(); ch = struct();
            evalc('[jq,~]=validateTracking6D(trajRaw,gtMat,''PlotResults'',false,''Verbose'',false);');
            ch = chenouardMetrics(trajRaw, gtMat, 'GateDist', cfg.chenouardGate);

            rowsBlink = [rowsBlink; struct( ...
                'blinkOff', bOff, 'ssOnFrac', ssOnFrac(vIdx), ...
                'seed', seed, 'wOrient', wO, ...
                'linkAcc', jq.linkingAccuracy, ...
                'meanPurity', jq.meanPurity, ...
                'meanCompleteness', jq.meanCompleteness, ...
                'nSplit', jq.nSplit, 'nMerged', jq.nMerged, ...
                'alpha', ch.alpha, 'beta', ch.beta, 'JSC', ch.JSC, ...
                'runtime_s', runtime)]; %#ok<AGROW>
            modeStr = '3D'; if wO > 0, modeStr = sprintf('6D w=%.1f', wO); end
            fprintf('   [seed=%d %s] acc=%.1f%% pur=%.1f%% merge=%d alpha=%.3f\n', ...
                seed, modeStr, 100*jq.linkingAccuracy, 100*jq.meanPurity, jq.nMerged, ch.alpha);
        end
    end
end

Tblink = struct2table(rowsBlink);
writetable(Tblink, fullfile(outDir, 'extras_blinking.csv'));

%% Blinking plot
fig = figure('Position', [50 50 1600 700], 'Color', 'w', 'Visible', 'off');
metricCols  = {'linkAcc','meanPurity','nMerged','alpha','beta','JSC'};
metricLabel = {'Linking accuracy','Mean purity','# merges','Chenouard \alpha','Chenouard \beta','JSC'};
metricIsPct = [true, true, false, false, false, false];

for iM = 1:6
    subplot(2,3,iM); hold on; grid on; box on;
    col = metricCols{iM};
    scale = 100 * metricIsPct(iM) + 1 * (~metricIsPct(iM));
    for wO = [0, 0.5]
        m = zeros(size(blinkOffVals)); s = zeros(size(blinkOffVals));
        for iv = 1:length(blinkOffVals)
            v = Tblink.(col)(Tblink.blinkOff==blinkOffVals(iv) & Tblink.wOrient==wO);
            m(iv) = mean(v) * scale; s(iv) = std(v)/sqrt(max(length(v),1)) * scale;
        end
        if wO == 0
            color = [0.20 0.45 0.70]; mk = 'o'; nm = '3D';
        else
            color = [0.85 0.33 0.10]; mk = 's'; nm = sprintf('6D wOrient=%g', wO);
        end
        errorbar(ssOnFrac, m, s, ['-' mk], 'Color', color, ...
            'LineWidth', 2, 'MarkerFaceColor', color, 'MarkerSize', 7, ...
            'DisplayName', nm);
    end
    set(gca, 'XDir', 'reverse');  % left=more on, right=more sparse
    xlabel('Steady-state on-fraction', 'FontWeight','bold','FontSize', 12);
    if metricIsPct(iM)
        ylabel([metricLabel{iM} ' (%)'], 'FontWeight','bold','FontSize', 12);
        ylim([0 105]);
    else
        ylabel(metricLabel{iM}, 'FontWeight','bold','FontSize', 12);
    end
    title(metricLabel{iM}, 'FontSize', 13, 'FontWeight','bold');
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'TickDir','in');
    if iM == 1, legend('Location','best','FontSize',10); end
end
sgtitle(sprintf('Blinking-density sweep (BlinkOnRate=%.2f fixed, density=%d, 3 seeds)', ...
    B.blinkOn, B.density), 'FontSize', 14, 'FontWeight','bold');
exportgraphics(fig, fullfile(outDir, 'extras_blinking.png'), 'Resolution', 300);
close(fig);
fprintf('Plot saved: extras_blinking.png\n');

%% =====================================================================
%  (B) wSpatial x wOrient 2D HEATMAP
%  =====================================================================
fprintf('\n=== (B) wSpatial x wOrient heatmap ===\n');
wSpatVals = [0, 0.3, 1, 3];
wOriVals  = [0, 0.3, 1, 3];
% (0,0) is degenerate (all-zero cost). Skip it but include rest of grid.

rowsGrid = [];
nCells = length(wSpatVals) * length(wOriVals) - 1;
cellCount = 0;
gridT0 = tic;
for iS = 1:length(wSpatVals)
    for iO = 1:length(wOriVals)
        wS = wSpatVals(iS); wO = wOriVals(iO);
        if wS == 0 && wO == 0, continue; end
        cellCount = cellCount + 1;
        fprintf(' -- cell %d/%d: wSpatial=%g, wOrient=%g --\n', cellCount, nCells, wS, wO);

        for seedIdx = 1:length(seeds)
            seed = seeds(seedIdx);
            [data, gtRaw] = makeSim(B, seed);
            movieInfo = convert6DSMOLMtoMovieInfo(data, ...
                'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);
            gtMat = gtFromRaw(gtRaw);

            tRun = tic;
            tracksFinal = runTracker(movieInfo, B, wS, wO);
            runtime = toc(tRun);
            trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
            if isempty(trajRaw)
                fprintf('   seed=%d: NO TRACKS\n', seed);
                continue;
            end

            jq = struct(); ch = struct();
            evalc('[jq,~]=validateTracking6D(trajRaw,gtMat,''PlotResults'',false,''Verbose'',false);');
            ch = chenouardMetrics(trajRaw, gtMat, 'GateDist', B.chenouardGate);

            rowsGrid = [rowsGrid; struct( ...
                'wSpatial', wS, 'wOrient', wO, 'seed', seed, ...
                'linkAcc', jq.linkingAccuracy, ...
                'meanPurity', jq.meanPurity, ...
                'meanCompleteness', jq.meanCompleteness, ...
                'nSplit', jq.nSplit, 'nMerged', jq.nMerged, ...
                'alpha', ch.alpha, 'beta', ch.beta, 'JSC', ch.JSC, ...
                'runtime_s', runtime)]; %#ok<AGROW>
        end
    end
end
fprintf('\n=== Grid done in %.1fs ===\n', toc(gridT0));

Tgrid = struct2table(rowsGrid);
writetable(Tgrid, fullfile(outDir, 'extras_weights.csv'));

%% Build heatmaps: linkAcc, # merges, alpha, completeness
heatMetrics = {'linkAcc','nMerged','alpha','meanCompleteness'};
heatLabels  = {'Linking accuracy','# merges','Chenouard \alpha','Mean completeness'};
heatIsPct   = [true, false, false, true];
heatCmap    = {'parula', 'parula', 'parula', 'parula'};
heatInvert  = [false, true, false, false];  % # merges: low=good

fig = figure('Position', [50 50 1500 1100], 'Color', 'w', 'Visible', 'off');
for iM = 1:4
    H = nan(length(wOriVals), length(wSpatVals));   % rows=wOrient (y), cols=wSpatial (x)
    for iS = 1:length(wSpatVals)
        for iO = 1:length(wOriVals)
            v = Tgrid.(heatMetrics{iM})( ...
                Tgrid.wSpatial==wSpatVals(iS) & Tgrid.wOrient==wOriVals(iO));
            if isempty(v), continue; end
            H(iO, iS) = mean(v);
            if heatIsPct(iM), H(iO, iS) = 100 * H(iO, iS); end
        end
    end

    subplot(2,2,iM);
    imAlpha = ~isnan(H);
    imagesc(1:length(wSpatVals), 1:length(wOriVals), H, 'AlphaData', imAlpha);
    set(gca, 'YDir', 'normal');
    colormap(gca, heatCmap{iM});
    if heatInvert(iM), colormap(gca, flipud(parula)); end
    cb = colorbar; cb.Label.String = [heatLabels{iM} compose('%s', '')];
    if heatIsPct(iM), cb.Label.String = [heatLabels{iM} ' (%)']; end
    cb.Label.FontWeight = 'bold';
    cb.Label.FontSize = 11;

    % Annotate each cell with its value
    for iS = 1:length(wSpatVals)
        for iO = 1:length(wOriVals)
            if isnan(H(iO, iS))
                text(iS, iO, 'N/A', 'HorizontalAlignment','center', ...
                    'FontSize', 10, 'Color', [0.5 0.5 0.5]);
            else
                if heatIsPct(iM)
                    txt = sprintf('%.1f%%', H(iO, iS));
                else
                    txt = sprintf('%.2f', H(iO, iS));
                end
                % Choose text color for contrast
                txtCol = 'w'; rng = caxis;
                if H(iO, iS) > mean(rng), txtCol = 'k'; end
                if heatInvert(iM)
                    txtCol = 'k'; if H(iO,iS) > mean(rng), txtCol = 'w'; end
                end
                text(iS, iO, txt, 'HorizontalAlignment','center', ...
                    'FontSize', 11, 'FontWeight','bold','Color', txtCol);
            end
        end
    end

    set(gca, 'XTick', 1:length(wSpatVals), 'XTickLabel', arrayfun(@num2str, wSpatVals, 'UniformOutput', false));
    set(gca, 'YTick', 1:length(wOriVals),  'YTickLabel', arrayfun(@num2str, wOriVals,  'UniformOutput', false));
    xlabel('wSpatial', 'FontWeight','bold','FontSize', 12);
    ylabel('wOrient',  'FontWeight','bold','FontSize', 12);
    title(heatLabels{iM}, 'FontSize', 13, 'FontWeight','bold');
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'TickDir','out');
end
sgtitle(sprintf('Cost-weight grid (baseline density=%d, normalizeSpatial=true, %d seeds)', ...
    B.density, length(seeds)), 'FontSize', 14, 'FontWeight','bold');
exportgraphics(fig, fullfile(outDir, 'extras_weights_heatmap.png'), 'Resolution', 300);
close(fig);
fprintf('Plot saved: extras_weights_heatmap.png\n');

fprintf('\nDone.\n');

%% ========================================================================
%  HELPERS
%  ========================================================================
function [data, gtRaw] = makeSim(cfg, seed)
    [data, gtRaw] = generate6DTrajectories( ...
        'NumParticles', cfg.density, 'NumFrames', cfg.nFrames, ...
        'FrameTime', cfg.frameTime, ...
        'FieldOfView', [-cfg.fovHalfWidth cfg.fovHalfWidth ...
                        -cfg.fovHalfWidth cfg.fovHalfWidth], ...
        'ZRange', [-500 500], ...
        'Dtrans', cfg.dtrans, 'Drot', cfg.drot, ...
        'SpatialModel', 'brownian', 'RotationalModel', 'brownian', ...
        'GammaMean', 0.7, 'GammaStd', 0.1, ...
        'Signal', cfg.signal, 'SignalStd', 300, ...
        'BlinkOffRate', cfg.blinkOff, 'BlinkOnRate', cfg.blinkOn, ...
        'BleachRate', 0.0005, ...
        'AddNoise', true, 'NoiseXY', cfg.noiseXY, 'NoiseZ', 50, ...
        'NoiseTheta', cfg.noiseTheta, 'NoisePhi', cfg.noiseTheta, ...
        'Seed', seed);
end

function gtMat = gtFromRaw(gtRaw)
    allGT = vertcat(gtRaw.trajectories{:});
    gtMat = zeros(size(allGT,1), 9);
    gtMat(:,1) = allGT(:,1);
    gtMat(:,2) = allGT(:,9);
    gtMat(:,3) = allGT(:,3);
    gtMat(:,4) = allGT(:,4);
    gtMat(:,5) = allGT(:,5);
    gtMat(:,6) = allGT(:,6) * (pi/180);
    gtMat(:,7) = allGT(:,7) * (pi/180);
    g = max(0, min(1, allGT(:,8)));
    gtMat(:,8) = pi * (3 - sqrt(1 + 8*g));
    gtMat(:,9) = allGT(:,2);
    if size(allGT,2) >= 10
        gtMat = gtMat(allGT(:,10) == 0, :);
    end
end

function tracksFinal = runTracker(movieInfo, cfg, wS, wO)
    % wSpatial=0 needs an expanded search radius so the spatial gate
    % doesn't kill all candidates before orientation can vote.
    if wS == 0
        msr_min = 10; msr_max = 10000; bsm = 20;
    else
        msr_min = cfg.minSearchRadius; msr_max = cfg.maxSearchRadius; bsm = cfg.brownStdMult;
    end

    costMatrices = struct([]);
    costMatrices(1).funcName = 'costMat6DSMOLMLink';
    costMatrices(1).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', msr_min, 'maxSearchRadius', msr_max, ...
        'brownStdMult', bsm, ...
        'useLocalDensity', 1, 'nnWindow', 4, 'diagnostics', [], ...
        'wSpatial', wS, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'orientKalmanLevel', cfg.orientKalmanLevel, ...
        'D_rot_prior', cfg.D_rot_prior, ...
        'orientCRLB', cfg.orientCRLB, ...
        'frameTime', cfg.frameTime, ...
        'useSignalWeighting', false, ...
        'normalizeSpatial', true, 'saveCostMatrix', false);

    costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
    costMatrices(2).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', msr_min, 'maxSearchRadius', msr_max, ...
        'brownStdMult', bsm * ones(cfg.timeWindow,1), ...
        'linStdMult', bsm * ones(cfg.timeWindow,1), ...
        'timeReachConfB', 3, 'timeReachConfL', 3, ...
        'brownScaling', [0.5 0.01], 'linScaling', [0.5 0.01], ...
        'useLocalDensity', 1, 'nnWindow', 4, 'lenForClassify', 5, ...
        'maxAngleVV', 45, 'gapPenalty', 1.5, 'resLimit', 10, ...
        'wSpatial', wS, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'useBoScaling', true, ...
        'normalizeSpatial', true, 'saveCostMatrix', false);

    gapCloseParam.timeWindow  = cfg.timeWindow;
    gapCloseParam.mergeSplit  = 0;
    gapCloseParam.minTrackLen = cfg.minTrackLen;
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

    evalc('[tracksFinal,~,~]=trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,1);');
end
