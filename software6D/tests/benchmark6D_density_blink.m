%% benchmark6D_density_blink.m
%  2D heatmap: density × blinking. At each cell, run 3D (wOrient=0) and
%  6D (wOrient=0.3) on the same GT, then plot 6D - 3D improvement for
%  each metric. Shows the regime where 6D's win is biggest.
%
%  Output (in tests/benchmark_output/sweeps_w0p3/):
%    density_blink.csv
%    density_blink_heatmap.png     (4 metrics × 6D-vs-3D delta)
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM: density x blinking heatmap\n');
fprintf('============================================\n\n');

%% Setup
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));
addpath(genpath('C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D\software'));

outDir = fullfile(thisDir, 'benchmark_output', 'sweeps_v3');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Baseline
B.nFrames        = 200;
B.fovHalfWidth   = 1500;
B.dtrans         = 500;
B.drot           = 0.05;
B.signal         = 1500;
B.noiseXY        = 20;
B.noiseTheta     = 5;
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
B.wOrient        = 0.3;

seeds = 101:110;  % 10 seeds for tighter CIs

%% Grid
densVals  = [5, 15, 30, 50];
blinkVals = [0, 0.15, 0.30, 0.60];   % BlinkOffRate. ssOn = onR/(onR+offR)
ssOn = B.blinkOn ./ (B.blinkOn + blinkVals);
fprintf('Density values: %s\n', mat2str(densVals));
fprintf('BlinkOff values: %s (ss-on: %s)\n\n', mat2str(blinkVals), mat2str(round(100*ssOn)));

rows = [];
totalT0 = tic;
nCells = length(densVals) * length(blinkVals);
cellCount = 0;
for iD = 1:length(densVals)
    for iB = 1:length(blinkVals)
        cellCount = cellCount + 1;
        cfg = B; cfg.density = densVals(iD); cfg.blinkOff = blinkVals(iB);
        fprintf(' cell %d/%d: density=%d, blinkOff=%.2f (%.0f%% on)\n', ...
            cellCount, nCells, cfg.density, cfg.blinkOff, 100*ssOn(iB));

        for seedIdx = 1:length(seeds)
            seed = seeds(seedIdx);
            [data, gtRaw] = makeSim(cfg, seed);
            movieInfo = convert6DSMOLMtoMovieInfo(data, ...
                'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);
            gtMat = gtFromRaw(gtRaw);

            for wO = [0, B.wOrient]
                tRun = tic;
                tracksFinal = runTracker(movieInfo, cfg, wO);
                runtime = toc(tRun);
                trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
                if isempty(trajRaw), continue; end
                jq = struct(); ch = struct();
                evalc('[jq,~]=validateTracking6D(trajRaw,gtMat,''PlotResults'',false,''Verbose'',false);');
                ch = chenouardMetrics(trajRaw, gtMat, 'GateDist', cfg.chenouardGate);
                oe = orientationErrorPerTrack(trajRaw, gtMat, 'DistThreshold', 50);

                rows = [rows; struct( ...
                    'density', cfg.density, 'blinkOff', cfg.blinkOff, ...
                    'ssOn', ssOn(iB), 'seed', seed, 'wOrient', wO, ...
                    'linkAcc', jq.linkingAccuracy, ...
                    'meanPurity', jq.meanPurity, ...
                    'meanCompleteness', jq.meanCompleteness, ...
                    'nSplit', jq.nSplit, 'nMerged', jq.nMerged, ...
                    'alpha', ch.alpha, 'beta', ch.beta, 'JSC', ch.JSC, ...
                    'oriMedianDeg', oe.medianErrDeg, 'oriP95Deg', oe.p95ErrDeg, ...
                    'runtime_s', runtime)]; %#ok<AGROW>
            end
        end
    end
end
fprintf('\n=== Grid done in %.1fs ===\n', toc(totalT0));

T = struct2table(rows);
writetable(T, fullfile(outDir, 'density_blink.csv'));
save(fullfile(outDir, 'density_blink.mat'), 'T', 'densVals', 'blinkVals', 'ssOn', 'B');
fprintf('Saved CSV / MAT.\n');

%% ========================================================================
%  HEATMAPS: 4 metrics × 3 panels each (3D / 6D / Δ)
%  ========================================================================
metricCols  = {'linkAcc','nMerged','alpha','meanCompleteness','oriMedianDeg'};
metricLabel = {'Linking accuracy','# merges','Chenouard \alpha','Mean completeness','Median orient. error (deg)'};
metricIsPct = [true, false, false, true, false];

fig = figure('Position', [40 40 1700 1450], 'Color', 'w', 'Visible', 'off');
tl = tiledlayout(5, 3, 'TileSpacing','compact', 'Padding','compact');

for iM = 1:5
    col = metricCols{iM};
    isPct = metricIsPct(iM);
    scale = 100 * isPct + 1 * (~isPct);

    H3 = nan(length(blinkVals), length(densVals));
    H6 = nan(length(blinkVals), length(densVals));
    for iD = 1:length(densVals)
        for iB = 1:length(blinkVals)
            v3 = T.(col)(T.density==densVals(iD) & T.blinkOff==blinkVals(iB) & T.wOrient==0);
            v6 = T.(col)(T.density==densVals(iD) & T.blinkOff==blinkVals(iB) & T.wOrient==B.wOrient);
            if ~isempty(v3), H3(iB,iD) = mean(v3)*scale; end
            if ~isempty(v6), H6(iB,iD) = mean(v6)*scale; end
        end
    end
    D = H6 - H3;

    % 3D panel
    nexttile;
    imagesc(1:length(densVals), 1:length(blinkVals), H3); set(gca,'YDir','normal');
    annotateHeatmap(H3, isPct);
    colormap(gca, parula); cb = colorbar; if isPct, cb.Label.String = '%'; end
    set(gca,'XTick',1:length(densVals),'XTickLabel',arrayfun(@num2str,densVals,'UniformOutput',false));
    set(gca,'YTick',1:length(blinkVals),'YTickLabel',arrayfun(@(x) sprintf('%.0f%%',100*x),ssOn,'UniformOutput',false));
    xlabel('# particles','FontWeight','bold'); ylabel('on-fraction','FontWeight','bold');
    title(sprintf('%s  (3D)', metricLabel{iM}), 'FontSize', 12, 'FontWeight','bold');
    set(gca,'FontSize',10,'LineWidth',1.2,'TickDir','out');

    % 6D panel
    nexttile;
    imagesc(1:length(densVals), 1:length(blinkVals), H6); set(gca,'YDir','normal');
    annotateHeatmap(H6, isPct);
    colormap(gca, parula); cb = colorbar; if isPct, cb.Label.String = '%'; end
    set(gca,'XTick',1:length(densVals),'XTickLabel',arrayfun(@num2str,densVals,'UniformOutput',false));
    set(gca,'YTick',1:length(blinkVals),'YTickLabel',arrayfun(@(x) sprintf('%.0f%%',100*x),ssOn,'UniformOutput',false));
    xlabel('# particles','FontWeight','bold'); ylabel('on-fraction','FontWeight','bold');
    title(sprintf('%s  (6D w=%.1f)', metricLabel{iM}, B.wOrient), 'FontSize', 12, 'FontWeight','bold');
    set(gca,'FontSize',10,'LineWidth',1.2,'TickDir','out');

    % Δ panel (6D - 3D)
    nexttile;
    imagesc(1:length(densVals), 1:length(blinkVals), D); set(gca,'YDir','normal');
    annotateHeatmap(D, isPct);
    cmax = max(abs(D(~isnan(D))));
    if isempty(cmax) || cmax == 0, cmax = 1; end
    caxis([-cmax cmax]);
    % red = bad for "lower is better" metrics (# merges); blue = good
    if strcmp(col, 'nMerged') || strcmp(col, 'oriMedianDeg')
        colormap(gca, flipud(redblue));   % blue = 6D better (lower is better)
    else
        colormap(gca, redblue);           % red = 6D better (higher is better)
    end
    cb = colorbar;
    if isPct, cb.Label.String = '\Delta %'; else, cb.Label.String = '\Delta'; end
    set(gca,'XTick',1:length(densVals),'XTickLabel',arrayfun(@num2str,densVals,'UniformOutput',false));
    set(gca,'YTick',1:length(blinkVals),'YTickLabel',arrayfun(@(x) sprintf('%.0f%%',100*x),ssOn,'UniformOutput',false));
    xlabel('# particles','FontWeight','bold'); ylabel('on-fraction','FontWeight','bold');
    title(sprintf('%s  (6D - 3D)', metricLabel{iM}), 'FontSize', 12, 'FontWeight','bold');
    set(gca,'FontSize',10,'LineWidth',1.2,'TickDir','out');
end
title(tl, sprintf('Density × blinking: 3D vs 6D (wOrient=%.1f, %d seeds)', B.wOrient, length(seeds)), ...
    'FontSize', 14, 'FontWeight','bold');
exportgraphics(fig, fullfile(outDir, 'density_blink_heatmap.png'), 'Resolution', 300);
close(fig);
fprintf('Plot saved: density_blink_heatmap.png\n');

fprintf('\nDone.\n');

%% ========================================================================
%  HELPERS
%  ========================================================================
function annotateHeatmap(H, isPct)
    for iR = 1:size(H,1)
        for iC = 1:size(H,2)
            v = H(iR,iC);
            if isnan(v), continue; end
            if isPct, txt = sprintf('%.1f%%', v);
            elseif abs(v) >= 10, txt = sprintf('%.0f', v);
            else, txt = sprintf('%.2f', v); end
            cax = caxis;
            txtCol = 'w'; if v > mean(cax), txtCol = 'k'; end
            text(iC, iR, txt, 'HorizontalAlignment','center', ...
                'FontSize', 10, 'FontWeight','bold','Color', txtCol);
        end
    end
end

function cm = redblue
    % Diverging red-white-blue colormap.
    n = 128;
    r = [linspace(0.95, 1, n)'; linspace(1, 0.05, n)'];
    g = [linspace(0.10, 1, n)'; linspace(1, 0.20, n)'];
    b = [linspace(0.20, 1, n)'; linspace(1, 0.55, n)'];
    cm = [r g b];
end

function [data, gtRaw] = makeSim(cfg, seed)
    [data, gtRaw] = generate6DTrajectories( ...
        'NumParticles', cfg.density, 'NumFrames', cfg.nFrames, ...
        'FrameTime', cfg.frameTime, ...
        'FieldOfView', [-cfg.fovHalfWidth cfg.fovHalfWidth ...
                        -cfg.fovHalfWidth cfg.fovHalfWidth], ...
        'ZRange', [-500 500], ...
        'Dtrans', cfg.dtrans, 'Drot', cfg.drot, ...
        'SpatialModel','brownian','RotationalModel','brownian', ...
        'GammaMean', 0.7, 'GammaStd', 0.1, ...
        'Signal', cfg.signal, 'SignalStd', 300, ...
        'BlinkOffRate', cfg.blinkOff, 'BlinkOnRate', cfg.blinkOn, ...
        'BleachRate', 0.0005, 'AddNoise', true, ...
        'NoiseXY', cfg.noiseXY, 'NoiseZ', 50, ...
        'NoiseTheta', cfg.noiseTheta, 'NoisePhi', cfg.noiseTheta, ...
        'Seed', seed);
end

function gtMat = gtFromRaw(gtRaw)
    allGT = vertcat(gtRaw.trajectories{:});
    gtMat = zeros(size(allGT,1), 9);
    gtMat(:,1)=allGT(:,1); gtMat(:,2)=allGT(:,9);
    gtMat(:,3)=allGT(:,3); gtMat(:,4)=allGT(:,4); gtMat(:,5)=allGT(:,5);
    gtMat(:,6)=allGT(:,6)*(pi/180); gtMat(:,7)=allGT(:,7)*(pi/180);
    g = max(0, min(1, allGT(:,8)));
    gtMat(:,8) = pi * (3 - sqrt(1 + 8*g));
    gtMat(:,9) = allGT(:,2);
    if size(allGT,2) >= 10, gtMat = gtMat(allGT(:,10) == 0, :); end
end

function tracksFinal = runTracker(movieInfo, cfg, wO)
    costMatrices = struct([]);
    costMatrices(1).funcName = 'costMat6DSMOLMLink';
    costMatrices(1).parameters = struct( ...
        'linearMotion', 0, 'minSearchRadius', cfg.minSearchRadius, ...
        'maxSearchRadius', cfg.maxSearchRadius, 'brownStdMult', cfg.brownStdMult, ...
        'useLocalDensity', 1, 'nnWindow', 4, 'diagnostics', [], ...
        'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'orientKalmanLevel', cfg.orientKalmanLevel, ...
        'D_rot_prior', cfg.D_rot_prior, 'orientCRLB', cfg.orientCRLB, ...
        'frameTime', cfg.frameTime, 'useSignalWeighting', false, ...
        'normalizeSpatial', true, 'saveCostMatrix', false);
    costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
    costMatrices(2).parameters = struct( ...
        'linearMotion', 0, 'minSearchRadius', cfg.minSearchRadius, ...
        'maxSearchRadius', cfg.maxSearchRadius, ...
        'brownStdMult', cfg.brownStdMult * ones(cfg.timeWindow,1), ...
        'linStdMult', cfg.brownStdMult * ones(cfg.timeWindow,1), ...
        'timeReachConfB', 3, 'timeReachConfL', 3, ...
        'brownScaling', [0.5 0.01], 'linScaling', [0.5 0.01], ...
        'useLocalDensity', 1, 'nnWindow', 4, 'lenForClassify', 5, ...
        'maxAngleVV', 45, 'gapPenalty', 1.5, 'resLimit', 10, ...
        'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'useBoScaling', true, 'normalizeSpatial', true, 'saveCostMatrix', false);
    gapCloseParam.timeWindow  = cfg.timeWindow;
    gapCloseParam.mergeSplit  = 0;
    gapCloseParam.minTrackLen = cfg.minTrackLen;
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    evalc('[tracksFinal,~,~]=trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,1);');
end
