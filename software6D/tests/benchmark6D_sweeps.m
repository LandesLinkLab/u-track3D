%% benchmark6D_sweeps.m
%  Multi-dimensional benchmark: at a baseline configuration, sweep ONE
%  parameter at a time across [wOrient, NoiseXY, maxSearchRadius,
%  timeWindow, Dtrans, density] and compare 3D vs 6D tracking.
%
%  Uses the NEW normalized cost function (normalizeSpatial=true), so
%  wOrient lives on a unit-ish scale (sweet spot ~0.3-1).
%
%  Computes both Jaqaman-style metrics (linkingAccuracy, splits, merges,
%  completeness, purity) AND Chenouard 2014 alpha / beta / JSC.
%
%  Output:
%    tests/benchmark_output/sweeps/sweeps_results.csv
%    tests/benchmark_output/sweeps/sweep_<axis>.png   (one per axis)
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM Multi-axis Sweep\n');
fprintf('============================================\n\n');

%% Setup
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));
addpath(genpath('C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D\software'));
if exist('trackCloseGapsKalmanSparse', 'file') ~= 2
    error('u-track3D not on path.');
end

outDir = fullfile(thisDir, 'benchmark_output', 'sweeps_v3');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Baseline configuration (kept fixed except for the swept axis)
B.density        = 30;
B.wOrient        = 0.3;          % sweet spot from extras_weights_heatmap (was 0.5)
B.nFrames        = 200;
B.fovHalfWidth   = 1500;          % nm  -> 3000x3000 nm
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
B.chenouardGate  = 100;          % nm matching threshold for Chenouard
B.gammaMean      = 0.7;          % wobble order parameter for baseline

seeds = 101:110;                 % 10 seeds for tighter CIs

%% Define sweeps: {label, fieldName, values}
sweeps = { ...
    {'wOrient',         'wOrient',         [0, 0.1, 0.3, 0.5, 1, 2, 5]}; ...
    {'density',         'density',         [5, 15, 30, 50]}; ...
    {'NoiseXY (nm)',    'noiseXY',         [10, 20, 50, 100]}; ...
    {'maxSearchRadius (nm)','maxSearchRadius', [200, 500, 1000, 2000]}; ...
    {'timeWindow (frames)', 'timeWindow', [1, 3, 5, 10, 20]}; ...
    {'Dtrans (nm^2/s)', 'dtrans',          [100, 500, 2000, 5000]}; ...
    {'\gamma (order param)', 'gammaMean',  [1.0, 0.85, 0.7, 0.5, 0.3, 0.1]}; ...
};

%% Aggregate results
allRows = [];
overallT0 = tic;

for iSweep = 1:length(sweeps)
    label    = sweeps{iSweep}{1};
    field    = sweeps{iSweep}{2};
    values   = sweeps{iSweep}{3};
    fprintf('\n========== Sweep: %s (%d values) ==========\n', label, length(values));

    for vIdx = 1:length(values)
        v = values(vIdx);

        % Determine which wOrients to test:
        %   wOrient sweep -> only that w value
        %   other sweeps  -> 3D (0) vs 6D (baseline w)
        if strcmp(field, 'wOrient')
            wList = v;
        else
            wList = [0, B.wOrient];
        end

        for seedIdx = 1:length(seeds)
            seed = seeds(seedIdx);

            % Apply the swept value to a fresh config
            cfg = B;
            cfg.(field) = v;

            %% Generate sim ONCE per (config, seed)
            [data, gtRaw] = generate6DTrajectories( ...
                'NumParticles', cfg.density, 'NumFrames', cfg.nFrames, ...
                'FrameTime', cfg.frameTime, ...
                'FieldOfView', [-cfg.fovHalfWidth cfg.fovHalfWidth ...
                                -cfg.fovHalfWidth cfg.fovHalfWidth], ...
                'ZRange', [-500 500], ...
                'Dtrans', cfg.dtrans, 'Drot', cfg.drot, ...
                'SpatialModel', 'brownian', 'RotationalModel', 'brownian', ...
                'GammaMean', cfg.gammaMean, 'GammaStd', 0.1, ...
                'Signal', cfg.signal, 'SignalStd', 300, ...
                'BlinkOffRate', cfg.blinkOff, 'BlinkOnRate', cfg.blinkOn, ...
                'BleachRate', 0.0005, ...
                'AddNoise', true, 'NoiseXY', cfg.noiseXY, 'NoiseZ', 50, ...
                'NoiseTheta', cfg.noiseTheta, 'NoisePhi', cfg.noiseTheta, ...
                'Seed', seed);

            movieInfo = convert6DSMOLMtoMovieInfo(data, ...
                'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);

            % Build GT matrix (frame, trackID, x, y, z, theta_rad, phi_rad, omega, signal)
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

            for wIdx = 1:length(wList)
                wO = wList(wIdx);

                %% --- Build cost matrices ---
                costMatrices = struct([]);
                costMatrices(1).funcName = 'costMat6DSMOLMLink';
                costMatrices(1).parameters = struct( ...
                    'linearMotion', 0, ...
                    'minSearchRadius', cfg.minSearchRadius, ...
                    'maxSearchRadius', cfg.maxSearchRadius, ...
                    'brownStdMult', cfg.brownStdMult, ...
                    'useLocalDensity', 1, 'nnWindow', 4, 'diagnostics', [], ...
                    'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
                    'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
                    'orientKalmanLevel', cfg.orientKalmanLevel, ...
                    'D_rot_prior', cfg.D_rot_prior, ...
                    'orientCRLB', cfg.orientCRLB, ...
                    'frameTime', cfg.frameTime, ...
                    'useSignalWeighting', false, ...
                    'normalizeSpatial', true, ...
                    'saveCostMatrix', false);

                costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
                costMatrices(2).parameters = struct( ...
                    'linearMotion', 0, ...
                    'minSearchRadius', cfg.minSearchRadius, ...
                    'maxSearchRadius', cfg.maxSearchRadius, ...
                    'brownStdMult', cfg.brownStdMult * ones(cfg.timeWindow,1), ...
                    'linStdMult', cfg.brownStdMult * ones(cfg.timeWindow,1), ...
                    'timeReachConfB', 3, 'timeReachConfL', 3, ...
                    'brownScaling', [0.5 0.01], 'linScaling', [0.5 0.01], ...
                    'useLocalDensity', 1, 'nnWindow', 4, 'lenForClassify', 5, ...
                    'maxAngleVV', 45, 'gapPenalty', 1.5, 'resLimit', 10, ...
                    'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
                    'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
                    'useBoScaling', true, ...
                    'normalizeSpatial', true, ...
                    'saveCostMatrix', false);

                gapCloseParam.timeWindow  = cfg.timeWindow;
                gapCloseParam.mergeSplit  = 0;
                gapCloseParam.minTrackLen = cfg.minTrackLen;
                kalmanFunctions.reserveMem  = 'kalmanResMemLM';
                kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
                kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
                kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

                tRun = tic;
                evalc('[tracksFinal,~,~]=trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,1);');
                runtime = toc(tRun);

                trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
                if isempty(trajRaw), continue; end

                % Jaqaman metrics
                jq = struct();
                evalc('[jq,~]=validateTracking6D(trajRaw,gtMat,''PlotResults'',false,''Verbose'',false);');
                % Chenouard metrics
                ch = chenouardMetrics(trajRaw, gtMat, 'GateDist', cfg.chenouardGate);
                % Conditional orientation error
                oe = orientationErrorPerTrack(trajRaw, gtMat, 'DistThreshold', 50);

                row = struct( ...
                    'sweep',          label, ...
                    'axisValue',      v, ...
                    'seed',           seed, ...
                    'wOrient',        wO, ...
                    'nGT',            jq.nGTParticles, ...
                    'nTracked',       jq.nTrackedTraj, ...
                    'linkAcc',        jq.linkingAccuracy, ...
                    'meanCompleteness', jq.meanCompleteness, ...
                    'meanPurity',     jq.meanPurity, ...
                    'nSplit',         jq.nSplit, ...
                    'nMerged',        jq.nMerged, ...
                    'alpha',          ch.alpha, ...
                    'beta',           ch.beta, ...
                    'JSC',            ch.JSC, ...
                    'nMatched_ch',    ch.nMatched, ...
                    'nMissed_ch',     ch.nMissed, ...
                    'nSpurious_ch',   ch.nSpurious, ...
                    'oriMeanDeg',     oe.meanErrDeg, ...
                    'oriMedianDeg',   oe.medianErrDeg, ...
                    'oriP95Deg',      oe.p95ErrDeg, ...
                    'oriNFrames',     oe.nLinkedFrames, ...
                    'runtime_s',      runtime);
                allRows = [allRows; row]; %#ok<AGROW>

                modeStr = '6D'; if wO == 0, modeStr = '3D'; end
                fprintf('  [%s=%g, seed=%d, %s w=%.2g] acc=%.1f%% pur=%.1f%% merge=%d alpha=%.2f oriMed=%.1f deg\n', ...
                    field, v, seed, modeStr, wO, ...
                    100*jq.linkingAccuracy, 100*jq.meanPurity, ...
                    jq.nMerged, ch.alpha, oe.medianErrDeg);
            end
        end
    end
end

fprintf('\n=== All sweeps complete in %.1fs ===\n', toc(overallT0));

%% Save CSV + MAT
T = struct2table(allRows);
csvPath = fullfile(outDir, 'sweeps_results.csv');
writetable(T, csvPath);
matPath = fullfile(outDir, 'sweeps_results.mat');
save(matPath, 'T', 'allRows', 'sweeps', 'B', 'seeds');
fprintf('Saved CSV: %s\n', csvPath);
fprintf('Saved MAT: %s\n', matPath);

%% ========================================================================
%  PLOTS — one figure per sweep with 6 metric panels
%  ========================================================================
metricCols  = {'linkAcc','meanPurity','nMerged','alpha','JSC','oriMedianDeg'};
metricLabel = {'Linking accuracy','Mean purity','# merges','Chenouard \alpha','JSC','Median orient. error (deg)'};
metricIsPct = [true, true, false, false, false, false];
metricRange = {[0 1.05], [0 1.05], [], [0 1.05], [0 1.05], []};

for iSweep = 1:length(sweeps)
    label    = sweeps{iSweep}{1};
    field    = sweeps{iSweep}{2};
    values   = sweeps{iSweep}{3};
    isWSweep = strcmp(field, 'wOrient');

    fig = figure('Position', [50 50 1600 700], 'Color', 'w', 'Visible', 'off');

    for iM = 1:6
        subplot(2, 3, iM); hold on; grid on; box on;
        colName = metricCols{iM};
        scale = 100 * metricIsPct(iM) + 1 * (~metricIsPct(iM));

        if isWSweep
            % Single curve: metric vs wOrient
            m = zeros(size(values)); s = zeros(size(values));
            for iv = 1:length(values)
                v = values(iv);
                msk = strcmp(T.sweep, label) & T.axisValue == v;
                xs = T.(colName)(msk);
                m(iv) = mean(xs) * scale;
                s(iv) = std(xs)/sqrt(max(length(xs),1)) * scale;
            end
            errorbar(values, m, s, '-o', 'Color', [0.6 0.15 0.5], ...
                'LineWidth', 2, 'MarkerFaceColor', [0.6 0.15 0.5], ...
                'MarkerSize', 8, 'DisplayName', 'normalized 6D');
            % Mark baseline w
            xline(B.wOrient, 'k--', sprintf('baseline w=%g', B.wOrient), ...
                'LabelOrientation','horizontal', 'Alpha', 0.5);
        else
            % Two curves: 3D (w=0) vs 6D (w=B.wOrient)
            m3 = zeros(size(values)); s3 = zeros(size(values));
            m6 = zeros(size(values)); s6 = zeros(size(values));
            for iv = 1:length(values)
                v = values(iv);
                m3rows = T.(colName)(strcmp(T.sweep,label) & T.axisValue==v & T.wOrient==0);
                m6rows = T.(colName)(strcmp(T.sweep,label) & T.axisValue==v & T.wOrient==B.wOrient);
                m3(iv) = mean(m3rows) * scale;
                s3(iv) = std(m3rows)/sqrt(max(length(m3rows),1)) * scale;
                m6(iv) = mean(m6rows) * scale;
                s6(iv) = std(m6rows)/sqrt(max(length(m6rows),1)) * scale;
            end
            errorbar(values, m3, s3, '-o', 'Color', [0.20 0.45 0.70], ...
                'LineWidth', 2, 'MarkerFaceColor', [0.20 0.45 0.70], ...
                'MarkerSize', 7, 'DisplayName', '3D (wOrient=0)');
            errorbar(values, m6, s6, '-s', 'Color', [0.85 0.33 0.10], ...
                'LineWidth', 2, 'MarkerFaceColor', [0.85 0.33 0.10], ...
                'MarkerSize', 7, 'DisplayName', sprintf('6D (wOrient=%g)', B.wOrient));
        end

        xlabel(label, 'FontWeight','bold', 'FontSize', 12);
        if metricIsPct(iM)
            ylabel([metricLabel{iM} ' (%)'], 'FontWeight','bold','FontSize', 12);
            ylim([0 105]);
        elseif ~isempty(metricRange{iM})
            ylabel(metricLabel{iM}, 'FontWeight','bold','FontSize', 12);
            ylim(metricRange{iM});
        else
            ylabel(metricLabel{iM}, 'FontWeight','bold','FontSize', 12);
        end
        title(metricLabel{iM}, 'FontSize', 13, 'FontWeight','bold');
        set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'TickDir','in');
        if values(end) > 10 * values(1)
            set(gca, 'XScale', 'log');
        end
        if iM == 1
            legend('Location','best','FontSize',10);
        end
    end
    sgtitle(sprintf('Sweep: %s   (baseline density=%d, wOrient=%g, %d seeds)', ...
        label, B.density, B.wOrient, length(seeds)), 'FontSize', 14, 'FontWeight','bold');

    safeName = regexprep(field, '\W', '_');
    pngPath = fullfile(outDir, sprintf('sweep_%s.png', safeName));
    exportgraphics(fig, pngPath, 'Resolution', 300);
    close(fig);
    fprintf('Plot saved: %s\n', pngPath);
end

%% Overall summary table
fprintf('\n=== Per-sweep summary (3D vs 6D averaged over axis) ===\n');
fprintf('  %-22s | %-10s %-10s %-8s %-8s\n', 'sweep', 'linkAcc 3D', 'linkAcc 6D', 'merge 3D', 'merge 6D');
fprintf('  %s\n', repmat('-', 1, 70));
for iSweep = 1:length(sweeps)
    label = sweeps{iSweep}{1};
    field = sweeps{iSweep}{2};
    if strcmp(field, 'wOrient'), continue; end
    msk3 = strcmp(T.sweep, label) & T.wOrient == 0;
    msk6 = strcmp(T.sweep, label) & T.wOrient == B.wOrient;
    fprintf('  %-22s | %-10.3f %-10.3f %-8.1f %-8.1f\n', label, ...
        mean(T.linkAcc(msk3)), mean(T.linkAcc(msk6)), ...
        mean(T.nMerged(msk3)), mean(T.nMerged(msk6)));
end
fprintf('\nDone.\n');
