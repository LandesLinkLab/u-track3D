%% benchmark6D_mc.m
%  Latin-hypercube Monte Carlo benchmark. Samples N points from a
%  15-dimensional joint space of (sim + algorithm) parameters, runs both
%  3D (wOrient=0) and 6D (wOrient sampled) on the SAME ground truth at
%  each point, and saves all metrics.
%
%  Output:
%    tests/benchmark_output/mc/mc_results.csv     (one row per (sample, mode))
%    tests/benchmark_output/mc/mc_design.csv      (LHS design matrix)
%    tests/benchmark_output/mc/mc_results.mat     (full MATLAB workspace)
%
%  Designed to run in background (~70 min for N=300).
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

function benchmark6D_mc()
clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM MC benchmark (Latin hypercube)\n');
fprintf('============================================\n\n');

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));
addpath(genpath('C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D\software'));
if exist('trackCloseGapsKalmanSparse', 'file') ~= 2
    error('u-track3D not on path.');
end

outDir = fullfile(thisDir, 'benchmark_output', 'mc');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ========================================================================
%  PARAMETER SPACE
%  Each row: {name, min, max, scale ('lin'|'log'|'int'|'choice'),
%             optional choice cell}
%  ========================================================================
P = { ...
    {'NumParticles',     5,    50,   'int'}, ...
    {'Dtrans',           100,  10000,'log'}, ...        % nm^2/s
    {'Drot',             0.005,0.5,  'log'}, ...        % rad^2/s
    {'GammaMean',        0.10, 1.00, 'lin'}, ...        % order param
    {'BlinkOffRate',     0.00, 0.80, 'lin'}, ...
    {'BlinkOnRate',      0.20, 1.00, 'lin'}, ...
    {'Signal',           200,  5000, 'log'}, ...        % photons
    {'NoiseXY',          5,    100,  'log'}, ...        % nm
    {'NoiseTheta',       1,    20,   'log'}, ...        % deg
    {'wOrient',          0.05, 5.0,  'log'}, ...        % 6D only; 3D forces 0
    {'maxSearchRadius',  200,  2000, 'log'}, ...        % nm
    {'brownStdMult',     2,    5,    'lin'}, ...
    {'timeWindow',       2,    20,   'int'}, ...        % frames
    {'orientKalmanLevel',1,    3,    'int'}, ...        % {1,2,3}
    {'orientCRLB',       2,    15,   'log'} ...         % deg
};
nVars = length(P);
nSamples = 300;
fixedSeed = 12345;   % seed for the LHS itself (reproducible design)

% Fixed sim parameters (anything not in P)
F.nFrames        = 200;
F.fovHalfWidth   = 1500;
F.frameTime      = 0.05;
F.gammaStd       = 0.1;
F.signalStd      = 300;
F.noiseZ         = 50;       % nm
F.bleachRate     = 0.0005;
F.spatialModel   = 'brownian';
F.rotationalModel= 'brownian';
F.D_rot_prior    = 0.05;
F.minSearchRadius= 50;
F.useLocalDensity= 1;
F.nnWindow       = 4;
F.mergeSplit     = 0;
F.minTrackLen    = 2;
F.gapPenalty     = 1.5;
F.normalizeSpatial = true;
F.useSignalWeighting = false;
F.chenouardGate  = 100;

% Each sim/run uses ONE seed for the simulator. Use the sample index as
% the simulation seed so re-running with a different LHS sample size still
% reuses comparable seeds where they align.
seedOffset = 7000;

%% ========================================================================
%  GENERATE LHS DESIGN
%  ========================================================================
rng(fixedSeed, 'twister');
U = lhsdesign(nSamples, nVars, 'criterion', 'maximin', 'iterations', 5);

design = cell(nSamples, nVars);
for v = 1:nVars
    pv = P{v};
    lo = pv{2}; hi = pv{3}; scale = pv{4};
    u  = U(:, v);
    switch scale
        case 'lin'
            x = lo + u .* (hi - lo);
        case 'log'
            x = exp(log(lo) + u .* (log(hi) - log(lo)));
        case 'int'
            x = round(lo + u .* (hi - lo));
        case 'choice'
            x = pv{5}(min(length(pv{5}), 1 + floor(u .* length(pv{5}))));
        otherwise
            error('Unknown scale %s for variable %s', scale, pv{1});
    end
    design(:, v) = num2cell(x);
end

% Save design matrix
varNames = cellfun(@(p) p{1}, P, 'UniformOutput', false);
Tdesign = cell2table(design, 'VariableNames', varNames);
Tdesign.sample = (1:nSamples)';
writetable(Tdesign, fullfile(outDir, 'mc_design.csv'));
fprintf('LHS design saved: %s  (N=%d, D=%d)\n\n', ...
    fullfile(outDir, 'mc_design.csv'), nSamples, nVars);

%% ========================================================================
%  RUN SAMPLES
%  ========================================================================
allRows = [];
overallT0 = tic;
for iS = 1:nSamples
    cfg = F;
    for v = 1:nVars
        cfg.(varNames{v}) = design{iS, v};
    end
    seed = seedOffset + iS;

    fprintf('--- Sample %d/%d (seed=%d) ---\n', iS, nSamples, seed);
    fprintf('   density=%d, blink=(%.2f,%.2f), signal=%.0f, noiseXY=%.1f, noiseTheta=%.2f\n', ...
        cfg.NumParticles, cfg.BlinkOffRate, cfg.BlinkOnRate, cfg.Signal, ...
        cfg.NoiseXY, cfg.NoiseTheta);
    fprintf('   Dtrans=%.1f, Drot=%.4f, gamma=%.2f, wOrient=%.3f, MSR=%.0f, tWin=%d, level=%d\n', ...
        cfg.Dtrans, cfg.Drot, cfg.GammaMean, cfg.wOrient, ...
        cfg.maxSearchRadius, cfg.timeWindow, cfg.orientKalmanLevel);

    try
        [data, gtRaw] = makeSim(cfg, seed);
        movieInfo = convert6DSMOLMtoMovieInfo(data, ...
            'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);
        if isempty(movieInfo) || all(arrayfun(@(s) isempty(s.xCoord), movieInfo))
            fprintf('   *** no detections, skip ***\n'); continue;
        end
        gtMat = gtFromRaw(gtRaw);

        for modeIdx = 1:2
            % mode 1 = 3D baseline (wOrient forced to 0)
            % mode 2 = 6D (wOrient as sampled)
            if modeIdx == 1
                wOrient = 0; modeName = '3D';
            else
                wOrient = cfg.wOrient; modeName = '6D';
            end
            tRun = tic;
            tracksFinal = runTracker(movieInfo, cfg, wOrient);
            runtime = toc(tRun);

            trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', false);
            if isempty(trajRaw)
                fprintf('   [%s] no tracks\n', modeName); continue;
            end
            jq = struct(); ch = struct(); oe = struct();
            evalc('[jq,~]=validateTracking6D(trajRaw,gtMat,''PlotResults'',false,''Verbose'',false);');
            ch = chenouardMetrics(trajRaw, gtMat, 'GateDist', cfg.chenouardGate);
            oe = orientationErrorPerTrack(trajRaw, gtMat, 'DistThreshold', 50);

            row = struct( ...
                'sample',          iS, ...
                'seed',            seed, ...
                'mode',            modeName, ...
                'wOrientApplied',  wOrient, ...
                'NumParticles',    cfg.NumParticles, ...
                'Dtrans',          cfg.Dtrans, ...
                'Drot',            cfg.Drot, ...
                'GammaMean',       cfg.GammaMean, ...
                'BlinkOffRate',    cfg.BlinkOffRate, ...
                'BlinkOnRate',     cfg.BlinkOnRate, ...
                'Signal',          cfg.Signal, ...
                'NoiseXY',         cfg.NoiseXY, ...
                'NoiseTheta',      cfg.NoiseTheta, ...
                'wOrientSampled',  cfg.wOrient, ...
                'maxSearchRadius', cfg.maxSearchRadius, ...
                'brownStdMult',    cfg.brownStdMult, ...
                'timeWindow',      cfg.timeWindow, ...
                'orientKalmanLevel', cfg.orientKalmanLevel, ...
                'orientCRLB',      cfg.orientCRLB, ...
                'nGT',             jq.nGTParticles, ...
                'nTracked',        jq.nTrackedTraj, ...
                'linkAcc',         jq.linkingAccuracy, ...
                'meanCompleteness',jq.meanCompleteness, ...
                'meanPurity',      jq.meanPurity, ...
                'nSplit',          jq.nSplit, ...
                'nMerged',         jq.nMerged, ...
                'alpha',           ch.alpha, ...
                'beta',            ch.beta, ...
                'JSC',             ch.JSC, ...
                'oriMedianDeg',    oe.medianErrDeg, ...
                'oriP95Deg',       oe.p95ErrDeg, ...
                'oriNFrames',      oe.nLinkedFrames, ...
                'runtime_s',       runtime);
            allRows = [allRows; row]; %#ok<AGROW>
            fprintf('   [%s] acc=%.1f%% pur=%.1f%% merge=%d alpha=%.2f oriMed=%.1f deg (%.1fs)\n', ...
                modeName, 100*jq.linkingAccuracy, 100*jq.meanPurity, ...
                jq.nMerged, ch.alpha, oe.medianErrDeg, runtime);
        end
    catch ME
        fprintf('   *** sample %d failed: %s ***\n', iS, ME.message);
    end

    % Save incrementally every 25 samples so we don't lose results if it crashes
    if mod(iS, 25) == 0 && ~isempty(allRows)
        T = struct2table(allRows);
        writetable(T, fullfile(outDir, 'mc_results.csv'));
        save(fullfile(outDir, 'mc_results.mat'), 'T', 'P', 'F', 'allRows', 'Tdesign', 'nSamples');
        fprintf('   ... checkpoint saved after sample %d (%.1f min elapsed)\n', iS, toc(overallT0)/60);
    end
end

%% ========================================================================
%  FINAL SAVE
%  ========================================================================
T = struct2table(allRows);
writetable(T, fullfile(outDir, 'mc_results.csv'));
save(fullfile(outDir, 'mc_results.mat'), 'T', 'P', 'F', 'allRows', 'Tdesign', 'nSamples');
fprintf('\n=== MC complete: %d samples × 2 modes in %.1f min ===\n', ...
    nSamples, toc(overallT0)/60);
fprintf('CSV: %s\n', fullfile(outDir, 'mc_results.csv'));
fprintf('MAT: %s\n', fullfile(outDir, 'mc_results.mat'));
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function [data, gtRaw] = makeSim(cfg, seed)
    [data, gtRaw] = generate6DTrajectories( ...
        'NumParticles', cfg.NumParticles, 'NumFrames', cfg.nFrames, ...
        'FrameTime', cfg.frameTime, ...
        'FieldOfView', [-cfg.fovHalfWidth cfg.fovHalfWidth ...
                        -cfg.fovHalfWidth cfg.fovHalfWidth], ...
        'ZRange', [-500 500], ...
        'Dtrans', cfg.Dtrans, 'Drot', cfg.Drot, ...
        'SpatialModel', cfg.spatialModel, 'RotationalModel', cfg.rotationalModel, ...
        'GammaMean', cfg.GammaMean, 'GammaStd', cfg.gammaStd, ...
        'Signal', cfg.Signal, 'SignalStd', cfg.signalStd, ...
        'BlinkOffRate', cfg.BlinkOffRate, 'BlinkOnRate', cfg.BlinkOnRate, ...
        'BleachRate', cfg.bleachRate, ...
        'AddNoise', true, 'NoiseXY', cfg.NoiseXY, 'NoiseZ', cfg.noiseZ, ...
        'NoiseTheta', cfg.NoiseTheta, 'NoisePhi', cfg.NoiseTheta, ...
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
    tw = cfg.timeWindow;
    costMatrices = struct([]);
    costMatrices(1).funcName = 'costMat6DSMOLMLink';
    costMatrices(1).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', cfg.minSearchRadius, ...
        'maxSearchRadius', cfg.maxSearchRadius, ...
        'brownStdMult', cfg.brownStdMult, ...
        'useLocalDensity', cfg.useLocalDensity, 'nnWindow', cfg.nnWindow, ...
        'diagnostics', [], ...
        'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'orientKalmanLevel', cfg.orientKalmanLevel, ...
        'D_rot_prior', cfg.D_rot_prior, ...
        'orientCRLB', cfg.orientCRLB * pi/180, ...
        'frameTime', cfg.frameTime, ...
        'useSignalWeighting', cfg.useSignalWeighting, ...
        'normalizeSpatial', cfg.normalizeSpatial, 'saveCostMatrix', false);
    costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
    costMatrices(2).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', cfg.minSearchRadius, ...
        'maxSearchRadius', cfg.maxSearchRadius, ...
        'brownStdMult', cfg.brownStdMult * ones(tw,1), ...
        'linStdMult',   cfg.brownStdMult * ones(tw,1), ...
        'timeReachConfB', 3, 'timeReachConfL', 3, ...
        'brownScaling', [0.5 0.01], 'linScaling', [0.5 0.01], ...
        'useLocalDensity', cfg.useLocalDensity, 'nnWindow', cfg.nnWindow, ...
        'lenForClassify', 5, 'maxAngleVV', 45, ...
        'gapPenalty', cfg.gapPenalty, 'resLimit', 10, ...
        'wSpatial', 1.0, 'wOrient', wO, 'wOmega', 0.0, ...
        'maxAngularDist', pi/2, 'useOrientation', wO > 0, ...
        'useBoScaling', true, ...
        'normalizeSpatial', cfg.normalizeSpatial, 'saveCostMatrix', false);
    gapCloseParam.timeWindow  = tw;
    gapCloseParam.mergeSplit  = cfg.mergeSplit;
    gapCloseParam.minTrackLen = cfg.minTrackLen;
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    evalc('[tracksFinal,~,~]=trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,1);');
end
