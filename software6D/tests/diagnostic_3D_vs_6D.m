%% diagnostic_3D_vs_6D.m
%  Why are 3D and 6D producing identical tracks in benchmark6D?
%  Save the cost matrix for one config in each mode and compare directly.

clear; clc;
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));
addpath(genpath('C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D\software'));

outDir = fullfile(thisDir, 'benchmark_output', 'diagnostic');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Same sim params as benchmark6D, density=15 seed=101
[data, gtRaw] = generate6DTrajectories( ...
    'NumParticles', 15, 'NumFrames', 200, 'Seed', 101, ...
    'FrameTime', 0.05, ...
    'FieldOfView', [-1500 1500 -1500 1500], 'ZRange', [-500 500], ...
    'Dtrans', 500, 'Drot', 0.05, ...
    'SpatialModel', 'brownian', 'RotationalModel', 'brownian', ...
    'GammaMean', 0.7, 'GammaStd', 0.1, ...
    'Signal', 1500, 'SignalStd', 300, ...
    'BlinkOffRate', 0.30, 'BlinkOnRate', 0.40, 'BleachRate', 0.0005, ...
    'AddNoise', true, 'NoiseXY', 20, 'NoiseZ', 50, ...
    'NoiseTheta', 5, 'NoisePhi', 5);

movieInfo = convert6DSMOLMtoMovieInfo(data, ...
    'AngleUnits','degrees','WobbleInput','gamma','Verbose',false);

fprintf('Generated %d frames with detections:\n', length(movieInfo));
nDetPerFrame = arrayfun(@(s) size(s.xCoord,1), movieInfo);
fprintf('  nDetections range: [%d, %d], mean %.1f\n', min(nDetPerFrame), max(nDetPerFrame), mean(nDetPerFrame));

%% Common params
ut.minSearchRadius = 50;  ut.maxSearchRadius = 500;
ut.brownStdMult = 3;  ut.timeWindow = 10;

linkBase = struct( ...
    'linearMotion', 0, 'minSearchRadius', ut.minSearchRadius, ...
    'maxSearchRadius', ut.maxSearchRadius, 'brownStdMult', ut.brownStdMult, ...
    'useLocalDensity', 1, 'nnWindow', 4, 'diagnostics', [], ...
    'wSpatial', 1.0, 'wOmega', 0.0, 'maxAngularDist', pi/2, ...
    'orientKalmanLevel', 2, 'D_rot_prior', 0.05, ...
    'orientCRLB', 5*pi/180, 'frameTime', 0.05, ...
    'useSignalWeighting', false, ...
    'saveCostMatrix', true, 'costMatSaveFrame', 5);

gapBase = struct( ...
    'linearMotion', 0, 'minSearchRadius', ut.minSearchRadius, ...
    'maxSearchRadius', ut.maxSearchRadius, ...
    'brownStdMult', ut.brownStdMult * ones(ut.timeWindow,1), ...
    'linStdMult', ut.brownStdMult * ones(ut.timeWindow,1), ...
    'timeReachConfB', 3, 'timeReachConfL', 3, ...
    'brownScaling', [0.5 0.01], 'linScaling', [0.5 0.01], ...
    'useLocalDensity', 1, 'nnWindow', 4, 'lenForClassify', 5, ...
    'maxAngleVV', 45, 'gapPenalty', 1.5, 'resLimit', 10, ...
    'wSpatial', 1.0, 'wOmega', 0.0, 'maxAngularDist', pi/2, ...
    'useBoScaling', true, 'saveCostMatrix', false);

gapCloseParam.timeWindow  = ut.timeWindow;
gapCloseParam.mergeSplit  = 0;
gapCloseParam.minTrackLen = 2;
kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% --- 3D run ---
fprintf('\n=== 3D RUN (wOrient=0) ===\n');
costMatrices3D = struct([]);
costMatrices3D(1).funcName = 'costMat6DSMOLMLink';
p3 = linkBase; p3.wOrient = 0; p3.useOrientation = false;
p3.costMatSavePath = fullfile(outDir, 'cm_3D');
if ~exist(p3.costMatSavePath, 'dir'), mkdir(p3.costMatSavePath); end
costMatrices3D(1).parameters = p3;
costMatrices3D(2).funcName = 'costMat6DSMOLMCloseGaps';
g3 = gapBase; g3.wOrient = 0; g3.useOrientation = false;
costMatrices3D(2).parameters = g3;
[tracks3D, ~, ~] = trackCloseGapsKalmanSparse(movieInfo, costMatrices3D, gapCloseParam, kalmanFunctions, 3, 0, 1);
fprintf('  -> %d tracks\n', length(tracks3D));

%% --- 6D run (wOrient=1, the benchmark default) ---
fprintf('\n=== 6D RUN (wOrient=1) ===\n');
costMatrices6D = struct([]);
costMatrices6D(1).funcName = 'costMat6DSMOLMLink';
p6 = linkBase; p6.wOrient = 1; p6.useOrientation = true;
p6.costMatSavePath = fullfile(outDir, 'cm_6D');
if ~exist(p6.costMatSavePath, 'dir'), mkdir(p6.costMatSavePath); end
costMatrices6D(1).parameters = p6;
costMatrices6D(2).funcName = 'costMat6DSMOLMCloseGaps';
g6 = gapBase; g6.wOrient = 1; g6.useOrientation = true;
costMatrices6D(2).parameters = g6;
[tracks6D, ~, ~] = trackCloseGapsKalmanSparse(movieInfo, costMatrices6D, gapCloseParam, kalmanFunctions, 3, 0, 1);
fprintf('  -> %d tracks\n', length(tracks6D));

%% --- 6D run with BIG weight (wOrient=10000) — was the workaround pre-normalization ---
fprintf('\n=== 6D RUN (wOrient=10000, pre-normalization workaround) ===\n');
costMatrices6Dbig = struct([]);
costMatrices6Dbig(1).funcName = 'costMat6DSMOLMLink';
pBig = linkBase; pBig.wOrient = 10000; pBig.useOrientation = true;
pBig.normalizeSpatial = false;  % old behavior
pBig.costMatSavePath = fullfile(outDir, 'cm_6Dbig');
if ~exist(pBig.costMatSavePath, 'dir'), mkdir(pBig.costMatSavePath); end
costMatrices6Dbig(1).parameters = pBig;
costMatrices6Dbig(2).funcName = 'costMat6DSMOLMCloseGaps';
gBig = gapBase; gBig.wOrient = 10000; gBig.useOrientation = true;
gBig.normalizeSpatial = false;
costMatrices6Dbig(2).parameters = gBig;
[tracks6Dbig, ~, ~] = trackCloseGapsKalmanSparse(movieInfo, costMatrices6Dbig, gapCloseParam, kalmanFunctions, 3, 0, 1);
fprintf('  -> %d tracks (unnormalized w=10000)\n', length(tracks6Dbig));

%% --- 6D run with NORMALIZED spatial cost + wOrient=1 ---
fprintf('\n=== 6D RUN (wOrient=1 + normalizeSpatial=true) ===\n');
costMatricesN = struct([]);
costMatricesN(1).funcName = 'costMat6DSMOLMLink';
pN = linkBase; pN.wOrient = 1; pN.useOrientation = true;
pN.normalizeSpatial = true;
pN.costMatSavePath = fullfile(outDir, 'cm_6Dnorm');
if ~exist(pN.costMatSavePath, 'dir'), mkdir(pN.costMatSavePath); end
costMatricesN(1).parameters = pN;
costMatricesN(2).funcName = 'costMat6DSMOLMCloseGaps';
gN = gapBase; gN.wOrient = 1; gN.useOrientation = true;
gN.normalizeSpatial = true;
costMatricesN(2).parameters = gN;
[tracksN, ~, ~] = trackCloseGapsKalmanSparse(movieInfo, costMatricesN, gapCloseParam, kalmanFunctions, 3, 0, 1);
fprintf('  -> %d tracks (normalized w=1)\n', length(tracksN));
fprintf('  isequaln(tracksN, tracks3D): %d  (should be 0 — confirms orientation acts)\n', isequaln(tracksN, tracks3D));

%% --- Compare cost matrices ---
fprintf('\n=== COST MATRIX COMPARISON (frame 5) ===\n');
cm3 = load(fullfile(outDir, 'cm_3D', 'costMatrix_linking.mat'));
cm6 = load(fullfile(outDir, 'cm_6D', 'costMatrix_linking.mat'));

C3 = full(cm3.costMatFinal); C3(C3==0) = NaN;
C6 = full(cm6.costMatFinal); C6(C6==0) = NaN;

fprintf('  3D cost matrix: size %dx%d, non-NaN cells %d\n', size(C3,1), size(C3,2), sum(~isnan(C3(:))));
fprintf('  3D cost range: [%.3f, %.3f], mean %.3f\n', ...
    min(C3(~isnan(C3))), max(C3(~isnan(C3))), mean(C3(~isnan(C3))));
fprintf('  6D cost matrix: size %dx%d, non-NaN cells %d\n', size(C6,1), size(C6,2), sum(~isnan(C6(:))));
fprintf('  6D cost range: [%.3f, %.3f], mean %.3f\n', ...
    min(C6(~isnan(C6))), max(C6(~isnan(C6))), mean(C6(~isnan(C6))));

% Element-wise difference
diff = C6 - C3;
fprintf('  Element-wise diff (6D - 3D): min %.3f, max %.3f, mean %.3f, |.|max %.3f\n', ...
    min(diff(~isnan(diff))), max(diff(~isnan(diff))), ...
    mean(diff(~isnan(diff))), max(abs(diff(~isnan(diff)))));

% Does orientation cost ever change the argmin per row?
[~, argmin3] = min(C3, [], 2, 'omitnan');
[~, argmin6] = min(C6, [], 2, 'omitnan');
nChanged = sum(argmin3 ~= argmin6);
fprintf('  Rows whose argmin changed: %d / %d (%.1f%%)\n', ...
    nChanged, length(argmin3), 100*nChanged/length(argmin3));

%% --- Compare track outputs directly ---
fprintf('\n=== TRACK OUTPUT COMPARISON ===\n');
fprintf('  3D nTracks: %d\n', length(tracks3D));
fprintf('  6D nTracks: %d\n', length(tracks6D));

n3 = sum(arrayfun(@(t) numel(t.tracksFeatIndxCG), tracks3D));
n6 = sum(arrayfun(@(t) numel(t.tracksFeatIndxCG), tracks6D));
fprintf('  3D total featIndx entries: %d\n', n3);
fprintf('  6D total featIndx entries: %d\n', n6);

% Compare track lengths
len3 = sort(arrayfun(@(t) size(t.tracksCoordAmpCG, 2)/8, tracks3D));
len6 = sort(arrayfun(@(t) size(t.tracksCoordAmpCG, 2)/8, tracks6D));
fprintf('  3D track length distribution: min=%d med=%d max=%d\n', min(len3), median(len3), max(len3));
fprintf('  6D track length distribution: min=%d med=%d max=%d\n', min(len6), median(len6), max(len6));

% Are the track structures byte-for-byte identical?
isIdent = isequaln(tracks3D, tracks6D);
fprintf('  isequaln(tracks3D, tracks6D): %d\n', isIdent);

if ~isIdent
    fprintf('\n  -> Tracks DIFFER but produce same validation metrics. Investigate validator.\n');
else
    fprintf('\n  -> Tracks are IDENTICAL. Investigate cost matrix (above).\n');
end
