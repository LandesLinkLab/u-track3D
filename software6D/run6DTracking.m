%% run6DTracking.m
%  Main wrapper script for the 6D-SMOLM tracking pipeline.
%
%  Runs all stages:
%    Stage 1: Preprocessing — convert fusion output to movieInfo
%    Stage 2: Frame-to-frame linking (Cost Matrix 1)
%    Stage 3: Gap closing with photoblinking awareness (Cost Matrix 2)
%    Stage 4: Post-processing — trajectory filtering and combination
%    Stage 5: Diffusion analysis (MLE with variable time lags)
%
%  USAGE:
%    1. Set the parameters in the USER CONFIGURATION section below.
%    2. Run this script.
%
%  REQUIREMENTS:
%    - u-track3D must be on the MATLAB path:
%      addpath(genpath('path/to/u-track3D/software'))
%    - software6D must be on the MATLAB path:
%      addpath(genpath('path/to/software6D'))
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('=== 6D-SMOLM Tracking Pipeline ===\n');
fprintf('Started: %s\n\n', datestr(now));

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

% --- Input data ---
% Path to fusion output (CSV or MAT file)
inputFile = '';  % <-- SET THIS: e.g., 'C:\...\fusion_output.csv'

% Column mapping for the input file (Fusion output format)
% Fusion columns: [frame, signal, x_nm, y_nm, z_nm, theta_deg, phi_deg, gamma]
colMap = struct('frame',1, 'intensity',2, 'x_nm',3, 'y_nm',4, ...
    'z_nm',5, 'theta',6, 'phi',7, 'gamma',8);

% --- Output directory ---
outputDir = '';  % <-- SET THIS: e.g., 'C:\...\tracking_results'

% --- Imaging parameters ---
frameTime = 0.05;  % Exposure time per frame (seconds). Adjust to your camera.

% --- Stage 1: Preprocessing parameters ---
uncertaintyX     = 20;    % nm (localization precision in x)
uncertaintyY     = 20;    % nm
uncertaintyZ     = 50;    % nm
uncertaintyTheta = 0.1;   % rad
uncertaintyPhi   = 0.2;   % rad
uncertaintyOmega = 0.3;   % sr

% --- Stage 2: Frame-to-frame linking parameters ---
linkParam = struct();
linkParam.linearMotion   = 0;     % 0 = Brownian, 1 = directed (no reversal), 2 = directed (with reversal)
linkParam.minSearchRadius = 50;   % nm — minimum spatial search radius
linkParam.maxSearchRadius = 500;  % nm — maximum spatial search radius
linkParam.brownStdMult   = 3;     % Multiplier for Brownian search radius
linkParam.useLocalDensity = 1;    % Use local density for adaptive search radius
linkParam.nnWindow       = 4;     % Frames for nearest-neighbor distance
linkParam.diagnostics    = [];    % Frame numbers for diagnostic output (e.g., [1 50 100])
% 6D parameters
linkParam.wSpatial       = 1.0;   % Weight for spatial cost
linkParam.wOrient        = 0.3;   % Weight for orientation cost (tune this!)
linkParam.wOmega         = 0.5;   % Weight for wobble cost
linkParam.maxAngularDist = pi/4;  % Max angular distance for linking (rad)
linkParam.useOrientation = true;  % Set false to revert to position-only linking

% --- Stage 3: Gap closing parameters ---
gapParam = struct();
gapParam.brownStdMult    = 3;
gapParam.linStdMult      = 3;
gapParam.minSearchRadius = 50;    % nm
gapParam.maxSearchRadius = 1000;  % nm
gapParam.timeReachConfB  = 3;     % Confidence time for Brownian (frames)
gapParam.timeReachConfL  = 4;
gapParam.brownScaling    = [0.5 0.01];
gapParam.linScaling      = [0.5 0.01];
gapParam.gapPenalty      = 1.5;   % Penalty per frame of gap
gapParam.resLimit        = 30;    % Resolution limit (nm)
% 6D parameters
gapParam.wSpatial        = 1.0;
gapParam.wOrient         = 0.3;
gapParam.wOmega          = 0.5;
gapParam.maxAngularDist  = pi/3;  % More permissive than linking (gaps allow more change)
gapParam.useOrientation  = true;
gapParam.useBoScaling    = true;  % Bo Shuang's sqrt(timeGap) normalization

% General gap closing parameters (u-track interface)
gapCloseParam.timeWindow  = 5;    % Maximum gap size in frames (tune to blinking stats!)
gapCloseParam.mergeSplit  = 0;    % 0 = no merge/split, 1 = both, 2 = merge only, 3 = split only
gapCloseParam.minTrackLen = 3;    % Minimum track segment length for gap closing

% --- Stage 4: Post-processing parameters ---
connectsize = 3;    % Max gap to interpolate (frames)
minTrajLen  = 5;    % Minimum trajectory length to keep (frames)

% --- Stage 5: Analysis parameters ---
minDisplacements = 3;  % Minimum displacement steps for diffusion estimation

% --- Pipeline control ---
runStage = [1 1 1 1 1];  % Set to 0 to skip a stage: [S1 S2&3 S2&3 S4 S5]

%% ========================================================================
%  VALIDATE CONFIGURATION
%  ========================================================================

if isempty(inputFile) || isempty(outputDir)
    error('run6DTracking:noConfig', ...
        'Please set inputFile and outputDir in the USER CONFIGURATION section.');
end

if ~exist(inputFile, 'file')
    error('run6DTracking:fileNotFound', 'Input file not found: %s', inputFile);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% ========================================================================
%  STAGE 1: PREPROCESSING
%  ========================================================================

if runStage(1)
    fprintf('\n--- STAGE 1: Preprocessing ---\n');
    tic;
    
    movieInfo = convert6DSMOLMtoMovieInfo(inputFile, ...
        'ColumnMap', colMap, ...
        'AngleUnits', 'degrees', ...
        'WobbleInput', 'gamma', ...
        'UncertaintyX', uncertaintyX, ...
        'UncertaintyY', uncertaintyY, ...
        'UncertaintyZ', uncertaintyZ, ...
        'UncertaintyTheta', uncertaintyTheta, ...
        'UncertaintyPhi', uncertaintyPhi, ...
        'UncertaintyOmega', uncertaintyOmega, ...
        'Verbose', true);
    
    % Save intermediate result
    save(fullfile(outputDir, 'movieInfo_6D.mat'), 'movieInfo', '-v7.3');
    fprintf('  Saved: movieInfo_6D.mat (%.1f s)\n', toc);
    
else
    fprintf('\n--- STAGE 1: Skipped (loading previous result) ---\n');
    load(fullfile(outputDir, 'movieInfo_6D.mat'), 'movieInfo');
end

%% ========================================================================
%  STAGES 2 & 3: TRACKING (Frame linking + Gap closing)
%  ========================================================================

if runStage(2)
    fprintf('\n--- STAGES 2 & 3: Tracking ---\n');
    tic;
    
    % --- Set up cost matrices for u-track ---
    % Cost Matrix 1: frame-to-frame linking
    costMatrices(1).funcName = 'costMat6DSMOLMLink';
    costMatrices(1).parameters = linkParam;
    
    % Cost Matrix 2: gap closing
    costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';
    costMatrices(2).parameters = gapParam;
    
    % --- Kalman filter functions ---
    % Use u-track's built-in linear motion Kalman filter for spatial coordinates.
    % Our cost functions handle orientation separately.
    kalmanFunctions.reserveMem    = 'kalmanResMemLM';
    kalmanFunctions.initialize    = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain      = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse   = 'kalmanReverseLinearMotion';
    
    % --- Problem dimensionality ---
    probDim = 3;  % Spatial dimensions for Kalman filter (x, y, z)
    
    % --- Run tracking ---
    fprintf('  Running trackCloseGapsKalmanSparse...\n');
    
    [tracksFinal, kalmanInfoLink, errFlag] = trackCloseGapsKalmanSparse(...
        movieInfo, costMatrices, gapCloseParam, kalmanFunctions, ...
        probDim, 0, 1);
    
    if errFlag
        warning('run6DTracking:trackingError', 'Tracking returned error flag.');
    end
    
    fprintf('  Tracking complete: %d tracks found (%.1f s)\n', length(tracksFinal), toc);
    
    % Save tracking results
    save(fullfile(outputDir, 'tracksFinal_6D.mat'), 'tracksFinal', 'kalmanInfoLink', '-v7.3');
    fprintf('  Saved: tracksFinal_6D.mat\n');
    
else
    fprintf('\n--- STAGES 2 & 3: Skipped (loading previous result) ---\n');
    load(fullfile(outputDir, 'tracksFinal_6D.mat'), 'tracksFinal');
end

%% ========================================================================
%  STAGE 4: POST-PROCESSING
%  ========================================================================

if runStage(4)
    fprintf('\n--- STAGE 4: Post-processing ---\n');
    tic;
    
    % Extract 6D trajectory data from u-track output
    trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', true);
    
    % Filter and interpolate
    [trajFiltered, filtStats] = traj_filt_6D(trajRaw, ...
        'ConnectSize', connectsize, ...
        'MinTrajLen', minTrajLen, ...
        'Verbose', true);
    
    % Save
    save(fullfile(outputDir, 'trajectories_6D_raw.mat'), 'trajRaw');
    save(fullfile(outputDir, 'trajectories_6D_filtered.mat'), 'trajFiltered', 'filtStats');
    
    % Also save as CSV for easy inspection
    if ~isempty(trajFiltered)
        header = 'frame,trackID,x_nm,y_nm,z_nm,theta,phi,omega,intensity';
        fid = fopen(fullfile(outputDir, 'trajectories_6D_filtered.csv'), 'w');
        fprintf(fid, '%s\n', header);
        fclose(fid);
        dlmwrite(fullfile(outputDir, 'trajectories_6D_filtered.csv'), ...
            trajFiltered, '-append', 'delimiter', ',', 'precision', '%.6f');
    end
    
    fprintf('  Post-processing complete (%.1f s)\n', toc);
    
else
    fprintf('\n--- STAGE 4: Skipped (loading previous result) ---\n');
    load(fullfile(outputDir, 'trajectories_6D_filtered.mat'), 'trajFiltered');
end

%% ========================================================================
%  STAGE 5: DIFFUSION ANALYSIS
%  ========================================================================

if runStage(5)
    fprintf('\n--- STAGE 5: Diffusion Analysis ---\n');
    tic;
    
    results = estimateDiffusion6D_MLE(trajFiltered, ...
        'FrameTime', frameTime, ...
        'LocPrecisionXY', uncertaintyX, ...
        'LocPrecisionZ', uncertaintyZ, ...
        'LocPrecisionTheta', uncertaintyTheta, ...
        'LocPrecisionPhi', uncertaintyPhi, ...
        'MinDisplacements', minDisplacements, ...
        'Verbose', true);
    
    % Save
    save(fullfile(outputDir, 'diffusion_results_6D.mat'), 'results');
    
    % Save summary as CSV
    if ~isempty(results)
        T = struct2table(results);
        writetable(T, fullfile(outputDir, 'diffusion_results_6D.csv'));
    end
    
    fprintf('  Analysis complete (%.1f s)\n', toc);
    
else
    fprintf('\n--- STAGE 5: Skipped ---\n');
end

%% ========================================================================
%  DONE
%  ========================================================================

fprintf('\n=== Pipeline complete ===\n');
fprintf('Output directory: %s\n', outputDir);
fprintf('Finished: %s\n', datestr(now));
