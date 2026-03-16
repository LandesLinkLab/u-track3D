%% run_demo.m
%  One-stop demo for the 6D-SMOLM tracking pipeline.
%
%  Generates simulated trajectories, runs all pipeline stages, and
%  compares tracking results against ground truth.
%
%  USAGE:
%    1. cd to the software6D folder (or add it to path)
%    2. Run this script
%
%  If u-track3D is on the MATLAB path, the full pipeline runs (Stages 1-5).
%  If u-track3D is NOT available, Stages 1, 4, and 5 run using the ground-
%  truth trajectories directly (skipping the LAP tracking in Stages 2&3).
%  This lets you test everything except the actual linking.
%
%  =========================================================================
%  COMPREHENSIVE PARAMETER REFERENCE
%  =========================================================================
%
%  --- TRACKING WEIGHT PARAMETERS (linkWeights / gapWeights) ---
%  wSpatial        : Weight for spatial (Kalman-predicted) distance cost
%                    Default: 1.0. Set to 0 for orientation-only tracking.
%  wOrient         : Weight for orientation angle cost
%                    Default: 0.3 (linking), 0.2 (gap closing)
%  wOmega          : Weight for wobble (Omega) parameter cost
%                    Default: 0.5 (linking), 0.1 (gap closing)
%
%  --- ORIENTATION KALMAN PARAMETERS (kalmanParams) ---
%  orientKalmanLevel : Controls sophistication of orientation tracking
%                      0 = No orientation cost (spatial only, like u-track3D)
%                      1 = Static variance (Δα²/maxAngDist², original method)
%                      2 = Adaptive pseudo-Kalman (Brownian: Δα²/S, S=4D·dt+2σ²)
%                      3 = Full 3-model Kalman (forward/backward/Brownian rotation)
%                      Default: 2
%  D_rot_prior     : Prior estimate of rotational diffusion coefficient [rad²/s]
%                    Used in Levels 2 & 3 for expected angular variance.
%                    Default: 0.05 rad²/s
%  orientCRLB      : Orientation measurement uncertainty (Cramér-Rao Lower Bound) [rad]
%                    Default: 5° (5 * pi/180 rad)
%  frameTime       : Time between frames [seconds]
%                    Default: 0.05 s (50 ms, 20 fps)
%  useSignalWeighting : Scale orientCRLB by 1/sqrt(signal) for each molecule
%                       Default: false
%
%  --- SIMULATION PARAMETERS (simParams for generate6DTrajectories) ---
%  NumParticles    : Number of molecules to simulate. Default: 5
%  NumFrames       : Total frames in simulation. Default: 300
%  FrameTime       : Time per frame [seconds]. Default: 0.05
%  FieldOfView     : [xmin xmax ymin ymax] in nm. Default: [0 6000 0 6000]
%  ZRange          : [zmin zmax] in nm. Default: [-500 500]
%
%  --- MOTION PARAMETERS ---
%  Dtrans          : Translational diffusion coefficient [nm²/s]. Default: 500
%  Drot            : Rotational diffusion coefficient [rad²/s]. Default: 0.05
%  SpatialModel    : 'brownian', 'confined', 'directed', 'hop', 'anomalous'. Default: 'brownian'
%  RotationalModel : 'brownian', 'confined', 'fixed'. Default: 'brownian'
%                    Rotational diffusion uses the geometric integration scheme of
%                    Höfling & Straube, Phys. Rev. Research 7, 043034 (2025).
%                    See generate6DTrajectories.m for implementation details.
%  ConfinementRadius : For confined model [nm]. Default: 200
%  DriftVelocity   : For directed model [vx, vy, vz] in nm/s. Default: [0 0 0]
%  HopRate         : For hop model, probability per frame. Default: 0.005
%  HopDistance     : For hop model [nm]. Default: 400
%
%  --- CHATTERJEE ET AL. MOTION PARAMETERS (Chem. Biomed. Imaging 2025) ---
%  These parameters are used by the auxiliary functions in tests/:
%    AnomalousDiffusion.m, ConfinedDiffusion.m, DirectedMotion.m, NormalDiffusion.m
%
%  AD_amin         : Minimum anomalous exponent α (sub-diffusion: 0 < α < 1). Default: 0.4
%  AD_amax         : Maximum anomalous exponent α. Default: 0.4
%                    α < 1 = sub-diffusion, α = 1 = normal, α > 1 = super-diffusion
%  CD_Bmin         : Minimum confinement B parameter (dimensionless). Default: 3
%  CD_Bmax         : Maximum confinement B parameter. Default: 3
%                    Larger B = stronger confinement. r = √(D·T·dt) / B^(1/3)
%  DM_kdm1         : Minimum speed coefficient k (speed = k·√D). Default: 10
%  DM_kdm2         : Maximum speed coefficient k. Default: 10
%                    k ≈ 1-5 = weak directed, k ≈ 10 = balanced, k > 20 = strong directed
%
%  --- ORIENTATION PARAMETERS ---
%  GammaMean       : Mean wobble parameter γ ∈ [0,1]. Default: 0.7
%                    γ=1 is fixed dipole, γ=0 is isotropic
%  GammaStd        : Std dev of gamma across particles. Default: 0.1
%  ThetaRange      : [min max] polar angle range [deg]. Default: [0 180]
%  PhiRange        : [min max] azimuthal angle range [deg]. Default: [0 360]
%
%  --- PHOTOPHYSICS PARAMETERS ---
%  Signal          : Mean photon count per detection. Default: 1500
%  SignalStd       : Std dev of signal. Default: 300
%  BlinkOffRate    : Probability of turning OFF per frame. Default: 0.10
%                    Related to τ_on: BlinkOffRate ≈ FrameTime / τ_on
%  BlinkOnRate     : Probability of turning ON per frame. Default: 0.70
%                    Related to τ_off: BlinkOnRate ≈ FrameTime / τ_off
%  BleachRate      : Probability of permanent bleaching per frame. Default: 0.0005
%
%  --- LOCALIZATION NOISE PARAMETERS ---
%  AddNoise        : Enable localization noise. Default: true
%  NoiseXY         : Lateral localization precision [nm]. Default: 20
%  NoiseZ          : Axial localization precision [nm]. Default: 50
%  NoiseTheta      : Polar angle precision [degrees]. Default: 5
%  NoisePhi        : Azimuthal angle precision [degrees]. Default: 5
%  NoiseGamma      : Wobble parameter noise. Default: 0.05
%
%  --- u-track3D COST MATRIX PARAMETERS ---
%  linearMotion    : 0=Brownian, 1=directed, 2=switching. Default: 0
%  minSearchRadius : Minimum linking distance [nm]. Default: 50
%  maxSearchRadius : Maximum linking distance [nm]. Default: 500
%  brownStdMult    : Multiplier for Brownian search radius. Default: 3
%  useLocalDensity : Adapt search radius to local density. Default: 1 (true)
%  nnWindow        : Frames for nearest-neighbor density. Default: 4
%  maxAngularDist  : Max angular distance for Level 1 normalization [rad]. Default: pi/2
%
%  --- GAP CLOSING PARAMETERS ---
%  timeWindow      : Max gap size to close [frames]. Default: 5
%  mergeSplit      : Allow merge/split events. Default: 0 (false)
%  minTrackLen     : Minimum track length to keep. Default: 2
%  gapPenalty      : Cost penalty per frame of gap. Default: 1.5
%
%  --- DIAGNOSTIC PARAMETERS ---
%  saveCostMatrix  : Save cost matrix for visualization. Default: true
%  costMatSavePath : Path to save cost matrix. Default: outputDir
%  costMatSaveFrame: Which frame to save. Default: 5
%
%  =========================================================================
%  BENCHMARKING QUICK REFERENCE
%  =========================================================================
%  
%  SNR Sweep (localization noise):
%    Low:  'Signal', 500,  'NoiseXY', 50, 'NoiseTheta', 15
%    Med:  'Signal', 1500, 'NoiseXY', 20, 'NoiseTheta', 5
%    High: 'Signal', 5000, 'NoiseXY', 10, 'NoiseTheta', 2
%
%  Density Sweep:
%    Low:  'NumParticles', 3
%    Med:  'NumParticles', 10
%    High: 'NumParticles', 30
%
%  Blinking Sweep:
%    Light:  'BlinkOffRate', 0.02, 'BlinkOnRate', 0.90
%    Medium: 'BlinkOffRate', 0.10, 'BlinkOnRate', 0.70
%    Heavy:  'BlinkOffRate', 0.30, 'BlinkOnRate', 0.40
%
%  Motion Speed Sweep:
%    Slow: 'Dtrans', 100,  'Drot', 0.01
%    Med:  'Dtrans', 500,  'Drot', 0.05
%    Fast: 'Dtrans', 2000, 'Drot', 0.20
%
%  =========================================================================
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

clear; clc;
fprintf('============================================\n');
fprintf('  6D-SMOLM Tracking Pipeline — Demo\n');
fprintf('============================================\n\n');

%% ========================================================================
%  SETUP
%  ========================================================================

% Add software6D and subfolders to path
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
software6DRoot = fullfile(thisDir, '..');
addpath(genpath(software6DRoot));

% Add u-track3D to path (adjust to your installation if needed)
utrack3DPath = 'C:\Users\emilg\OneDrive - University of Illinois - Urbana\2026\u-track3D';
if exist(utrack3DPath, 'dir')
    addpath(genpath(fullfile(utrack3DPath, 'software')));
end

% Check if u-track3D is available
hasUtrack = exist('trackCloseGapsKalmanSparse', 'file') == 2;
if hasUtrack
    fprintf('  u-track3D: FOUND on path\n');
    fprintf('  -> Full pipeline will run (Stages 1-5)\n\n');
else
    fprintf('  u-track3D: NOT FOUND on path\n');
    fprintf('  -> Will skip Stages 2&3 (tracking) and use ground truth instead\n');
    fprintf('  -> To enable full tracking, run:\n');
    fprintf('       addpath(genpath(''path/to/u-track3D/software''))\n\n');
end

% Output directory
outputDir = fullfile(thisDir, 'demo_output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% ========================================================================
%  TRACKING WEIGHT PARAMETERS (Easy to adjust!)
%  ========================================================================
%  Modify these to test different weighting schemes:
%    - Set wSpatial=0 to use ONLY orientation for linking
%    - Increase wOrient to emphasize orientation matching
%    - Set wOmega>0 to also penalize wobble differences

% Frame-to-frame linking weights
linkWeights.wSpatial = 1.0;      % Weight for spatial (Kalman-predicted) cost
linkWeights.wOrient  = 1.0;      % Weight for orientation cost
linkWeights.wOmega   = 0.0;      % Weight for wobble (Omega) cost

% Gap closing weights
gapWeights.wSpatial = 1.0;       % Weight for spatial cost
gapWeights.wOrient  = 1.0;       % Weight for orientation cost  
gapWeights.wOmega   = 0.0;       % Weight for wobble cost

%% ========================================================================
%  ORIENTATION KALMAN PARAMETERS (NEW!)
%  ========================================================================
%  Choose orientation Kalman level:
%    Level 0: No orientation cost (spatial only, equivalent to u-track3D)
%    Level 1: Static variance (simple Δα²/maxAngDist², original method)
%    Level 2: Adaptive pseudo-Kalman (Brownian, Δα²/S where S=4D·dt+2σ²)
%    Level 3: Full 3-model Kalman (forward/backward/Brownian rotation)
%
%  Parameters for Levels 2 & 3:
%    - D_rot_prior: expected rotational diffusion coefficient (rad²/s)
%    - orientCRLB: measurement uncertainty for orientation (radians)
%    - useSignalWeighting: scale uncertainty by 1/sqrt(signal)

kalmanParams.orientKalmanLevel = 3;        % 0=none, 1=static, 2=adaptive, 3=full
kalmanParams.D_rot_prior = 0.05;           % Prior D_rot estimate (rad²/s)
kalmanParams.orientCRLB = 5 * pi/180;      % Orientation CRLB (~5 degrees)
kalmanParams.frameTime = 0.05;             % Frame time (seconds)
kalmanParams.useSignalWeighting = false;   % Weight by signal strength

% Display current settings
fprintf('=== Tracking Weight Settings ===\n');
fprintf('  Frame-to-frame: wSpatial=%.1f, wOrient=%.2f, wOmega=%.2f\n', ...
    linkWeights.wSpatial, linkWeights.wOrient, linkWeights.wOmega);
fprintf('  Gap closing:    wSpatial=%.1f, wOrient=%.2f, wOmega=%.2f\n', ...
    gapWeights.wSpatial, gapWeights.wOrient, gapWeights.wOmega);
if linkWeights.wSpatial == 0
    fprintf('  *** ORIENTATION-ONLY MODE: Spatial cost disabled! ***\n');
    fprintf('      -> Search radius will be expanded to avoid spatial gating\n');
end
levelNames = {'None (spatial only)', 'Static', 'Adaptive (Brownian)', 'Full 3-Model'};
fprintf('  *** ORIENTATION KALMAN LEVEL: %d (%s) ***\n', ...
    kalmanParams.orientKalmanLevel, levelNames{kalmanParams.orientKalmanLevel + 1});
if kalmanParams.orientKalmanLevel >= 2
    fprintf('      D_rot_prior=%.3f rad²/s, orientCRLB=%.1f°, signalWeight=%d\n', ...
        kalmanParams.D_rot_prior, kalmanParams.orientCRLB * 180/pi, kalmanParams.useSignalWeighting);
end
fprintf('\n');

% Search radius settings - use scenario-defined values if available, else defaults
% (scenarioSearchRadius is defined within each scenario block above)
if exist('scenarioSearchRadius', 'var') && isstruct(scenarioSearchRadius)
    searchRadiusSettings = scenarioSearchRadius;
    fprintf('  Using scenario-defined search radius: min=%d, max=%d, brownStdMult=%d\n', ...
        searchRadiusSettings.minSearchRadius, searchRadiusSettings.maxSearchRadius, ...
        searchRadiusSettings.brownStdMult);
else
    % Default search radius settings
    searchRadiusSettings.minSearchRadius = 50;
    searchRadiusSettings.maxSearchRadius = 500;
    searchRadiusSettings.brownStdMult = 3;
end

% If orientation-only mode (wSpatial=0), expand search radius to avoid spatial gating
% This overrides scenario settings to ensure orientation-only tracking works
if linkWeights.wSpatial == 0
    searchRadiusSettings.minSearchRadius = 10;
    searchRadiusSettings.maxSearchRadius = 10000;  % Very large - effectively no spatial gating
    searchRadiusSettings.brownStdMult = 20;        % Large multiplier
    fprintf('  *** Orientation-only mode: search radius expanded ***\n');
    fprintf('      min=%d, max=%d, brownStdMult=%d\n', ...
        searchRadiusSettings.minSearchRadius, searchRadiusSettings.maxSearchRadius, ...
        searchRadiusSettings.brownStdMult);
end

%% ========================================================================
%  CHOOSE SIMULATION SCENARIO
%  ========================================================================
%  Uncomment ONE scenario block, or define your own.
%
%  Motion model parameters from Chatterjee et al., Chem. Biomed. Imaging (2025):
%    - AD (Anomalous):  amin, amax  (anomalous exponent α, typically 0.4)
%    - CD (Confined):   Bmin, Bmax  (dimensionless confinement, typically 3)
%    - DM (Directed):   kdm1, kdm2  (speed coefficient k, speed = k·√D, typically 10)
%    - ND (Normal):     Standard Brownian diffusion (no extra params)

% --- Scenario A: Simple Brownian (default) ---
scenarioName = 'brownian';
simParams = {
    'NumParticles', 5, ...
    'NumFrames', 300, ...
    'Dtrans', 500, ...
    'Drot', 0.05, ...
    'SpatialModel', 'brownian', ...
    'RotationalModel', 'brownian', ...
    'GammaMean', 0.7, ...
    'Signal', 1500, ...
    'BlinkOffRate', 0.0, ...
    'BlinkOnRate', 1.0, ...
    'StartSpread', 500, ...
    'Seed', 42 ...
};
% Search radius settings for Scenario A (defaults)
scenarioSearchRadius.minSearchRadius = 50;
scenarioSearchRadius.maxSearchRadius = 500;
scenarioSearchRadius.brownStdMult = 3;

% % --- Scenario B: Confined diffusion (Chatterjee B-parameter method) ---
% scenarioName = 'confined';
% % Confinement parameters from Chatterjee et al.:
% CD_Bmin = 3;    % Minimum B parameter (dimensionless confinement)
% CD_Bmax = 3;    % Maximum B parameter (larger = stronger confinement)
% simParams = {
%     'NumParticles', 4, ...
%     'NumFrames', 500, ...
%     'Dtrans', 300, ...
%     'Drot', 0.02, ...
%     'SpatialModel', 'confined', ...
%     'ConfinementRadius', 200, ...  % Fallback for standard method
%     'CD_Bmin', CD_Bmin, ...        % Chatterjee B-parameter (if using aux function)
%     'CD_Bmax', CD_Bmax, ...
%     'RotationalModel', 'confined', ...
%     'GammaMean', 0.85, ...
%     'Signal', 2000, ...
%     'Seed', 123 ...
% };
% % Search radius settings for Scenario B
% scenarioSearchRadius.minSearchRadius = 50;
% scenarioSearchRadius.maxSearchRadius = 400;  % Smaller for confined motion
% scenarioSearchRadius.brownStdMult = 3;

% % --- Scenario C: Directed transport (Chatterjee fluctuating direction) ---
% scenarioName = 'directed';
% % Directed motion parameters from Chatterjee et al.:
% DM_kdm1 = 10;   % Minimum speed coefficient (speed = k·√D)
% DM_kdm2 = 10;   % Maximum speed coefficient
% simParams = {
%     'NumParticles', 3, ...
%     'NumFrames', 400, ...
%     'Dtrans', 200, ...
%     'Drot', 0.01, ...
%     'SpatialModel', 'directed', ...
%     'DriftVelocity', [100 50 0], ...  % Fallback for standard method
%     'DM_kdm1', DM_kdm1, ...           % Chatterjee k-parameter (if using aux function)
%     'DM_kdm2', DM_kdm2, ...
%     'RotationalModel', 'fixed', ...
%     'GammaMean', 0.95, ...
%     'Signal', 1800, ...
%     'Seed', 7 ...
% };
% % Search radius settings for Scenario C
% scenarioSearchRadius.minSearchRadius = 50;
% scenarioSearchRadius.maxSearchRadius = 800;  % Larger for directed motion
% scenarioSearchRadius.brownStdMult = 4;

% % --- Scenario D: Hopping between traps ---
% scenarioName = 'hop';
% simParams = {
%     'NumParticles', 3, ...
%     'NumFrames', 800, ...
%     'Dtrans', 400, ...
%     'Drot', 0.08, ...
%     'SpatialModel', 'hop', ...
%     'ConfinementRadius', 100, ...
%     'HopRate', 0.005, ...
%     'HopDistance', 400, ...
%     'RotationalModel', 'brownian', ...
%     'GammaMean', 0.6, ...
%     'Seed', 99 ...
% };
% % Search radius settings for Scenario D
% scenarioSearchRadius.minSearchRadius = 50;
% scenarioSearchRadius.maxSearchRadius = 600;  % Medium for hopping
% scenarioSearchRadius.brownStdMult = 3;

% % --- Scenario E: Anomalous (sub-diffusive) motion ---
% scenarioName = 'anomalous';
% % Anomalous diffusion parameters from Chatterjee et al.:
% AD_amin = 0.4;  % Minimum anomalous exponent (α < 1 = sub-diffusion)
% AD_amax = 0.4;  % Maximum anomalous exponent
% simParams = {
%     'NumParticles', 5, ...
%     'NumFrames', 300, ...
%     'Dtrans', 500, ...
%     'Drot', 0.03, ...
%     'SpatialModel', 'anomalous', ...
%     'AD_amin', AD_amin, ...
%     'AD_amax', AD_amax, ...
%     'RotationalModel', 'brownian', ...
%     'GammaMean', 0.7, ...
%     'Signal', 1500, ...
%     'BlinkOffRate', 0.0, ...
%     'BlinkOnRate', 1.0, ...
%     'StartSpread', 500, ...
%     'Seed', 55 ...
% };
% % Search radius settings for Scenario E
% scenarioSearchRadius.minSearchRadius = 30;
% scenarioSearchRadius.maxSearchRadius = 400;  % Smaller for sub-diffusion
% scenarioSearchRadius.brownStdMult = 3;

% % --- Scenario F: Orientation-only tracking (wSpatial=0) ---
% % Use this to test pure orientation-based linking
% scenarioName = 'orient_only';
% simParams = {
%     'NumParticles', 5, ...
%     'NumFrames', 300, ...
%     'Dtrans', 500, ...
%     'Drot', 0.05, ...
%     'SpatialModel', 'brownian', ...
%     'RotationalModel', 'brownian', ...
%     'GammaMean', 0.7, ...
%     'Signal', 1500, ...
%     'BlinkOffRate', 0.0, ...
%     'BlinkOnRate', 1.0, ...
%     'StartSpread', 500, ...
%     'Seed', 42 ...
% };
% % Search radius settings for Scenario F (expanded for orientation-only)
% scenarioSearchRadius.minSearchRadius = 10;
% scenarioSearchRadius.maxSearchRadius = 10000;  % Very large - no spatial gating
% scenarioSearchRadius.brownStdMult = 20;
% % Also set weights for orientation-only:
% % linkWeights.wSpatial = 0.0;
% % linkWeights.wOrient = 1.0;

fprintf('--- Scenario: %s ---\n\n', scenarioName);

%% ========================================================================
%  STAGE 0: GENERATE SIMULATED DATA
%  ========================================================================
fprintf('=== STAGE 0: Generating simulated trajectories ===\n');

simFile = fullfile(outputDir, sprintf('sim_%s.mat', scenarioName));
[data, groundTruth] = generate6DTrajectories(simParams{:}, ...
    'OutputFile', simFile);

fprintf('\n');

%% ========================================================================
%  STAGE 1: PREPROCESSING
%  ========================================================================
fprintf('=== STAGE 1: Converting to movieInfo ===\n');
tic;

movieInfo = convert6DSMOLMtoMovieInfo(data, ...
    'AngleUnits', 'degrees', ...
    'WobbleInput', 'gamma', ...
    'Verbose', true);

save(fullfile(outputDir, 'movieInfo_6D.mat'), 'movieInfo', '-v7.3');
fprintf('  Stage 1 complete (%.1f s)\n\n', toc);

%% ========================================================================
%  STAGES 2 & 3: TRACKING
%  ========================================================================
if hasUtrack
    fprintf('=== STAGES 2 & 3: Tracking with u-track3D ===\n');
    tic;
    
    % --- Cost Matrix 1: frame-to-frame linking ---
    costMatrices(1).funcName = 'costMat6DSMOLMLink';
    costMatrices(1).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', searchRadiusSettings.minSearchRadius, ...
        'maxSearchRadius', searchRadiusSettings.maxSearchRadius, ...
        'brownStdMult', searchRadiusSettings.brownStdMult, ...
        'useLocalDensity', 1, ...
        'nnWindow', 4, ...
        'diagnostics', [], ...
        'wSpatial', linkWeights.wSpatial, ...       % Use variable from top of script
        'wOrient', linkWeights.wOrient, ...         % Use variable from top of script
        'wOmega', linkWeights.wOmega, ...           % Use variable from top of script
        'maxAngularDist', pi/2, ...
        'useOrientation', true, ...
        'orientKalmanLevel', kalmanParams.orientKalmanLevel, ...  % 0=none, 1=static, 2=adaptive, 3=full
        'D_rot_prior', kalmanParams.D_rot_prior, ...              % Rotational diffusion prior
        'orientCRLB', kalmanParams.orientCRLB, ...                % Orientation measurement uncertainty
        'frameTime', kalmanParams.frameTime, ...                  % Frame time for diffusion calc
        'useSignalWeighting', kalmanParams.useSignalWeighting, ...% Signal-dependent CRLB
        'saveCostMatrix', true, ...         % 6D-SMOLM: save cost matrix for visualization
        'costMatSavePath', outputDir, ...   % 6D-SMOLM: path to save cost matrix
        'costMatSaveFrame', 5);             % 6D-SMOLM: which frame to save (middle of early frames)
    
    % --- Cost Matrix 2: gap closing ---
    % *** KEY: brownStdMult and linStdMult must be COLUMN vectors! ***
    timeWindow = 5;
    costMatrices(2).funcName = 'costMat6DSMOLMCloseGaps';  % Use our 6D wrapper
    costMatrices(2).parameters = struct( ...
        'linearMotion', 0, ...
        'minSearchRadius', searchRadiusSettings.minSearchRadius, ...      
        'maxSearchRadius', searchRadiusSettings.maxSearchRadius, ...     
        'brownStdMult', searchRadiusSettings.brownStdMult * ones(timeWindow, 1), ... % COLUMN vector of length timeWindow
        'linStdMult', searchRadiusSettings.brownStdMult * ones(timeWindow, 1), ...   % COLUMN vector of length timeWindow
        'timeReachConfB', 3, ...
        'timeReachConfL', 3, ...
        'brownScaling', [0.5 0.01], ...
        'linScaling', [0.5 0.01], ...
        'useLocalDensity', 1, ...
        'nnWindow', 4, ...
        'lenForClassify', 5, ...
        'maxAngleVV', 45, ...
        'gapPenalty', 1.5, ...
        'resLimit', 10, ...
        'wSpatial', gapWeights.wSpatial, ...        % Use variable from top of script
        'wOrient', gapWeights.wOrient, ...          % Use variable from top of script
        'wOmega', gapWeights.wOmega, ...            % Use variable from top of script
        'maxAngularDist', pi/2, ...     % 6D-SMOLM: max angular distance
        'useOrientation', true, ...     % 6D-SMOLM: enable orientation matching
        'useBoScaling', true, ...       % 6D-SMOLM: sqrt(timeGap) normalization
        'saveCostMatrix', true, ...     % 6D-SMOLM: save cost matrix for visualization
        'costMatSavePath', outputDir);  % 6D-SMOLM: path to save cost matrix
    
    % --- Gap closing parameters ---
    gapCloseParam.timeWindow = timeWindow;
    gapCloseParam.mergeSplit = 0;
    gapCloseParam.minTrackLen = 2;  % Reduced from 3 - keep shorter segments for gap closing
    
    % --- Kalman filter ---
    kalmanFunctions.reserveMem  = 'kalmanResMemLM';
    kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
    kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
    kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
    
    probDim = 3;
    
    fprintf('  Running trackCloseGapsKalmanSparse...\n');
    [tracksFinal, kalmanInfoLink, errFlag] = trackCloseGapsKalmanSparse( ...
        movieInfo, costMatrices, gapCloseParam, kalmanFunctions, probDim, 0, 1);
    
    fprintf('  Found %d tracks (%.1f s)\n', length(tracksFinal), toc);
    save(fullfile(outputDir, 'tracksFinal_6D.mat'), 'tracksFinal', 'kalmanInfoLink', '-v7.3');
    
    % Extract to trajectory matrix
    trajRaw = extractTracks6D(tracksFinal, movieInfo, 'Verbose', true);
    trackingSource = 'u-track3D';
    
else
    fprintf('=== STAGES 2 & 3: SKIPPED (no u-track3D) ===\n');
    fprintf('  Using ground-truth trajectories instead.\n');
    
    % Build trajectory matrix from ground truth
    % groundTruth.trajectories{i} columns:
    %   [frame, signal, x, y, z, theta_deg, phi_deg, gamma, trackID, isInterp]
    allGT = vertcat(groundTruth.trajectories{:});
    
    % Convert to pipeline format: [frame, trackID, x, y, z, theta_rad, phi_rad, omega, signal]
    trajRaw = zeros(size(allGT, 1), 9);
    trajRaw(:, 1) = allGT(:, 1);                     % frame
    trajRaw(:, 2) = allGT(:, 9);                     % trackID
    trajRaw(:, 3) = allGT(:, 3);                     % x_nm
    trajRaw(:, 4) = allGT(:, 4);                     % y_nm
    trajRaw(:, 5) = allGT(:, 5);                     % z_nm
    trajRaw(:, 6) = allGT(:, 6) * (pi/180);          % theta: deg -> rad
    trajRaw(:, 7) = allGT(:, 7) * (pi/180);          % phi: deg -> rad
    % gamma -> omega conversion
    gammaVals = max(0, min(1, allGT(:, 8)));
    trajRaw(:, 8) = pi * (3 - sqrt(1 + 8 * gammaVals));  % omega (sr)
    trajRaw(:, 9) = allGT(:, 2);                     % signal
    
    trackingSource = 'ground truth';
end

save(fullfile(outputDir, 'trajectories_6D_raw.mat'), 'trajRaw');
fprintf('  Raw trajectories: %d points, %d tracks (source: %s)\n\n', ...
    size(trajRaw,1), length(unique(trajRaw(:,2))), trackingSource);

%% ========================================================================
%  STAGE 4: POST-PROCESSING
%  ========================================================================
fprintf('=== STAGE 4: Filtering trajectories ===\n');
tic;

[trajFiltered, filtStats] = traj_filt_6D(trajRaw, ...
    'ConnectSize', 10, ...   % Increased from 3 - allow larger gaps
    'MinTrajLen', 5, ...
    'Verbose', true);

save(fullfile(outputDir, 'trajectories_6D_filtered.mat'), 'trajFiltered', 'filtStats');
fprintf('  Stage 4 complete (%.1f s)\n\n', toc);

%% ========================================================================
%  STAGE 5: DIFFUSION ANALYSIS (Dual MLE)
%  ========================================================================
fprintf('=== STAGE 5: Diffusion analysis (Dual MLE) ===\n');
tic;

frameTime = 0.05;  % Must match simulation FrameTime

% Get noise parameters from simulation (if AddNoise was true)
% These MUST match what was used in generate6DTrajectories
if isfield(groundTruth.params, 'AddNoise') && groundTruth.params.AddNoise
    locPrecXY = groundTruth.params.NoiseXY;
    locPrecZ = groundTruth.params.NoiseZ;
    locPrecTheta = groundTruth.params.NoiseTheta * pi/180;  % Convert deg to rad
    locPrecPhi = groundTruth.params.NoisePhi * pi/180;      % Convert deg to rad
    fprintf('  Using noise parameters from simulation: XY=%.0f nm, Z=%.0f nm, theta=%.1f°, phi=%.1f°\n', ...
        locPrecXY, locPrecZ, groundTruth.params.NoiseTheta, groundTruth.params.NoisePhi);
else
    locPrecXY = 0;
    locPrecZ = 0;
    locPrecTheta = 0;
    locPrecPhi = 0;
    fprintf('  No localization noise in simulation\n');
end

% === OPTION: Use ground truth trajectories for MLE validation ===
% Set to true to bypass tracking errors and validate MLE algorithms directly
useGroundTruthForMLE = true;

% === OPTION: Show MLE(2) results in plots ===
% Set to false to only show True vs MLE(1) comparison (MLE(2) is still broken)
showMLE2inPlots = false;

% === OPTION: Save figures to disk ===
% Set to false to skip saving figures (faster execution)
saveFigures = false;

if useGroundTruthForMLE
    fprintf('  *** Using GROUND TRUTH trajectories for MLE (bypassing tracking) ***\n');
    
    % Build trajectory matrix from ground truth with NOISY positions
    % 'data' contains the noisy observations from generate6DTrajectories
    % data columns: [frame, signal, x, y, z, theta_deg, phi_deg, gamma]
    trajForMLE = [];
    for iPart = 1:length(groundTruth.trajectories)
        gt = groundTruth.trajectories{iPart};
        if isempty(gt), continue; end
        
        % gt columns: [frame, signal, x, y, z, theta_deg, phi_deg, gamma, trackID, isInterp]
        nPts = size(gt, 1);
        
        % Use the OBSERVED (noisy) positions from 'data' variable
        trajPart = zeros(nPts, 9);
        trajPart(:, 1) = gt(:, 1);              % frame
        trajPart(:, 2) = iPart;                 % trackID = particle ID
        
        % Get noisy observations for this particle from 'data'
        for iPt = 1:nPts
            f = gt(iPt, 1);
            % Find matching frame in data (there may be multiple detections per frame)
            % Match by frame AND by proximity to ground truth position
            dataFrame = data(data(:,1) == f, :);
            if ~isempty(dataFrame)
                % Find closest detection to this GT position
                gtPos = [gt(iPt, 3), gt(iPt, 4), gt(iPt, 5)];
                dists = sqrt(sum((dataFrame(:, 3:5) - gtPos).^2, 2));
                [~, closest] = min(dists);
                trajPart(iPt, 3:5) = dataFrame(closest, 3:5);  % x, y, z (noisy)
                trajPart(iPt, 6) = dataFrame(closest, 6) * pi/180;  % theta (rad)
                trajPart(iPt, 7) = dataFrame(closest, 7) * pi/180;  % phi (rad)
                % gamma -> omega
                gamma = max(0, min(1, dataFrame(closest, 8)));
                trajPart(iPt, 8) = pi * (3 - sqrt(1 + 8 * gamma));
                trajPart(iPt, 9) = dataFrame(closest, 2);  % signal
            else
                % Fallback to ground truth if no match (shouldn't happen)
                trajPart(iPt, 3:5) = gt(iPt, 3:5);
                trajPart(iPt, 6) = gt(iPt, 6) * pi/180;
                trajPart(iPt, 7) = gt(iPt, 7) * pi/180;
                gamma = max(0, min(1, gt(iPt, 8)));
                trajPart(iPt, 8) = pi * (3 - sqrt(1 + 8 * gamma));
                trajPart(iPt, 9) = gt(iPt, 2);
            end
        end
        trajForMLE = [trajForMLE; trajPart];
    end
    fprintf('  Ground truth trajectories: %d points, %d particles\n', ...
        size(trajForMLE, 1), length(unique(trajForMLE(:,2))));
else
    % Use tracked trajectories (may have tracking errors)
    trajForMLE = trajFiltered;
    fprintf('  Using TRACKED trajectories for MLE\n');
end

% Run Dual MLE analysis (both MLE(1) and MLE(2) from Bo Shuang 2013)
results = estimateDiffusion6D_DualMLE(trajForMLE, ...
    'FrameTime', frameTime, ...
    'LocPrecisionXY', locPrecXY, ...
    'LocPrecisionZ', locPrecZ, ...
    'LocPrecisionTheta', locPrecTheta, ...
    'LocPrecisionPhi', locPrecPhi, ...
    'MinDisplacements', 5, ...
    'Verbose', true);

save(fullfile(outputDir, 'diffusion_results_6D.mat'), 'results');
fprintf('  Stage 5 complete (%.1f s)\n\n', toc);

%% ========================================================================
%  VALIDATION: COMPARE TO GROUND TRUTH
%  ========================================================================
fprintf('=== VALIDATION ===\n');

trueDtrans = groundTruth.Dtrans;
trueDrot   = groundTruth.Drot;

% Initialize for plotting
nCompare = 0;
estDtrans_MLE1 = [];
estDtrans_MLE2 = [];
estDrot_MLE1 = [];
estDrot_MLE2 = [];

if ~isempty(results)
    estDtrans_MLE1 = [results.D_trans_MLE1]';
    estDtrans_MLE2 = [results.D_trans_MLE2]';
    estDrot_MLE1 = [results.D_rot_MLE1]';
    estDrot_MLE2 = [results.D_rot_MLE2]';
    recommended = [results.recommended_MLE]';
    
    fprintf('\n  %-10s  %10s  %12s  %12s\n', 'Particle', 'True D_t', 'MLE(1) D_t', 'MLE(2) D_t');
    fprintf('  %-10s  %10s  %12s  %12s\n', '--------', '--------', '----------', '----------');
    
    nCompare = min(length(trueDtrans), length(estDtrans_MLE1));
    for i = 1:nCompare
        fprintf('  %-10d  %8.0f    %8.0f      %8.0f\n', ...
            i, trueDtrans(i), estDtrans_MLE1(i), estDtrans_MLE2(i));
    end
    
    fprintf('\n  %-10s  %10s  %12s  %12s\n', 'Particle', 'True D_r', 'MLE(1) D_r', 'MLE(2) D_r');
    fprintf('  %-10s  %10s  %12s  %12s\n', '--------', '--------', '----------', '----------');
    for i = 1:nCompare
        fprintf('  %-10d  %8.4f    %8.4f      %8.4f\n', ...
            i, trueDrot(i), estDrot_MLE1(i), estDrot_MLE2(i));
    end
    
    fprintf('\n  Summary (median ratio est/true):\n');
    fprintf('    D_trans MLE(1): %.2f    D_trans MLE(2): %.2f\n', ...
        median(estDtrans_MLE1(1:nCompare) ./ trueDtrans(1:nCompare)), ...
        median(estDtrans_MLE2(1:nCompare) ./ trueDtrans(1:nCompare)));
    fprintf('    D_rot MLE(1):   %.2f    D_rot MLE(2):   %.2f\n', ...
        median(estDrot_MLE1(1:nCompare) ./ trueDrot(1:nCompare)), ...
        median(estDrot_MLE2(1:nCompare) ./ trueDrot(1:nCompare)));
    fprintf('    Recommended MLE: MLE(%d) for %d%% of trajectories\n', ...
        mode(recommended), round(100*sum(recommended==mode(recommended))/length(recommended)));
else
    fprintf('  No diffusion results to compare.\n');
end

%% ========================================================================
%  FIGURE 0: TRACKING ERROR DIAGNOSIS
%  ========================================================================
fprintf('\n=== Tracking Error Diagnosis ===\n');

% Run error diagnosis comparing tracked trajectories to ground truth
% This generates a figure with:
%   - Error classification (fragmentation, merging, ID swaps)
%   - Spatial displacement histogram
%   - Angular displacement histogram
%   - GT assignment over time

[errorReport, fig0] = diagnoseTrackingErrors(trajFiltered, groundTruth, ...
    'SpatialJumpThreshold', 200, ...   % nm - flag large spatial jumps
    'AngularJumpThreshold', 30, ...    % degrees - flag large orientation jumps
    'MatchRadius', 100, ...            % nm - max distance to match to GT
    'Verbose', true, ...
    'PlotResults', true);

% Save Figure 0
if saveFigures && ~isempty(fig0) && ishandle(fig0)
    saveas(fig0, fullfile(outputDir, sprintf('demo_%s_errorDiagnosis.png', scenarioName)));
    fprintf('  Figure 0 saved: demo_%s_errorDiagnosis.png\n', scenarioName);
    
    % Export individual Figure 0 panels at 600 DPI
    if isfield(errorReport, 'figAxes') && isfield(errorReport, 'figPanelNames')
        fig0PanelDir = fullfile(outputDir, 'Demo_Figure_0_Panels');
        if ~exist(fig0PanelDir, 'dir')
            mkdir(fig0PanelDir);
        end
        fprintf('  Exporting Figure 0 individual panels to: %s\n', fig0PanelDir);
        
        % Local formatting constants for Figure 0 (defined here since global ones come later)
        FMT_EXPORT_DPI_F0 = 600;
        FMT_FONTSIZE_TITLE_F0 = 20;
        FMT_FONTSIZE_LABEL_F0 = 18;
        FMT_FONTSIZE_TICK_F0 = 16;
        FMT_FONTSIZE_LEGEND_F0 = 14;
        
        for i = 1:length(errorReport.figAxes)
            if ~isempty(errorReport.figAxes{i}) && isvalid(errorReport.figAxes{i})
                panelName = errorReport.figPanelNames{i};
                panelAx = errorReport.figAxes{i};
                
                % Create invisible figure for export
                figTemp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 650]);
                
                % Copy the axes to the new figure
                newAx = copyobj(panelAx, figTemp);
                
                % Check if there's a legend associated with this axes
                origLeg = findobj(panelAx.Parent, 'Type', 'Legend');
                hasLegend = false;
                legStrings = {};
                legLocation = 'best';
                for legIdx = 1:length(origLeg)
                    % Check if this legend belongs to our axes
                    try
                        if isequal(origLeg(legIdx).Axes, panelAx)
                            hasLegend = true;
                            legStrings = origLeg(legIdx).String;
                            legLocation = origLeg(legIdx).Location;
                            break;
                        end
                    catch
                        % Legend might not have Axes property in older MATLAB
                    end
                end
                
                if hasLegend && ~isempty(legStrings)
                    newAx.Position = [0.12 0.15 0.70 0.75];
                    % Recreate legend
                    legend(newAx, legStrings, 'Location', legLocation, ...
                        'FontSize', FMT_FONTSIZE_LEGEND_F0, 'FontWeight', 'bold');
                else
                    newAx.Position = [0.15 0.15 0.75 0.75];
                end
                
                % Apply larger font sizes
                newAx.FontSize = FMT_FONTSIZE_TICK_F0;
                newAx.FontWeight = 'bold';
                if ~isempty(newAx.Title)
                    newAx.Title.FontSize = FMT_FONTSIZE_TITLE_F0;
                    newAx.Title.FontWeight = 'bold';
                end
                if ~isempty(newAx.XLabel)
                    newAx.XLabel.FontSize = FMT_FONTSIZE_LABEL_F0;
                    newAx.XLabel.FontWeight = 'bold';
                end
                if ~isempty(newAx.YLabel)
                    newAx.YLabel.FontSize = FMT_FONTSIZE_LABEL_F0;
                    newAx.YLabel.FontWeight = 'bold';
                end
                
                % Export at high DPI
                exportFile = fullfile(fig0PanelDir, sprintf('%s_%s.png', scenarioName, panelName));
                exportgraphics(figTemp, exportFile, 'Resolution', FMT_EXPORT_DPI_F0);
                
                close(figTemp);
            end
        end
        fprintf('  Figure 0 panels exported (%d panels at %d DPI)\n', length(errorReport.figAxes), FMT_EXPORT_DPI_F0);
    end
end

% Save error report
save(fullfile(outputDir, 'errorReport.mat'), 'errorReport');
fprintf('  Error report saved: errorReport.mat\n');

%% ========================================================================
%  QUICK VISUALIZATION
%  ========================================================================
fprintf('\n=== Plotting ===\n');

% === STANDARD FORMATTING CONSTANTS ===
% These are used throughout all figures for consistency
FMT_FONTSIZE_TITLE = 20;
FMT_FONTSIZE_LABEL = 18;
FMT_FONTSIZE_TICK = 16;
FMT_FONTSIZE_LEGEND = 14;
FMT_LINEWIDTH_AXES = 2.0;
FMT_LINEWIDTH_PLOT = 1.5;
FMT_FONTWEIGHT = 'bold';
FMT_EXPORT_DPI = 600;  % DPI for exported panels

% Helper function to apply standard formatting to an axes (includes grid)
applyStandardFormatting = @(ax) set(ax, 'FontSize', FMT_FONTSIZE_TICK, ...
    'FontWeight', FMT_FONTWEIGHT, 'LineWidth', FMT_LINEWIDTH_AXES, ...
    'TickLength', [0.02 0.02], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');

uniqueTracks = unique(trajFiltered(:,2));
nTracks = length(uniqueTracks);

% Get global frame range for consistent colormap scaling
minFrame = min(trajFiltered(:,1));
maxFrame = max(trajFiltered(:,1));

% Helper function to plot trajectory with time-encoded jet colormap
plotTrajWithTime = @(ax, xData, yData, frameData) plotTimeColoredTraj(ax, xData, yData, frameData, minFrame, maxFrame);

%% === FIGURE 1: 3x4 Summary (Raw | Linked | Orientation | Diffusion) ===
fig1 = figure('Position', [50 50 1600 1000], 'Color', 'w', 'Name', ...
    sprintf('6D-SMOLM Demo — %s (Figure 1)', scenarioName));

% Extract raw detection data from movieInfo for comparison
rawX = []; rawY = []; rawZ = []; rawFrames = [];
for f = 1:length(movieInfo)
    nDet = size(movieInfo(f).xCoord, 1);
    if nDet > 0
        rawX = [rawX; movieInfo(f).xCoord(:,1)];
        rawY = [rawY; movieInfo(f).yCoord(:,1)];
        rawZ = [rawZ; movieInfo(f).zCoord(:,1)];
        rawFrames = [rawFrames; repmat(f, nDet, 1)];
    end
end

% =========================================================================
% COLUMN 1: RAW (Unlinked) Detections
% =========================================================================

% --- (1,1): XY raw detections ---
ax11 = subplot(3, 4, 1);
scatter(rawX, rawY, 8, rawFrames, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('XY Detections (Raw)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
axis equal; 
applyStandardFormatting(gca);
cb11 = colorbar; ylabel(cb11, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax11, jet); caxis([minFrame maxFrame]);

% --- (2,1): XZ raw detections ---
ax21 = subplot(3, 4, 5);
scatter(rawX, rawZ, 8, rawFrames, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('z (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('XZ Detections (Raw)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
applyStandardFormatting(gca);
cb21 = colorbar; ylabel(cb21, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax21, jet); caxis([minFrame maxFrame]);

% --- (3,1): 3D raw detections ---
ax31 = subplot(3, 4, 9);
scatter3(rawX, rawY, rawZ, 8, rawFrames, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
zlabel('z (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('3D Detections (Raw)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
axis equal; grid on;
applyStandardFormatting(gca);
view(45, 30);
cb31 = colorbar; ylabel(cb31, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax31, jet); caxis([minFrame maxFrame]);

% =========================================================================
% COLUMN 2: LINKED Trajectories
% =========================================================================

% --- (1,2): XY linked trajectories ---
ax12 = subplot(3, 4, 2);
hold on;
for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    plotTimeColoredTraj(ax12, tData(:,3), tData(:,4), tData(:,1), minFrame, maxFrame);
end
hold off;
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('XY Trajectories (Linked)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
axis equal;
applyStandardFormatting(gca);
cb12 = colorbar; ylabel(cb12, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax12, jet); caxis([minFrame maxFrame]);

% --- (2,2): XZ linked trajectories ---
ax22 = subplot(3, 4, 6);
hold on;
for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    plotTimeColoredTraj(ax22, tData(:,3), tData(:,5), tData(:,1), minFrame, maxFrame);
end
hold off;
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('z (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('XZ Trajectories (Linked)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
applyStandardFormatting(gca);
cb22 = colorbar; ylabel(cb22, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax22, jet); caxis([minFrame maxFrame]);

% --- (3,2): 3D linked trajectories with LINES ---
ax32 = subplot(3, 4, 10);
hold on;
for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    x = tData(:,3); y = tData(:,4); z = tData(:,5);
    frames = tData(:,1);
    
    % Plot line connecting points
    plot3(x, y, z, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    % Overlay colored scatter points
    scatter3(x, y, z, 15, frames, 'filled');
end
hold off;
xlabel('x (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
zlabel('z (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('3D Trajectories (Linked)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
axis equal; grid on;
applyStandardFormatting(gca);
view(45, 30);
cb32 = colorbar; ylabel(cb32, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
colormap(ax32, jet); caxis([minFrame maxFrame]);

% =========================================================================
% COLUMN 3: Orientation vs Time
% =========================================================================

% --- (1,3): Theta (polar angle) vs frame ---
ax13 = subplot(3, 4, 3);
hold on;
for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    plot(tData(:,1), tData(:,6) * 180/pi, '-k', 'LineWidth', 1.2);
end
hold off;
xlabel('Frame', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('\theta (deg)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('\theta (Polar Angle) — Tracked', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
applyStandardFormatting(gca);

% --- (2,3): Phi (azimuthal angle) as POLAR PLOT ---
% Using polar coordinates: angle = phi, radius = frame number
% This avoids the -180/+180 wrap-around discontinuity
ax23 = subplot(3, 4, 7);

% Get position and adjust to avoid clashing with xlabel above
pos23 = get(ax23, 'Position');
pos23(2) = pos23(2) - 0.02;  % Move down slightly
pos23(4) = pos23(4) - 0.02;  % Reduce height slightly

% Convert to polar axes
pax = polaraxes('Position', pos23);
delete(ax23);  % Remove the Cartesian axes
hold(pax, 'on');

for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    
    % Phi angle (already in radians in column 7, range -pi to pi)
    phiVals = tData(:,7);
    
    % Radius = actual frame number (not normalized)
    rVals = tData(:,1);
    
    % Plot as black polar line
    polarplot(pax, phiVals, rVals, '-k', 'LineWidth', FMT_LINEWIDTH_PLOT);
end
hold(pax, 'off');

% Configure polar axes for -180 to 180 range
pax.ThetaZeroLocation = 'right';  % 0° on the right (standard math convention)
pax.ThetaDir = 'counterclockwise';
pax.ThetaLim = [-180 180];
pax.ThetaTick = [-180 -135 -90 -45 0 45 90 135 180];
pax.ThetaTickLabel = {'-180°', '-135°', '-90°', '-45°', '0°', '45°', '90°', '135°', '180°'};
pax.FontSize = FMT_FONTSIZE_TICK;
pax.FontWeight = FMT_FONTWEIGHT;
pax.LineWidth = FMT_LINEWIDTH_AXES;

% Set radial limits and show frame number labels
pax.RLim = [0 maxFrame];

% Title
title(pax, '\phi (Azimuthal) — Tracked', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);

% Add "r = frame" annotation outside plot area (bottom left corner)
annotation('textbox', [pos23(1), pos23(2)-0.025, 0.1, 0.025], ...
    'String', 'r = frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

% --- (3,3): Omega (wobble) vs frame ---
ax33 = subplot(3, 4, 11);
hold on;
for i = 1:nTracks
    idx = trajFiltered(:,2) == uniqueTracks(i);
    tData = trajFiltered(idx, :);
    [~, sortIdx] = sort(tData(:,1));
    tData = tData(sortIdx, :);
    plot(tData(:,1), tData(:,8), '-k', 'LineWidth', FMT_LINEWIDTH_PLOT);
end
hold off;
xlabel('Frame', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
ylabel('\Omega (sr)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
title('\Omega (Wobble) — Tracked', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
applyStandardFormatting(gca);

% =========================================================================
% COLUMN 4: Diffusion estimates and parameters
% =========================================================================

% --- (1,4): D_trans comparison (True vs MLE1, optionally MLE2) ---
ax14 = subplot(3, 4, 4);
if ~isempty(results) && exist('nCompare', 'var') && nCompare > 0 && ~isempty(estDtrans_MLE1)
    trueVals = trueDtrans(1:nCompare);
    est1Vals = estDtrans_MLE1(1:nCompare);
    est2Vals = estDtrans_MLE2(1:nCompare);
    if isrow(trueVals), trueVals = trueVals'; end
    if isrow(est1Vals), est1Vals = est1Vals'; end
    if isrow(est2Vals), est2Vals = est2Vals'; end
    
    % Build bar data based on showMLE2inPlots flag
    if showMLE2inPlots
        barData = [trueVals, est1Vals, est2Vals];
        legendLabels = {'True', 'MLE(1)', 'MLE(2)'};
        barColors = {[0.2 0.4 0.8], [0.9 0.5 0.1], [0.4 0.8 0.4]};
    else
        barData = [trueVals, est1Vals];
        legendLabels = {'True', 'MLE(1)'};
        barColors = {[0.2 0.4 0.8], [0.9 0.5 0.1]};
    end
    
    % Handle single particle case - bar() needs special handling
    if nCompare == 1
        bh = bar(1, barData, 'grouped');
        set(gca, 'XTick', 1, 'XTickLabel', {'P1'});
    else
        bh = bar(barData, 'grouped');
        set(gca, 'XTickLabel', arrayfun(@(i) sprintf('P%d', i), 1:nCompare, 'Uni', false));
    end
    
    % Apply colors
    for iBar = 1:min(numel(bh), numel(barColors))
        bh(iBar).FaceColor = barColors{iBar};
    end
    
    legend(legendLabels, 'Location', 'best', 'FontSize', FMT_FONTSIZE_LEGEND, 'FontWeight', FMT_FONTWEIGHT);
    ylabel('D_{trans} (nm^2/s)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('Translational Diffusion', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    applyStandardFormatting(gca);
else
    text(0.5, 0.5, 'No diffusion data', 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
        'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('Translational Diffusion', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    axis off;
end

% --- (2,4): D_rot comparison (True vs MLE1, optionally MLE2) ---
ax24 = subplot(3, 4, 8);
if ~isempty(results) && exist('nCompare', 'var') && nCompare > 0 && ~isempty(estDrot_MLE1)
    trueVals = trueDrot(1:nCompare);
    est1Vals = estDrot_MLE1(1:nCompare);
    est2Vals = estDrot_MLE2(1:nCompare);
    if isrow(trueVals), trueVals = trueVals'; end
    if isrow(est1Vals), est1Vals = est1Vals'; end
    if isrow(est2Vals), est2Vals = est2Vals'; end
    
    % Build bar data based on showMLE2inPlots flag
    if showMLE2inPlots
        barData = [trueVals, est1Vals, est2Vals];
        legendLabels = {'True', 'MLE(1)', 'MLE(2)'};
        barColors = {[0.2 0.4 0.8], [0.9 0.5 0.1], [0.4 0.8 0.4]};
    else
        barData = [trueVals, est1Vals];
        legendLabels = {'True', 'MLE(1)'};
        barColors = {[0.2 0.4 0.8], [0.9 0.5 0.1]};
    end
    
    % Handle single particle case
    if nCompare == 1
        bh = bar(1, barData, 'grouped');
        set(gca, 'XTick', 1, 'XTickLabel', {'P1'});
    else
        bh = bar(barData, 'grouped');
        set(gca, 'XTickLabel', arrayfun(@(i) sprintf('P%d', i), 1:nCompare, 'Uni', false));
    end
    
    % Apply colors
    for iBar = 1:min(numel(bh), numel(barColors))
        bh(iBar).FaceColor = barColors{iBar};
    end
    
    legend(legendLabels, 'Location', 'best', 'FontSize', FMT_FONTSIZE_LEGEND, 'FontWeight', FMT_FONTWEIGHT);
    ylabel('D_{rot} (rad^2/s)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('Rotational Diffusion', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    applyStandardFormatting(gca);
else
    text(0.5, 0.5, 'No diffusion data', 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
        'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('Rotational Diffusion', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    axis off;
end

% --- (3,4): Simulation parameters text box ---
ax34 = subplot(3, 4, 12);
axis off;

% Parameters as single column (avoid overlap)
paramStr = {
    sprintf('\\bfSimulation\\rm')
    sprintf('  Particles: %d', groundTruth.params.NumParticles)
    sprintf('  Frames: %d', groundTruth.params.NumFrames)
    sprintf('  Frame time: %.0f ms', groundTruth.params.FrameTime * 1000)
    sprintf('  D_{trans}: %.0f nm^2/s', groundTruth.params.Dtrans)
    sprintf('  D_{rot}: %.3f rad^2/s', groundTruth.params.Drot)
    ''
    sprintf('\\bfNoise\\rm')
    sprintf('  Add noise: %s', mat2str(groundTruth.params.AddNoise))
    sprintf('  \\sigma_{xy}: %.1f nm', groundTruth.params.NoiseXY)
    sprintf('  \\sigma_z: %.1f nm', groundTruth.params.NoiseZ)
    sprintf('  \\sigma_\\theta: %.1f°', groundTruth.params.NoiseTheta)
    sprintf('  \\sigma_\\phi: %.1f°', groundTruth.params.NoisePhi)
};

text(0.05, 0.95, paramStr, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'Interpreter', 'tex');
title('Parameters', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);

% Build title with emitter count and weight info
nEmitters = groundTruth.params.NumParticles;
if linkWeights.wSpatial == 0
    weightStr = 'ORIENT-ONLY';
else
    weightStr = sprintf('w_s=%.1f, w_o=%.2f', linkWeights.wSpatial, linkWeights.wOrient);
end
sgtitle(sprintf('6D-SMOLM Demo: %s  |  %d emitters, %d tracks  |  %s  |  source: %s', ...
    scenarioName, nEmitters, nTracks, weightStr, trackingSource), ...
    'FontWeight', FMT_FONTWEIGHT, 'FontSize', FMT_FONTSIZE_TITLE);

% Save Figure 1
if saveFigures
    saveas(fig1, fullfile(outputDir, sprintf('demo_%s_fig1.png', scenarioName)));
    fprintf('  Figure 1 saved: demo_%s_fig1.png\n', scenarioName);
end

%% === Export Individual Figure 1 Panels at 600 DPI ===
if saveFigures
fig1PanelDir = fullfile(outputDir, 'Demo_Figure_1_Panels');
if ~exist(fig1PanelDir, 'dir')
    mkdir(fig1PanelDir);
end
fprintf('  Exporting Figure 1 individual panels to: %s\n', fig1PanelDir);

% Panel names for Figure 1 (3x4 grid)
fig1PanelNames = {
    'XY_Detections_Raw', 'XY_Trajectories_Linked', 'Theta_Tracked', 'Dtrans_Comparison';
    'XZ_Detections_Raw', 'XZ_Trajectories_Linked', 'Phi_Azimuthal_Tracked', 'Drot_Comparison';
    '3D_Detections_Raw', '3D_Trajectories_Linked', 'Omega_Wobble_Tracked', 'Parameters'
};

% Panels that have colorbars (need special handling)
panelsWithColorbar = {'XY_Detections_Raw', 'XZ_Detections_Raw', '3D_Detections_Raw', ...
                      'XY_Trajectories_Linked', 'XZ_Trajectories_Linked', '3D_Trajectories_Linked'};

% Panels that have legends (need special handling)
panelsWithLegend = {'Dtrans_Comparison', 'Drot_Comparison'};

% Export each panel as standalone figure
fig1Axes = {ax11, ax12, ax13, ax14; ax21, ax22, pax, ax24; ax31, ax32, ax33, ax34};
for row = 1:3
    for col = 1:4
        panelName = fig1PanelNames{row, col};
        panelAx = fig1Axes{row, col};
        
        % Create invisible figure for export
        figTemp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 650]);
        
        % Copy the axes to the new figure
        if strcmp(class(panelAx), 'matlab.graphics.axis.PolarAxes')
            % Handle polar axes - need to also add the r=frame annotation
            newAx = copyobj(panelAx, figTemp);
            newAx.Position = [0.15 0.18 0.7 0.7];
            newAx.FontSize = FMT_FONTSIZE_TICK;
            newAx.FontWeight = FMT_FONTWEIGHT;
            
            % Add "r = frame" annotation for polar plot
            annotation(figTemp, 'textbox', [0.15, 0.08, 0.2, 0.05], ...
                'String', 'r = frame', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'left');
        else
            newAx = copyobj(panelAx, figTemp);
            
            % Check if this panel needs a colorbar
            if any(strcmp(panelName, panelsWithColorbar))
                newAx.Position = [0.12 0.15 0.65 0.75];
                % Add colorbar
                cb = colorbar(newAx);
                ylabel(cb, 'Frame', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
                cb.FontSize = FMT_FONTSIZE_TICK;
                cb.FontWeight = FMT_FONTWEIGHT;
                colormap(newAx, jet);
                caxis(newAx, [minFrame maxFrame]);
            elseif any(strcmp(panelName, panelsWithLegend))
                % Panels with legends need room for the legend
                newAx.Position = [0.15 0.15 0.75 0.75];
                % Copy the legend if it exists
                origLeg = findobj(panelAx.Parent, 'Type', 'Legend');
                if ~isempty(origLeg)
                    % Recreate the legend
                    legLabels = origLeg(1).String;
                    leg = legend(newAx, legLabels, 'Location', 'best', ...
                        'FontSize', FMT_FONTSIZE_LEGEND, 'FontWeight', FMT_FONTWEIGHT);
                end
            else
                newAx.Position = [0.15 0.15 0.75 0.75];
            end
            
            % Ensure font sizes are applied
            newAx.FontSize = FMT_FONTSIZE_TICK;
            newAx.FontWeight = FMT_FONTWEIGHT;
        end
        
        % Export at high DPI
        exportFile = fullfile(fig1PanelDir, sprintf('%s_%s.png', scenarioName, panelName));
        exportgraphics(figTemp, exportFile, 'Resolution', FMT_EXPORT_DPI);
        
        close(figTemp);
    end
end
fprintf('  Figure 1 panels exported (%d panels at %d DPI)\n', 12, FMT_EXPORT_DPI);
end  % end saveFigures block for Figure 1

%% === Figure 2: Linking Cost Matrix Visualization ===
costMatLinkFile = fullfile(outputDir, 'costMatrix_linking.mat');
if exist(costMatLinkFile, 'file')
    fprintf('  Loading linking cost matrix for visualization...\n');
    costMatLinkData = load(costMatLinkFile);
    
    % Create linking cost matrix visualization figure
    fig2 = figure('Position', [100, 100, 1500, 450], 'Color', 'w', 'Name', ...
        'Cost Matrix Visualization (Frame-to-Frame Linking)');
    
    % --- Panel 1: Final cost matrix (spatial + orientation) ---
    ax1 = subplot(1, 4, 1);
    if ~isempty(costMatLinkData.costMatFinal) && nnz(costMatLinkData.costMatFinal) > 0
        [ii, jj, vv] = find(costMatLinkData.costMatFinal);
        scatter(jj, ii, 30, log10(vv), 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax1, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Detections (frame t+1)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Tracks (frame t)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatLinkData.costMatFinal, 2)+1]);
        ylim([0 size(costMatLinkData.costMatFinal, 1)+1]);
        applyStandardFormatting(gca);
        
        % Update title based on wSpatial
        if isfield(costMatLinkData, 'parameters') && costMatLinkData.parameters.wSpatial == 0
            title(sprintf('Final Cost (ORIENT-ONLY)\n%d valid links', nnz(costMatLinkData.costMatFinal)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        else
            title(sprintf('Final Cost Matrix (6D)\n%d valid links', nnz(costMatLinkData.costMatFinal)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        end
        grid on;
    else
        text(0.5, 0.5, 'No valid links', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Final Cost Matrix (6D)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 2: Spatial-only cost matrix (unweighted, for reference) ---
    ax2 = subplot(1, 4, 2);
    if isfield(costMatLinkData, 'costMatSpatial') && ~isempty(costMatLinkData.costMatSpatial) && nnz(costMatLinkData.costMatSpatial) > 0
        [ii, jj, vv] = find(costMatLinkData.costMatSpatial);
        scatter(jj, ii, 30, log10(vv), 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax2, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Detections (frame t+1)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Tracks (frame t)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatLinkData.costMatSpatial, 2)+1]);
        ylim([0 size(costMatLinkData.costMatSpatial, 1)+1]);
        applyStandardFormatting(gca);
        
        % Note if wSpatial = 0, this is just for reference
        if isfield(costMatLinkData, 'parameters') && costMatLinkData.parameters.wSpatial == 0
            title(sprintf('Spatial Cost (UNUSED, w_s=0)\n%d valid links', nnz(costMatLinkData.costMatSpatial)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        else
            title(sprintf('Spatial-Only Cost Matrix\n%d valid links', nnz(costMatLinkData.costMatSpatial)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        end
        grid on;
    else
        text(0.5, 0.5, 'No spatial cost matrix saved', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Spatial-Only Cost Matrix', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 3: Orientation cost contribution (log scale) ---
    ax3 = subplot(1, 4, 3);
    if isfield(costMatLinkData, 'costMatSpatial') && ~isempty(costMatLinkData.costMatSpatial) && ...
       ~isempty(costMatLinkData.costMatFinal) && nnz(costMatLinkData.costMatFinal) > 0
        
        % Get the weight parameters
        if isfield(costMatLinkData, 'parameters')
            wSpatialViz = costMatLinkData.parameters.wSpatial;
        else
            wSpatialViz = 1.0;
        end
        
        % Compute orientation contribution correctly by operating on full matrices
        % finalCost = wSpatial * spatialCost + orientCost
        % Therefore: orientCost = finalCost - wSpatial * spatialCost
        costFinal = full(costMatLinkData.costMatFinal);
        costSpatial = full(costMatLinkData.costMatSpatial);
        
        % Ensure same size
        [nRows, nCols] = size(costFinal);
        if ~isequal(size(costSpatial), [nRows, nCols])
            % Pad or trim to match
            costSpatialPadded = zeros(nRows, nCols);
            [rS, cS] = size(costSpatial);
            costSpatialPadded(1:min(rS,nRows), 1:min(cS,nCols)) = ...
                costSpatial(1:min(rS,nRows), 1:min(cS,nCols));
            costSpatial = costSpatialPadded;
        end
        
        % Compute orientation contribution element-wise
        orientContribMat = costFinal - wSpatialViz * costSpatial;
        
        % Find non-zero entries in the final cost matrix for plotting
        [ii, jj, vvFinal] = find(costMatLinkData.costMatFinal);
        
        % Get corresponding orientation contributions
        orientContribVals = zeros(length(ii), 1);
        for k = 1:length(ii)
            orientContribVals(k) = orientContribMat(ii(k), jj(k));
        end
        
        % Use log10 for visualization (handle zeros/negatives)
        % For negative values (orientation helping), show as negative log
        orientContribVals_log = sign(orientContribVals) .* log10(max(abs(orientContribVals), 1e-6));
        
        scatter(jj, ii, 30, orientContribVals_log, 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax3, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Orient. Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Detections (frame t+1)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Tracks (frame t)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatLinkData.costMatFinal, 2)+1]);
        ylim([0 size(costMatLinkData.costMatFinal, 1)+1]);
        applyStandardFormatting(gca);
        title(sprintf('Orientation Contribution\nmean=%.4f, max=%.4f', ...
            mean(orientContribVals), max(orientContribVals)), 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        grid on;
    else
        text(0.5, 0.5, 'Cannot compute difference', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Orientation Contribution', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 4: Cost distributions ---
    if ~isempty(costMatLinkData.costMatFinal) && nnz(costMatLinkData.costMatFinal) > 0 && ...
       isfield(costMatLinkData, 'costMatSpatial') && ~isempty(costMatLinkData.costMatSpatial)
        
        % Get the weight parameters
        if isfield(costMatLinkData, 'parameters')
            wSpatialViz = costMatLinkData.parameters.wSpatial;
        else
            wSpatialViz = 1.0;
        end
        
        % Compute orientation contribution correctly by operating on full matrices
        costFinal = full(costMatLinkData.costMatFinal);
        costSpatial = full(costMatLinkData.costMatSpatial);
        
        % Ensure same size
        [nRows, nCols] = size(costFinal);
        if ~isequal(size(costSpatial), [nRows, nCols])
            costSpatialPadded = zeros(nRows, nCols);
            [rS, cS] = size(costSpatial);
            costSpatialPadded(1:min(rS,nRows), 1:min(cS,nCols)) = ...
                costSpatial(1:min(rS,nRows), 1:min(cS,nCols));
            costSpatial = costSpatialPadded;
        end
        
        % Find non-zero entries in final cost matrix
        [ii, jj] = find(costMatLinkData.costMatFinal);
        costsFinalVec = zeros(length(ii), 1);
        costsSpatialVec = zeros(length(ii), 1);
        orientContribVec = zeros(length(ii), 1);
        for k = 1:length(ii)
            costsFinalVec(k) = costFinal(ii(k), jj(k));
            costsSpatialVec(k) = costSpatial(ii(k), jj(k));
            orientContribVec(k) = costsFinalVec(k) - wSpatialViz * costsSpatialVec(k);
        end
        
        % Top: Spatial costs (raw, unweighted)
        ax4a = subplot(2, 4, 4);
        histogram(costsSpatialVec, 30, 'FaceColor', [0.2 0.6 0.8], ...
            'FaceAlpha', 0.8, 'EdgeColor', 'w');
        xlabel('Spatial Cost', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Spatial: %.1f ± %.1f', mean(costsSpatialVec), std(costsSpatialVec)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyStandardFormatting(gca);
        grid on;
        
        % Bottom: Orientation contribution
        ax4b = subplot(2, 4, 8);
        histogram(orientContribVec, 30, 'FaceColor', [0.9 0.3 0.3], ...
            'FaceAlpha', 0.8, 'EdgeColor', 'w');
        xlabel('Orientation Cost Added', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Orient: %.4f ± %.4f', mean(orientContribVec), std(orientContribVec)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyStandardFormatting(gca);
        grid on;
        
        % Safe xlim - handle cases where orientContribVec has no positive values
        maxOrient = max(orientContribVec);
        minOrient = min(orientContribVec);
        if maxOrient > 0 && isfinite(maxOrient)
            xlim([max(0, minOrient * 1.1), maxOrient * 1.1]);
        elseif minOrient < maxOrient
            xlim([minOrient * 1.1, maxOrient * 1.1]);
        end
    end
    
    % Add parameters text (including orientKalmanLevel)
    if isfield(costMatLinkData, 'parameters')
        p = costMatLinkData.parameters;
        % Check for orientKalmanLevel or legacy useOrientKalman
        if isfield(p, 'orientKalmanLevel')
            levelStr = sprintf('Level %d', p.orientKalmanLevel);
        elseif isfield(p, 'useOrientKalman') && p.useOrientKalman
            levelStr = 'Pseudo-Kalman';
        else
            levelStr = 'Static';
        end
        paramStr = sprintf('Frame %d→%d | w_{spatial}=%.1f, w_{orient}=%.2f, w_{\\omega}=%.2f, %s', ...
            costMatLinkData.iFrame, costMatLinkData.iFrame+1, ...
            p.wSpatial, p.wOrient, p.wOmega, levelStr);
        sgtitle(sprintf('Frame-to-Frame Linking Cost Matrix | %s', paramStr), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % Save Figure 2
    if saveFigures
        saveas(fig2, fullfile(outputDir, sprintf('demo_%s_costMatrix_linking.png', scenarioName)));
        fprintf('  Figure 2 saved: demo_%s_costMatrix_linking.png\n', scenarioName);
    end
    
    %% === Export Individual Figure 2 Panels at 600 DPI ===
    if saveFigures
    fig2PanelDir = fullfile(outputDir, 'Demo_Figure_2_Panels');
    if ~exist(fig2PanelDir, 'dir')
        mkdir(fig2PanelDir);
    end
    fprintf('  Exporting Figure 2 individual panels to: %s\n', fig2PanelDir);
    
    % Panel names for Figure 2
    fig2PanelNames = {'Final_Cost_Matrix', 'Spatial_Cost_Matrix', 'Orientation_Contribution', ...
                      'Histogram_Spatial', 'Histogram_Orientation'};
    fig2Axes = {ax1, ax2, ax3};
    if exist('ax4a', 'var'), fig2Axes{4} = ax4a; end
    if exist('ax4b', 'var'), fig2Axes{5} = ax4b; end
    
    % Panels that need colorbars (the scatter plot cost matrices)
    fig2PanelsWithColorbar = [1, 2, 3];  % First three panels have colorbars
    fig2ColorbarLabels = {'log_{10}(Cost)', 'log_{10}(Cost)', 'log_{10}(Orient. Cost)'};
    
    for i = 1:length(fig2Axes)
        if i <= length(fig2PanelNames) && ~isempty(fig2Axes{i}) && isvalid(fig2Axes{i})
            panelName = fig2PanelNames{i};
            panelAx = fig2Axes{i};
            
            % Create invisible figure for export
            figTemp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 650]);
            
            % Copy the axes to the new figure
            newAx = copyobj(panelAx, figTemp);
            
            % Check if this panel needs a colorbar
            if ismember(i, fig2PanelsWithColorbar)
                newAx.Position = [0.12 0.15 0.62 0.75];
                
                % Set colormap and preserve color limits from original
                colormap(newAx, 'turbo');
                if isprop(panelAx, 'CLim')
                    newAx.CLim = panelAx.CLim;
                end
                
                % Create colorbar and set its properties
                cb = colorbar(newAx);
                cb.Label.String = fig2ColorbarLabels{i};
                cb.Label.FontSize = FMT_FONTSIZE_TICK;
                cb.Label.FontWeight = FMT_FONTWEIGHT;
                cb.FontSize = FMT_FONTSIZE_TICK;
                cb.FontWeight = FMT_FONTWEIGHT;
            else
                newAx.Position = [0.15 0.15 0.75 0.75];
            end
            
            % Ensure font sizes are applied
            newAx.FontSize = FMT_FONTSIZE_TICK;
            newAx.FontWeight = FMT_FONTWEIGHT;
            if ~isempty(newAx.Title)
                newAx.Title.FontSize = FMT_FONTSIZE_TITLE;
            end
            if ~isempty(newAx.XLabel)
                newAx.XLabel.FontSize = FMT_FONTSIZE_LABEL;
            end
            if ~isempty(newAx.YLabel)
                newAx.YLabel.FontSize = FMT_FONTSIZE_LABEL;
            end
            
            % Export at high DPI
            exportFile = fullfile(fig2PanelDir, sprintf('%s_%s.png', scenarioName, panelName));
            exportgraphics(figTemp, exportFile, 'Resolution', FMT_EXPORT_DPI);
            
            close(figTemp);
        end
    end
    fprintf('  Figure 2 panels exported at %d DPI\n', FMT_EXPORT_DPI);
    end  % end saveFigures block for Figure 2
    
else
    fprintf('  Linking cost matrix file not found (tracking may have been skipped).\n');
end

%% === Figure 3: Gap Closing Cost Matrix Visualization ===
costMatFile = fullfile(outputDir, 'costMatrix_gapClosing.mat');
if exist(costMatFile, 'file')
    fprintf('  Loading gap closing cost matrix for visualization...\n');
    costMatData = load(costMatFile);
    
    % Create gap closing cost matrix visualization figure
    fig3 = figure('Position', [100, 100, 1500, 450], 'Color', 'w', 'Name', ...
        'Cost Matrix Visualization (Gap Closing)');
    
    % --- Panel 1: Final cost matrix (spatial + orientation) ---
    ax1 = subplot(1, 4, 1);
    if ~isempty(costMatData.costMatFinal) && nnz(costMatData.costMatFinal) > 0
        [ii, jj, vv] = find(costMatData.costMatFinal);
        scatter(jj, ii, 20, log10(vv), 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax1, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Track Starts (j)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Track Ends (i)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatData.costMatFinal, 2)+1]);
        ylim([0 size(costMatData.costMatFinal, 1)+1]);
        applyStandardFormatting(gca);
        
        % Update title based on wSpatial
        if isfield(costMatData, 'parameters') && costMatData.parameters.wSpatial == 0
            title(sprintf('Final Cost (ORIENT-ONLY)\n%d valid links', nnz(costMatData.costMatFinal)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        else
            title(sprintf('Final Cost Matrix (6D)\n%d valid links', nnz(costMatData.costMatFinal)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        end
        grid on;
    else
        text(0.5, 0.5, 'No valid links', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Final Cost Matrix (6D)', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 2: Spatial-only cost matrix (unweighted, for reference) ---
    ax2 = subplot(1, 4, 2);
    if isfield(costMatData, 'costMatSpatial') && ~isempty(costMatData.costMatSpatial) && nnz(costMatData.costMatSpatial) > 0
        [ii, jj, vv] = find(costMatData.costMatSpatial);
        scatter(jj, ii, 20, log10(vv), 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax2, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Track Starts (j)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Track Ends (i)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatData.costMatSpatial, 2)+1]);
        ylim([0 size(costMatData.costMatSpatial, 1)+1]);
        applyStandardFormatting(gca);
        
        % Note if wSpatial = 0, this is just for reference
        if isfield(costMatData, 'parameters') && costMatData.parameters.wSpatial == 0
            title(sprintf('Spatial Cost (UNUSED, w_s=0)\n%d valid links', nnz(costMatData.costMatSpatial)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        else
            title(sprintf('Spatial-Only Cost Matrix\n%d valid links', nnz(costMatData.costMatSpatial)), ...
                'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        end
        grid on;
    else
        text(0.5, 0.5, 'No spatial cost matrix saved', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Spatial-Only Cost Matrix', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 3: Orientation cost contribution (log scale) ---
    ax3 = subplot(1, 4, 3);
    if isfield(costMatData, 'costMatSpatial') && ~isempty(costMatData.costMatSpatial) && ...
       ~isempty(costMatData.costMatFinal) && nnz(costMatData.costMatFinal) > 0
        
        % Get the weight parameters
        if isfield(costMatData, 'parameters')
            wSpatialViz = costMatData.parameters.wSpatial;
        else
            wSpatialViz = 1.0;
        end
        
        % Compute orientation contribution correctly by operating on full matrices
        % finalCost = wSpatial * spatialCost + orientCost
        % Therefore: orientCost = finalCost - wSpatial * spatialCost
        costFinal = full(costMatData.costMatFinal);
        costSpatial = full(costMatData.costMatSpatial);
        
        % Ensure same size
        [nRows, nCols] = size(costFinal);
        if ~isequal(size(costSpatial), [nRows, nCols])
            % Pad or trim to match
            costSpatialPadded = zeros(nRows, nCols);
            [rS, cS] = size(costSpatial);
            costSpatialPadded(1:min(rS,nRows), 1:min(cS,nCols)) = ...
                costSpatial(1:min(rS,nRows), 1:min(cS,nCols));
            costSpatial = costSpatialPadded;
        end
        
        % Compute orientation contribution element-wise
        orientContribMat = costFinal - wSpatialViz * costSpatial;
        
        % Find non-zero entries in the final cost matrix for plotting
        [ii, jj, vvFinal] = find(costMatData.costMatFinal);
        
        % Get corresponding orientation contributions
        orientContribVals = zeros(length(ii), 1);
        for k = 1:length(ii)
            orientContribVals(k) = orientContribMat(ii(k), jj(k));
        end
        
        % Use log10 for visualization (handle zeros/negatives)
        % For negative values (orientation helping), show as negative log
        orientContribVals_log = sign(orientContribVals) .* log10(max(abs(orientContribVals), 1e-6));
        
        scatter(jj, ii, 20, orientContribVals_log, 'filled', 'MarkerFaceAlpha', 0.8);
        set(gca, 'YDir', 'normal');
        colormap(ax3, 'turbo');
        cb = colorbar; ylabel(cb, 'log_{10}(Orient. Cost)', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        xlabel('Track Starts (j)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Track Ends (i)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        xlim([0 size(costMatData.costMatFinal, 2)+1]);
        ylim([0 size(costMatData.costMatFinal, 1)+1]);
        applyStandardFormatting(gca);
        title(sprintf('Orientation Contribution\nmean=%.4f, max=%.4f', ...
            mean(orientContribVals), max(orientContribVals)), 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        grid on;
    else
        text(0.5, 0.5, 'Cannot compute difference', 'HorizontalAlignment', 'center', 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Orientation Contribution', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % --- Panel 4: Cost distributions (separate subplots) ---
    if ~isempty(costMatData.costMatFinal) && nnz(costMatData.costMatFinal) > 0 && ...
       isfield(costMatData, 'costMatSpatial') && ~isempty(costMatData.costMatSpatial)
        
        % Get the weight parameters
        if isfield(costMatData, 'parameters')
            wSpatialViz = costMatData.parameters.wSpatial;
        else
            wSpatialViz = 1.0;
        end
        
        % Compute orientation contribution correctly by operating on full matrices
        costFinal = full(costMatData.costMatFinal);
        costSpatial = full(costMatData.costMatSpatial);
        
        % Ensure same size
        [nRows, nCols] = size(costFinal);
        if ~isequal(size(costSpatial), [nRows, nCols])
            costSpatialPadded = zeros(nRows, nCols);
            [rS, cS] = size(costSpatial);
            costSpatialPadded(1:min(rS,nRows), 1:min(cS,nCols)) = ...
                costSpatial(1:min(rS,nRows), 1:min(cS,nCols));
            costSpatial = costSpatialPadded;
        end
        
        % Find non-zero entries in final cost matrix
        [ii, jj] = find(costMatData.costMatFinal);
        costsFinalVec = zeros(length(ii), 1);
        costsSpatialVec = zeros(length(ii), 1);
        orientContribVec = zeros(length(ii), 1);
        for k = 1:length(ii)
            costsFinalVec(k) = costFinal(ii(k), jj(k));
            costsSpatialVec(k) = costSpatial(ii(k), jj(k));
            orientContribVec(k) = costsFinalVec(k) - wSpatialViz * costsSpatialVec(k);
        end
        
        % Top: Spatial costs (raw, unweighted)
        ax4a = subplot(2, 4, 4);
        histogram(costsSpatialVec, 30, 'FaceColor', [0.2 0.6 0.8], ...
            'FaceAlpha', 0.8, 'EdgeColor', 'w');
        xlabel('Spatial Cost', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Spatial: %.1f ± %.1f', mean(costsSpatialVec), std(costsSpatialVec)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyStandardFormatting(gca);
        grid on;
        
        % Bottom: Orientation contribution
        ax4b = subplot(2, 4, 8);
        histogram(orientContribVec, 30, 'FaceColor', [0.9 0.3 0.3], ...
            'FaceAlpha', 0.8, 'EdgeColor', 'w');
        xlabel('Orientation Cost Added', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Orient: %.4f ± %.4f', mean(orientContribVec), std(orientContribVec)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyStandardFormatting(gca);
        grid on;
        
        % Safe xlim - handle cases where orientContribVec has no positive values
        maxOrient = max(orientContribVec);
        minOrient = min(orientContribVec);
        if maxOrient > 0 && isfinite(maxOrient)
            xlim([max(0, minOrient * 1.1), maxOrient * 1.1]);
        elseif minOrient < maxOrient
            xlim([minOrient * 1.1, maxOrient * 1.1]);
        end
    end
    
    % Add parameters text as sgtitle (including orientKalmanLevel info)
    if isfield(costMatData, 'parameters')
        p = costMatData.parameters;
        % Check for orientKalmanLevel or legacy useOrientKalman
        if isfield(p, 'orientKalmanLevel')
            levelStr = sprintf('Level %d', p.orientKalmanLevel);
        elseif isfield(p, 'useOrientKalman') && p.useOrientKalman
            levelStr = 'Pseudo-Kalman';
        else
            levelStr = 'Static';
        end
        paramStr = sprintf('w_{spatial}=%.1f, w_{orient}=%.2f, w_{\\omega}=%.2f, %s', ...
            p.wSpatial, p.wOrient, p.wOmega, levelStr);
        sgtitle(sprintf('Gap Closing Cost Matrix | %s', paramStr), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    
    % Save Figure
    if saveFigures
        saveas(fig3, fullfile(outputDir, sprintf('demo_%s_costMatrix.png', scenarioName)));
        fprintf('  Figure 3 saved: demo_%s_costMatrix.png\n', scenarioName);
    end
    
    %% === Export Individual Figure 3 Panels at 600 DPI ===
    if saveFigures
    fig3PanelDir = fullfile(outputDir, 'Demo_Figure_3_Panels');
    if ~exist(fig3PanelDir, 'dir')
        mkdir(fig3PanelDir);
    end
    fprintf('  Exporting Figure 3 individual panels to: %s\n', fig3PanelDir);
    
    % Panel names for Figure 3
    fig3PanelNames = {'Final_Cost_Matrix', 'Spatial_Cost_Matrix', 'Orientation_Contribution', ...
                      'Histogram_Spatial', 'Histogram_Orientation'};
    fig3Axes = {ax1, ax2, ax3};
    if exist('ax4a', 'var'), fig3Axes{4} = ax4a; end
    if exist('ax4b', 'var'), fig3Axes{5} = ax4b; end
    
    % Panels that need colorbars (the scatter plot cost matrices)
    fig3PanelsWithColorbar = [1, 2, 3];  % First three panels have colorbars
    fig3ColorbarLabels = {'log_{10}(Cost)', 'log_{10}(Cost)', 'log_{10}(Orient. Cost)'};
    
    for i = 1:length(fig3Axes)
        if i <= length(fig3PanelNames) && ~isempty(fig3Axes{i}) && isvalid(fig3Axes{i})
            panelName = fig3PanelNames{i};
            panelAx = fig3Axes{i};
            
            % Create invisible figure for export
            figTemp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 650]);
            
            % Copy the axes to the new figure
            newAx = copyobj(panelAx, figTemp);
            
            % Check if this panel needs a colorbar
            if ismember(i, fig3PanelsWithColorbar)
                newAx.Position = [0.12 0.15 0.62 0.75];
                
                % Set colormap and preserve color limits from original
                colormap(newAx, 'turbo');
                if isprop(panelAx, 'CLim')
                    newAx.CLim = panelAx.CLim;
                end
                
                % Create colorbar and set its properties
                cb = colorbar(newAx);
                cb.Label.String = fig3ColorbarLabels{i};
                cb.Label.FontSize = FMT_FONTSIZE_TICK;
                cb.Label.FontWeight = FMT_FONTWEIGHT;
                cb.FontSize = FMT_FONTSIZE_TICK;
                cb.FontWeight = FMT_FONTWEIGHT;
            else
                newAx.Position = [0.15 0.15 0.75 0.75];
            end
            
            % Ensure font sizes are applied
            newAx.FontSize = FMT_FONTSIZE_TICK;
            newAx.FontWeight = FMT_FONTWEIGHT;
            if ~isempty(newAx.Title)
                newAx.Title.FontSize = FMT_FONTSIZE_TITLE;
            end
            if ~isempty(newAx.XLabel)
                newAx.XLabel.FontSize = FMT_FONTSIZE_LABEL;
            end
            if ~isempty(newAx.YLabel)
                newAx.YLabel.FontSize = FMT_FONTSIZE_LABEL;
            end
            
            % Export at high DPI
            exportFile = fullfile(fig3PanelDir, sprintf('%s_%s.png', scenarioName, panelName));
            exportgraphics(figTemp, exportFile, 'Resolution', FMT_EXPORT_DPI);
            
            close(figTemp);
        end
    end
    fprintf('  Figure 3 panels exported at %d DPI\n', FMT_EXPORT_DPI);
    end  % end saveFigures block for Figure 3
    
else
    fprintf('  Cost matrix file not found (tracking may have been skipped).\n');
end

%% ========================================================================
%  DONE
%  ========================================================================
fprintf('\n============================================\n');
fprintf('  Demo complete!\n');
fprintf('  Output: %s\n', outputDir);
fprintf('  Files:\n');
d = dir(fullfile(outputDir, '*.*'));
for i = 1:length(d)
    if ~d(i).isdir
        fprintf('    %s\n', d(i).name);
    end
end
fprintf('============================================\n');

%% ========================================================================
%  LOCAL HELPER FUNCTION
%  ========================================================================
function plotTimeColoredTraj(ax, xData, yData, frameData, minFrame, maxFrame)
%PLOTTIMECOLOREDTRAJ Plot trajectory segments with jet colormap encoding time.
%   Each segment between consecutive points is colored by its frame number.

    if isempty(xData) || length(xData) < 2
        return;
    end
    
    % Sort by frame
    [frameData, sortIdx] = sort(frameData);
    xData = xData(sortIdx);
    yData = yData(sortIdx);
    
    % Normalize frames to [0, 1] for colormap
    normFrames = (frameData - minFrame) / max(maxFrame - minFrame, 1);
    
    % Get jet colormap
    cmap = jet(256);
    
    % Plot each segment with its own color
    for j = 1:(length(xData) - 1)
        % Color index based on average frame of segment
        avgNorm = (normFrames(j) + normFrames(j+1)) / 2;
        colorIdx = max(1, min(256, round(avgNorm * 255) + 1));
        segColor = cmap(colorIdx, :);
        
        plot(ax, xData(j:j+1), yData(j:j+1), '-', 'Color', segColor, 'LineWidth', 1.5);
    end
end
