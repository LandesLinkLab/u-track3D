function [results] = estimateDiffusion6D_MLE(trajData, varargin)
%ESTIMATEDIFFUSION6D_MLE Estimate diffusion coefficients using MLE for 6D trajectories.
%
%   Extends Bo Shuang's MLE method (Langmuir 2013) to estimate both
%   translational and rotational diffusion coefficients from 6D-SMOLM
%   trajectories. Properly handles variable time lags caused by
%   photoblinking.
%
%   Standard MSD analysis assumes uniform time steps and gives biased
%   estimates when gaps exist. This MLE approach weights each displacement
%   by its actual time lag, giving unbiased D estimates even with heavy
%   blinking.
%
%   SYNOPSIS:
%       results = estimateDiffusion6D_MLE(trajData)
%       results = estimateDiffusion6D_MLE(trajData, 'Name', Value, ...)
%
%   INPUT:
%       trajData : [N x 9] trajectory matrix with columns:
%                  [frame, trackID, x_nm, y_nm, z_nm, theta, phi, omega, intensity]
%                  OR a cell array where each cell is one trajectory.
%
%   NAME-VALUE PARAMETERS:
%       'FrameTime'      : Time per frame in seconds. Default: 0.05 (50 ms).
%       'LocPrecisionXY' : Localization precision in nm (x,y). Default: 20.
%       'LocPrecisionZ'  : Localization precision in nm (z). Default: 50.
%       'LocPrecisionTheta': Angular precision in rad. Default: 0.1.
%       'LocPrecisionPhi'  : Angular precision in rad. Default: 0.2.
%       'MinDisplacements' : Minimum number of displacements per trajectory. Default: 3.
%       'Verbose'        : Print results. Default: true.
%
%   OUTPUT:
%       results : Struct array (one per trajectory) with fields:
%           .trackID        : Track ID
%           .nDisplacements : Number of displacement steps
%           .D_trans        : Translational diffusion coefficient (nm^2/s)
%           .D_trans_err    : Estimated uncertainty in D_trans
%           .D_rot          : Rotational diffusion coefficient (rad^2/s)
%           .D_rot_err      : Estimated uncertainty in D_rot
%           .D_xy           : 2D lateral diffusion coefficient (nm^2/s)
%           .meanTimeLag    : Mean time lag between steps (s)
%           .fracGaps       : Fraction of steps with time lag > 1 frame
%
%   ALGORITHM:
%       This implements MLE(1) from Bo Shuang's Langmuir 2013 paper (Eq. 3-4).
%
%       For TRANSLATIONAL diffusion, the log-likelihood is (Eq. 3):
%
%           L(Δs) = Σ ln{ 1/√(4πD_s·n_i·Δt) · exp(-Δs_i²/(4D_s·n_i·Δt)) }
%
%       Maximizing gives the closed-form MLE(1) estimator (Eq. 4):
%
%           D_s^(1) = (1/(2N·Δt)) · Σ (Δs_i²/n_i)
%
%       where:
%           Δs_i = spatial displacement for step i (nm)
%           n_i  = number of frames for step i (1 normally, k+1 if blinked k frames)
%           N    = total number of measured displacements
%           Δt   = time per frame (seconds)
%
%       For ROTATIONAL diffusion, the same form applies:
%
%           L(Δα) = Σ ln{ 1/√(4πD_α·n_i·Δt) · exp(-Δα_i²/(4D_α·n_i·Δt)) }
%
%           D_α^(1) = (1/(2N·Δt)) · Σ (Δα_i²/n_i)
%
%       where:
%           Δα_i = angular displacement for step i (rad)
%                  computed as α = arccos(|μ̂₁·μ̂₂|) ∈ [0, π/2]
%
%       NOTE: MLE(1) is BIASED when localization error is significant.
%       For noisy data, use estimateDiffusion6D_DualMLE.m which implements
%       MLE(2) with covariance matrix to jointly estimate D and σ².
%
%   REFERENCES:
%       Shuang, B. et al. "Improved Analysis for Determining Diffusion
%       Coefficients from Short Single-Molecule Trajectories with
%       Photoblinking." Langmuir 29, 228-234 (2013).
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'trajData');
addParameter(p, 'FrameTime', 0.05, @isnumeric);        % 50 ms
addParameter(p, 'LocPrecisionXY', 20, @isnumeric);      % nm
addParameter(p, 'LocPrecisionZ', 50, @isnumeric);       % nm
addParameter(p, 'LocPrecisionTheta', 0.1, @isnumeric);  % rad
addParameter(p, 'LocPrecisionPhi', 0.2, @isnumeric);    % rad
addParameter(p, 'MinDisplacements', 3, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
parse(p, trajData, varargin{:});

dt       = p.Results.FrameTime;
sigmaXY  = p.Results.LocPrecisionXY;
sigmaZ   = p.Results.LocPrecisionZ;
sigmaT   = p.Results.LocPrecisionTheta;
sigmaP   = p.Results.LocPrecisionPhi;
minDisp  = p.Results.MinDisplacements;
verbose  = p.Results.Verbose;

%% --- Convert to cell array of individual trajectories ---
if istable(trajData)
    trajData = table2array(trajData);
end

if iscell(trajData)
    trajs = trajData;
else
    uniqueIDs = unique(trajData(:, 2));
    trajs = cell(length(uniqueIDs), 1);
    for i = 1:length(uniqueIDs)
        trajs{i} = trajData(trajData(:, 2) == uniqueIDs(i), :);
    end
end

nTrajs = length(trajs);

%% --- Process each trajectory ---
results = struct('trackID', {}, 'nDisplacements', {}, ...
    'D_trans', {}, 'D_trans_err', {}, ...
    'D_rot', {}, 'D_rot_err', {}, ...
    'D_xy', {}, 'meanTimeLag', {}, 'fracGaps', {});

nValid = 0;

for iTraj = 1:nTrajs
    traj = trajs{iTraj};
    
    if isempty(traj) || size(traj, 1) < minDisp + 1
        continue;
    end
    
    % Sort by frame
    traj = sortrows(traj, 1);
    
    frames = traj(:, 1);
    x      = traj(:, 3);
    y      = traj(:, 4);
    z      = traj(:, 5);
    theta  = traj(:, 6);
    phi    = traj(:, 7);
    % omega column 8 is not used for rotational D estimation
    
    nSteps = size(traj, 1) - 1;
    
    if nSteps < minDisp
        continue;
    end
    
    % --- Compute displacements and time lags ---
    dFrames = diff(frames);           % Time lag in frames
    dtVec   = dFrames * dt;           % Time lag in seconds
    
    % Skip zero time lags (shouldn't happen, but safety check)
    validSteps = dFrames > 0;
    if sum(validSteps) < minDisp
        continue;
    end
    
    % Spatial displacements
    dx = diff(x);
    dy = diff(y);
    dz = diff(z);
    dr2_3d = dx.^2 + dy.^2 + dz.^2;  % 3D squared displacement
    dr2_2d = dx.^2 + dy.^2;           % 2D squared displacement
    
    % Angular displacements (using angularDistance for proper symmetry)
    dAlpha = zeros(nSteps, 1);
    for iStep = 1:nSteps
        [dAlpha(iStep), ~, ~] = angularDistance(...
            theta(iStep), phi(iStep), [], ...
            theta(iStep+1), phi(iStep+1), []);
    end
    dAlpha2 = dAlpha.^2;
    
    % Apply validity mask
    dtVec   = dtVec(validSteps);
    dr2_3d  = dr2_3d(validSteps);
    dr2_2d  = dr2_2d(validSteps);
    dAlpha2 = dAlpha2(validSteps);
    dFramesValid = dFrames(validSteps);
    N = sum(validSteps);
    
    % --- MLE(1) for translational diffusion (3D) ---
    % Bo Shuang Eq. 4: D = (1/(2N·dt)) · Σ(Δ_i²/n_i)
    %
    % For 3D motion, we have d=3 dimensions, so:
    %   D_s^(1) = (1/(2·3·N·dt)) · Σ(Δs_i²/n_i)
    %
    % where n_i = dFrames (number of frames for step i, handles blinking)
    %
    % NOTE: This is the UNCORRECTED MLE(1). For data with significant
    % localization error, the estimate will be biased high by ~σ²/dt.
    % Use MLE(2) or subtract correction term if σ is known.
    
    meanDt = mean(dtVec);  % Mean time lag for noise correction
    
    % MLE(1) estimate - matches Bo's Eq. 4
    sum_r2_over_n = sum(dr2_3d ./ dFramesValid);
    D_trans_raw = sum_r2_over_n / (2 * 3 * N * dt);
    
    % Optional correction for known localization precision
    % Localization variance per displacement (2 endpoints, each with noise)
    locVar3D = 2 * (2*sigmaXY^2 + sigmaZ^2);  % nm^2
    D_trans = D_trans_raw - locVar3D / (2 * 3 * meanDt);
    D_trans = max(D_trans, 0);  % Floor at 0
    
    % Uncertainty estimate (simplified - RSD ~ 1/√N for MLE)
    D_trans_err = abs(D_trans) / sqrt(N);
    
    % 2D lateral diffusion (same MLE(1) approach, d=2)
    sum_r2_over_n_2d = sum(dr2_2d ./ dFramesValid);
    D_xy_raw = sum_r2_over_n_2d / (2 * 2 * N * dt);
    locVar2D = 2 * (2*sigmaXY^2);  % nm^2
    D_xy = D_xy_raw - locVar2D / (2 * 2 * meanDt);
    D_xy = max(D_xy, 0);
    
    % --- MLE(1) for rotational diffusion ---
    % Same form as translational but on sphere (effectively d=1 for angle):
    %   D_α^(1) = (1/(2N·dt)) · Σ(Δα_i²/n_i)
    %
    % Angular localization variance (2 endpoints)
    locVarRot = 2 * (sigmaT^2 + sigmaP^2);  % rad^2
    
    sum_alpha2_over_n = sum(dAlpha2 ./ dFramesValid);
    D_rot_raw = sum_alpha2_over_n / (2 * N * dt);
    D_rot = D_rot_raw - locVarRot / (2 * meanDt);
    D_rot = max(D_rot, 0);
    D_rot_err = abs(D_rot) / sqrt(N);
    
    % --- Fraction of gaps ---
    fracGaps = sum(dFramesValid > 1) / N;
    
    % --- Store results ---
    nValid = nValid + 1;
    results(nValid).trackID        = traj(1, 2);
    results(nValid).nDisplacements = N;
    results(nValid).D_trans        = D_trans;
    results(nValid).D_trans_err    = D_trans_err;
    results(nValid).D_rot          = D_rot;
    results(nValid).D_rot_err      = D_rot_err;
    results(nValid).D_xy           = D_xy;
    results(nValid).meanTimeLag    = meanDt;
    results(nValid).fracGaps       = fracGaps;
end

%% --- Summary ---
if verbose && nValid > 0
    allDtrans = [results.D_trans];
    allDrot   = [results.D_rot];
    allFracGap = [results.fracGaps];
    
    fprintf('--- estimateDiffusion6D_MLE summary ---\n');
    fprintf('  Trajectories analyzed: %d / %d\n', nValid, nTrajs);
    fprintf('  D_trans (nm^2/s): median = %.1f, mean = %.1f, range = [%.1f, %.1f]\n', ...
        median(allDtrans), mean(allDtrans), min(allDtrans), max(allDtrans));
    fprintf('  D_rot (rad^2/s):  median = %.4f, mean = %.4f, range = [%.4f, %.4f]\n', ...
        median(allDrot), mean(allDrot), min(allDrot), max(allDrot));
    fprintf('  Mean fraction gaps: %.2f\n', mean(allFracGap));
    fprintf('---------------------------------------\n');
elseif verbose
    fprintf('  No valid trajectories for diffusion analysis.\n');
end

end
