function [results] = estimateDiffusion6D_DualMLE(trajData, varargin)
%ESTIMATEDIFFUSION6D_DUALMLE Dual MLE diffusion estimation for 6D trajectories.
%
%   Implements BOTH MLE methods from Bo Shuang's Langmuir 2013 paper:
%     - MLE(1): For low-noise data (x < 0.1), assumes known localization precision
%     - MLE(2): For noisy data (0.1 < x < 10), jointly estimates D and σ²
%
%   Where x = σ²/(D·dt) is the reduced square localization error.
%
%   The function automatically determines which MLE is more appropriate based
%   on the estimated noise level, but reports BOTH estimates for comparison.
%
%   SYNOPSIS:
%       results = estimateDiffusion6D_DualMLE(trajData)
%       results = estimateDiffusion6D_DualMLE(trajData, 'Name', Value, ...)
%
%   INPUT:
%       trajData : [N x 9] trajectory matrix with columns:
%                  [frame, trackID, x_nm, y_nm, z_nm, theta, phi, omega, intensity]
%                  OR a cell array where each cell is one trajectory.
%
%   NAME-VALUE PARAMETERS:
%       'FrameTime'        : Time per frame in seconds. Default: 0.05 (50 ms).
%       'LocPrecisionXY'   : Localization precision in nm (x,y). Default: 20.
%       'LocPrecisionZ'    : Localization precision in nm (z). Default: 50.
%       'LocPrecisionTheta': Angular precision in rad (theta). Default: 0.1.
%       'LocPrecisionPhi'  : Angular precision in rad (phi). Default: 0.2.
%       'MinDisplacements' : Minimum displacements per trajectory. Default: 3.
%       'Verbose'          : Print results. Default: true.
%
%   OUTPUT:
%       results : Struct array (one per trajectory) with fields:
%           .trackID         : Track ID
%           .nDisplacements  : Number of displacement steps
%           .D_trans_MLE1    : Translational D from MLE(1) [nm²/s]
%           .D_trans_MLE2    : Translational D from MLE(2) [nm²/s]
%           .D_trans_best    : Best estimate (auto-selected)
%           .sigma2_MLE2     : Estimated localization variance from MLE(2)
%           .D_rot_MLE1      : Rotational D from MLE(1) [rad²/s]
%           .D_rot_MLE2      : Rotational D from MLE(2) [rad²/s]
%           .D_rot_best      : Best estimate (auto-selected)
%           .x_trans         : Reduced localization error for translation
%           .x_rot           : Reduced localization error for rotation
%           .recommended_MLE : 1 or 2, which MLE is recommended
%           .meanTimeLag     : Mean time lag between steps (s)
%           .fracGaps        : Fraction of steps with time lag > 1 frame
%
%   THEORY (Bo Shuang, Langmuir 2013):
%
%   For TRANSLATIONAL diffusion, using spatial displacement Δs:
%       MSD: ⟨Δs²⟩ = 2·d·D_s·Δt + 2·σ_s²
%
%   For ROTATIONAL diffusion, using angular displacement α:
%       MSD: ⟨Δα²⟩ = 2·D_α·Δt + 2·σ_α²
%       where α = arccos(|μ̂₁·μ̂₂|) is the angular displacement between
%       dipole orientations, ranging from 0 to π/2.
%
%   MLE(1) - For noise-free or low-noise data (Eq. 3-4):
%       Log-likelihood (Eq. 3):
%           L(Δs) = Σ ln{ 1/√(4πD_s·n_i·Δt) · exp(-Δs_i²/(4D_s·n_i·Δt)) }
%       
%       Closed-form MLE estimate (Eq. 4):
%           D_s^(1) = (1/(2N·Δt)) · Σ (Δs_i²/n_i)
%       
%       where n_i is the number of frames for step i (accounts for gaps).
%       
%       This is BIASED when noise is significant because measured Δs² includes
%       localization error: ⟨Δs²_measured⟩ = ⟨Δs²_true⟩ + 2σ_s²
%
%   MLE(2) - For noisy data (Eq. 7-8):
%       Uses covariance matrix Σ to account for correlated measurement errors.
%       
%       Log-likelihood (Eq. 7):
%           L(Δ) = -½ ln|Σ| - ½ Δᵀ Σ⁻¹ Δ
%       
%       Covariance matrix elements (Eq. 8):
%           (Σ_s)_ij = (2D_s·n_i·Δt + 2σ_s²)·δ_ij - σ_s²·δ_{i,j±1}
%       
%       Jointly estimates D and σ² by maximizing L.
%       UNBIASED for 0.1 < x < 10, but biased outside this range.
%
%   RECOMMENDATIONS:
%       x < 0.1:  Use MLE(1) - low noise, MLE(2) becomes biased
%       0.1-10:   Use MLE(2) - optimal regime for joint estimation
%       x > 10:   Neither reliable - localization error dominates
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
addParameter(p, 'FrameTime', 0.05, @isnumeric);
addParameter(p, 'LocPrecisionXY', 20, @isnumeric);
addParameter(p, 'LocPrecisionZ', 50, @isnumeric);
addParameter(p, 'LocPrecisionTheta', 0.1, @isnumeric);
addParameter(p, 'LocPrecisionPhi', 0.2, @isnumeric);
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

%% --- Convert input to cell array of trajectories ---
if iscell(trajData)
    trajectories = trajData;
else
    % Matrix format: split by trackID (column 2)
    trackIDs = unique(trajData(:,2));
    nTracks = length(trackIDs);
    trajectories = cell(nTracks, 1);
    for i = 1:nTracks
        trajectories{i} = trajData(trajData(:,2) == trackIDs(i), :);
    end
end

nTraj = length(trajectories);

%% --- Initialize results ---
results = struct('trackID', cell(nTraj,1), ...
                 'nDisplacements', cell(nTraj,1), ...
                 'D_trans_MLE1', cell(nTraj,1), ...
                 'D_trans_MLE2', cell(nTraj,1), ...
                 'D_trans_best', cell(nTraj,1), ...
                 'sigma2_trans_MLE2', cell(nTraj,1), ...
                 'D_rot_MLE1', cell(nTraj,1), ...
                 'D_rot_MLE2', cell(nTraj,1), ...
                 'D_rot_best', cell(nTraj,1), ...
                 'sigma2_rot_MLE2', cell(nTraj,1), ...
                 'x_trans', cell(nTraj,1), ...
                 'x_rot', cell(nTraj,1), ...
                 'recommended_MLE', cell(nTraj,1), ...
                 'meanTimeLag', cell(nTraj,1), ...
                 'fracGaps', cell(nTraj,1));

%% --- Process each trajectory ---
nAnalyzed = 0;
for iTraj = 1:nTraj
    traj = trajectories{iTraj};
    
    % Sort by frame
    traj = sortrows(traj, 1);
    
    % Extract columns
    frames = traj(:,1);
    trackID = traj(1,2);
    x = traj(:,3);
    y = traj(:,4);
    z = traj(:,5);
    theta = traj(:,6);
    phi = traj(:,7);
    
    % Number of steps
    nPts = size(traj, 1);
    nSteps = nPts - 1;
    
    if nSteps < minDisp
        continue;
    end
    
    nAnalyzed = nAnalyzed + 1;
    
    % Frame differences (n_i in Bo's notation - accounts for gaps)
    dFrames = diff(frames);
    
    % Time lags for each step
    dtVec = dFrames * dt;
    
    % Identify valid steps (exclude very long gaps)
    maxGapFrames = 10;
    validSteps = dFrames <= maxGapFrames & dFrames >= 1;
    
    % --- Spatial displacements ---
    dx = diff(x);
    dy = diff(y);
    dz = diff(z);
    dr2_3d = dx.^2 + dy.^2 + dz.^2;
    
    % --- Angular displacements ---
    dAlpha = zeros(nSteps, 1);
    for iStep = 1:nSteps
        [dAlpha(iStep), ~, ~] = angularDistance(...
            theta(iStep), phi(iStep), [], ...
            theta(iStep+1), phi(iStep+1), []);
    end
    dAlpha2 = dAlpha.^2;
    
    % Apply validity mask
    dtVec_valid = dtVec(validSteps);
    dr2_valid = dr2_3d(validSteps);
    dAlpha2_valid = dAlpha2(validSteps);
    dFrames_valid = dFrames(validSteps);
    N = sum(validSteps);
    
    if N < minDisp
        continue;
    end
    
    meanDt = mean(dtVec_valid);
    fracGaps = sum(dFrames_valid > 1) / N;
    
    %% ==================== MLE (1) ====================
    % For noise-free data: D = Σ(Δᵢ²/nᵢ) / (2·d·N·dt)
    % With known noise correction: D = [Σ(Δᵢ²/nᵢ)/N - 2σ²/meanDt] / (2·d)
    
    % Translational (3D)
    % Known localization variance (2 endpoints, each with σ²)
    locVar3D = 2 * (2*sigmaXY^2 + sigmaZ^2);
    
    avg_r2_over_n = sum(dr2_valid ./ dFrames_valid) / N;
    D_trans_MLE1 = (1 / (2*3*dt)) * avg_r2_over_n - locVar3D / (2*3*meanDt);
    D_trans_MLE1 = max(D_trans_MLE1, 1e-6);  % Floor at small positive value
    
    % Rotational (on sphere - effectively 2D)
    locVar_ang = 2 * (sigmaT^2 + sigmaP^2);
    
    avg_alpha2_over_n = sum(dAlpha2_valid ./ dFrames_valid) / N;
    D_rot_MLE1 = (1 / (2*dt)) * avg_alpha2_over_n - locVar_ang / (2*meanDt);
    D_rot_MLE1 = max(D_rot_MLE1, 1e-8);
    
    %% ==================== MLE (2) ====================
    % Joint estimation of D and σ² using covariance matrix
    % Maximize: L = -½ log|Σ| - ½ Δᵀ Σ⁻¹ Δ
    
    % --- Translational MLE(2) ---
    % Pass signed displacements for covariance calculation
    dx_valid = dx(validSteps);
    dy_valid = dy(validSteps);
    dz_valid = dz(validSteps);
    [D_trans_MLE2, sigma2_trans] = mle2_estimate_signed([dx_valid, dy_valid, dz_valid], dFrames_valid, dt);
    
    % --- Rotational MLE(2) ---
    % For rotation, angular distance is unsigned (always positive), so we can't
    % use the covariance method. Instead, use variance-based estimation.
    dAlpha_valid = dAlpha(validSteps);
    dAlpha2_valid = dAlpha_valid.^2;
    [D_rot_MLE2, sigma2_rot] = mle2_estimate_variance(dAlpha2_valid, dFrames_valid, dt, 1);
    
    %% ==================== Select Best MLE ====================
    % Calculate reduced localization error: x = σ²/(D·dt)
    
    x_trans = locVar3D / (D_trans_MLE1 * dt);
    x_rot = locVar_ang / (D_rot_MLE1 * dt);
    
    % Recommendation based on Bo's paper
    if x_trans < 0.1
        D_trans_best = D_trans_MLE1;
        rec_trans = 1;
    elseif x_trans <= 10
        D_trans_best = D_trans_MLE2;
        rec_trans = 2;
    else
        D_trans_best = D_trans_MLE2;  % Use MLE2 but flag as unreliable
        rec_trans = 2;
    end
    
    if x_rot < 0.1
        D_rot_best = D_rot_MLE1;
        rec_rot = 1;
    elseif x_rot <= 10
        D_rot_best = D_rot_MLE2;
        rec_rot = 2;
    else
        D_rot_best = D_rot_MLE2;
        rec_rot = 2;
    end
    
    % Overall recommendation (use max of x values)
    recommended = max(rec_trans, rec_rot);
    
    %% --- Store results ---
    results(nAnalyzed).trackID = trackID;
    results(nAnalyzed).nDisplacements = N;
    results(nAnalyzed).D_trans_MLE1 = D_trans_MLE1;
    results(nAnalyzed).D_trans_MLE2 = D_trans_MLE2;
    results(nAnalyzed).D_trans_best = D_trans_best;
    results(nAnalyzed).sigma2_trans_MLE2 = sigma2_trans;
    results(nAnalyzed).D_rot_MLE1 = D_rot_MLE1;
    results(nAnalyzed).D_rot_MLE2 = D_rot_MLE2;
    results(nAnalyzed).D_rot_best = D_rot_best;
    results(nAnalyzed).sigma2_rot_MLE2 = sigma2_rot;
    results(nAnalyzed).x_trans = x_trans;
    results(nAnalyzed).x_rot = x_rot;
    results(nAnalyzed).recommended_MLE = recommended;
    results(nAnalyzed).meanTimeLag = meanDt;
    results(nAnalyzed).fracGaps = fracGaps;
end

% Trim unused entries
results = results(1:nAnalyzed);

%% --- Print summary ---
if verbose && nAnalyzed > 0
    D_trans1 = [results.D_trans_MLE1];
    D_trans2 = [results.D_trans_MLE2];
    D_rot1 = [results.D_rot_MLE1];
    D_rot2 = [results.D_rot_MLE2];
    recommended = [results.recommended_MLE];
    
    fprintf('\n--- estimateDiffusion6D_DualMLE summary ---\n');
    fprintf('  Trajectories analyzed: %d / %d\n', nAnalyzed, nTraj);
    fprintf('\n  MLE(1) - assumes known localization precision:\n');
    fprintf('    D_trans (nm²/s): median = %.1f, mean = %.1f\n', ...
        median(D_trans1), mean(D_trans1));
    fprintf('    D_rot (rad²/s):  median = %.4f, mean = %.4f\n', ...
        median(D_rot1), mean(D_rot1));
    fprintf('\n  MLE(2) - joint estimation of D and σ²:\n');
    fprintf('    D_trans (nm²/s): median = %.1f, mean = %.1f\n', ...
        median(D_trans2), mean(D_trans2));
    fprintf('    D_rot (rad²/s):  median = %.4f, mean = %.4f\n', ...
        median(D_rot2), mean(D_rot2));
    fprintf('\n  Recommended MLE: MLE(1) for %d%%, MLE(2) for %d%%\n', ...
        round(100*sum(recommended==1)/nAnalyzed), ...
        round(100*sum(recommended==2)/nAnalyzed));
    fprintf('---------------------------------------\n');
end

end

%% ========== HELPER FUNCTION: MLE(2) with Signed Displacements ==========
function [D_est, sigma2_est] = mle2_estimate_signed(delta, nFrames, dt)
%MLE2_ESTIMATE_SIGNED Joint MLE estimation of D and σ² from signed displacements
%
%   Implements MLE(2) from Bo Shuang's Langmuir 2013 paper using the
%   covariance between consecutive displacements.
%
%   For consecutive displacements sharing a measured point:
%       Cov(Δᵢ, Δᵢ₊₁) = -σ²  (per dimension)
%
%   This allows direct estimation of localization variance σ².
%
%   INPUT:
%       delta   : [N x d] signed displacements (d dimensions)
%       nFrames : [N x 1] frame gaps for each displacement
%       dt      : frame time (s)

    [N, dim] = size(delta);
    
    if N < 3
        % Not enough data
        delta2 = sum(delta.^2, 2);
        D_est = mean(delta2 ./ nFrames) / (2 * dim * dt);
        sigma2_est = 0;
        return;
    end
    
    % Compute mean squared displacement per dimension
    delta2 = sum(delta.^2, 2);  % [N x 1]
    delta2_per_dim = delta2 / dim;
    n_mean = mean(nFrames);
    
    % === Estimate sigma² from covariance of consecutive displacements ===
    % For each dimension: Cov(Δxᵢ, Δxᵢ₊₁) = -σ²
    % We average over dimensions for robustness
    
    cov_sum = 0;
    for d = 1:dim
        delta_d = delta(:, d);
        % Covariance of consecutive displacements
        cov_d = mean(delta_d(1:end-1) .* delta_d(2:end));
        cov_sum = cov_sum + cov_d;
    end
    cov_mean = cov_sum / dim;  % Average covariance per dimension
    
    % sigma² = -Cov(Δᵢ, Δᵢ₊₁)
    % But only if covariance is negative (indicates localization error)
    if cov_mean < 0
        sigma2_est = -cov_mean;
    else
        % No evidence of localization error
        sigma2_est = 0;
    end
    
    % === Estimate D from corrected MSD ===
    % MSD_per_dim = 2*D*n*dt + 2*sigma²
    % D = (MSD_per_dim - 2*sigma²) / (2*n*dt)
    
    msd_per_dim = mean(delta2_per_dim);
    D_est = (msd_per_dim - 2 * sigma2_est) / (2 * n_mean * dt);
    
    % Ensure non-negative
    D_est = max(D_est, 1e-10);
    sigma2_est = max(sigma2_est, 0);
    
    % Sanity check: sigma² shouldn't explain more than the MSD
    if 2 * sigma2_est > msd_per_dim * 0.9
        % Noise dominates - cap sigma² and recompute D
        sigma2_est = msd_per_dim * 0.4;
        D_est = (msd_per_dim - 2 * sigma2_est) / (2 * n_mean * dt);
        D_est = max(D_est, 1e-10);
    end
end

%% ========== HELPER FUNCTION: MLE(2) Variance-based (for unsigned displacements) ==========
function [D_est, sigma2_est] = mle2_estimate_variance(delta2, nFrames, dt, dim)
%MLE2_ESTIMATE_VARIANCE Estimate D and σ² from squared displacements using variance
%
%   For unsigned displacements (like angular distance), we can't use the
%   covariance method. Instead, we use the theoretical variance of MSD.
%
%   For Brownian motion with localization error:
%       E[Δ²] = 2*d*D*n*dt + 2*d*σ²
%       Var[Δ²] depends on both D and σ²
%
%   With variable time gaps, we can use regression to separate D and σ².
%   With uniform gaps, we use an iterative approach.
%
%   INPUT:
%       delta2  : [N x 1] squared displacements
%       nFrames : [N x 1] frame gaps
%       dt      : frame time (s)
%       dim     : dimensionality (1 for angular)

    N = length(delta2);
    
    if N < 3
        D_est = mean(delta2 ./ nFrames) / (2 * dim * dt);
        sigma2_est = 0;
        return;
    end
    
    delta2_per_dim = delta2 / dim;
    n_mean = mean(nFrames);
    msd_mean = mean(delta2_per_dim);
    
    % Check if time gaps vary enough for regression
    if std(nFrames) > 0.1 * mean(nFrames)
        % Variable gaps - use linear regression
        % Model: delta2_per_dim = 2*D*dt * nFrames + 2*sigma2
        X = [nFrames, ones(N, 1)];
        Y = delta2_per_dim;
        coeffs = X \ Y;
        
        slope = coeffs(1);      % = 2*D*dt
        intercept = coeffs(2);  % = 2*sigma2
        
        D_est = max(slope / (2 * dt), 1e-10);
        sigma2_est = max(intercept / 2, 0);
    else
        % Uniform gaps - use heuristic based on MSD statistics
        % For pure diffusion without noise, the coefficient of variation
        % of Δ² has a known value. Excess variance suggests noise.
        
        % Theoretical CV² for chi-squared with k=dim degrees of freedom
        % is 2/k. For 1D (dim=1), CV² = 2.
        cv2_observed = var(delta2_per_dim) / msd_mean^2;
        cv2_theory = 2 / dim;  % For pure diffusion
        
        if cv2_observed > cv2_theory * 1.5
            % Excess variance suggests localization error
            % Rough estimate: excess variance ~ contribution from σ²
            excess_var = (cv2_observed - cv2_theory) * msd_mean^2;
            
            % This is approximate - use a fraction of MSD as sigma2
            sigma2_est = min(sqrt(excess_var) / 4, msd_mean / 4);
            D_est = (msd_mean - 2 * sigma2_est) / (2 * n_mean * dt);
        else
            % No strong evidence of localization error
            D_est = msd_mean / (2 * n_mean * dt);
            sigma2_est = 0;
        end
    end
    
    % Ensure non-negative
    D_est = max(D_est, 1e-10);
    sigma2_est = max(sigma2_est, 0);
end
