function [costMat, propagationScheme, kalmanFilterInfoTmp, nonlinkMarker, ...
    errFlag, orientPropagationScheme] = costMat6DSMOLMLink(movieInfo, kalmanFilterInfoTmp, ...
    costMatParam, nnDistFeatures, probDim, prevCostStruct, ...
    featLifetime, trackedFeatureIndx, iFrame, varargin)
%COSTMAT6DSMOLMLINK 6D-SMOLM frame-to-frame linking cost function.
%
%   WRAPPER approach: calls the original costMatRandomDirectedSwitchingMotionLink
%   to get a properly formatted cost matrix, then modifies costs by adding
%   orientation distance terms for 6D-SMOLM tracking.
%
%   =========================================================================
%   ORIENTATION KALMAN LEVELS (orientKalmanLevel parameter)
%   =========================================================================
%   
%   LEVEL 0: NO ORIENTATION COST (spatial only)
%       - Equivalent to standard u-track3D
%       - Use when orientation data is unreliable
%   
%   LEVEL 1: STATIC VARIANCE (original method)
%       - Prediction: μ' = μ (no change predicted)
%       - Cost: C = Δα² / maxAngDist²
%       - Simple normalized angular distance
%       - No physics-based weighting
%   
%   LEVEL 2: ADAPTIVE KALMAN (Brownian, LEARNS D_rot)
%       - Prediction: μ' = μ (Brownian assumption)
%       - Variance: S = 4·D_rot·Δt + 2·σ_meas²
%         (factor of 4 because angular distance involves two independent
%          tangent-plane components, each with variance 2·D_rot·Δt)
%       - D_rot is LEARNED from track history (like u-track3D spatial)
%       - Falls back to D_rot_prior for new tracks
%   
%   LEVEL 3: FULL ADAPTIVE KALMAN (3 motion models, LEARNS D_rot)
%       - Three competing models like spatial Kalman:
%           1. Forward rotation:  μ' = R(+ω)·μ
%           2. Backward rotation: μ' = R(-ω)·μ  
%           3. Brownian:          μ' = μ
%       - Angular velocity ω estimated from track history
%       - D_rot LEARNED from angular innovation variance
%       - Best model stored in orientPropagationScheme
%   
%   ADAPTIVE LEARNING (Levels 2 & 3):
%       After linking, angular innovations are stored in kalmanFilterInfoTmp.
%       D_rot is estimated from mean squared angular displacements:
%           E[Δα²] = 4·D_rot·Δt + 2·σ²
%           D_rot_estimated = (mean(Δα²) - 2·σ²) / (4·Δt)
%
%   See docs/tracking_theory.html for full details.
%   =========================================================================
%
%   COST FORMULA:
%       C_total = wSpatial * C_spatial + wOrient * C_orient + wOmega * C_omega
%   
%   ADDITIONAL costMatParam FIELDS (beyond standard u-track):
%     --- Basic orientation parameters ---
%     wSpatial         : Weight for spatial cost (default: 1.0).
%     wOrient          : Weight for orientation cost (default: 1.0).
%     wOmega           : Weight for wobble cost (default: 0.0).
%     maxAngularDist   : Max angular distance for Level 1 normalization (default: pi/2).
%     useOrientation   : Enable orientation cost (default: true).
%
%     --- Orientation Kalman parameters ---
%     orientKalmanLevel: 0=none, 1=static, 2=adaptive, 3=full (default: 2).
%     D_rot_prior      : Prior rotational diffusion coeff, rad²/s (default: 0.1).
%                        Used only when track has insufficient history.
%     orientCRLB       : Orientation measurement uncertainty, rad (default: 5°).
%     frameTime        : Time between frames, seconds (default: 0.05).
%     useSignalWeighting : Scale orientCRLB by 1/sqrt(signal) (default: false).
%     useOmegaWeighting  : Scale measurement noise by 1/gamma(omega) so that
%                          wobbling molecules get a larger angular search
%                          window. gamma is the order parameter derived from
%                          the wobble solid angle omega (default: true).
%     gammaFloor         : Minimum gamma value to avoid divergent noise for
%                          highly wobbling molecules (default: 0.15).
%     minHistoryForAdaptive : Min track length to use adaptive D_rot (default: 3).
%
%     --- Diagnostic/save parameters ---
%     saveCostMatrix   : Save cost matrix for visualization (default: false).
%     costMatSavePath  : Path to save cost matrix (default: pwd).
%     costMatSaveFrame : Which frame to save (default: 5).
%
%   KALMAN FILTER INFO EXTENSIONS:
%     This function adds the following fields to kalmanFilterInfoTmp:
%       .orient6D.mu          : [nFeatures x 3] unit vector orientations
%       .orient6D.omega       : [nFeatures x 3] angular velocity (axis-angle)
%       .orient6D.D_rot       : [nFeatures x 1] learned D_rot per track
%       .orient6D.noiseHistory: {nFeatures x 1} cell of angular innovations
%
%   OUTPUT:
%     costMat          : Combined cost matrix
%     propagationScheme: Spatial propagation scheme from u-track3D
%     kalmanFilterInfoTmp: Kalman filter info (with orientation extensions)
%     nonlinkMarker    : Value indicating forbidden links
%     errFlag          : Error flag
%     orientPropagationScheme: (Level 3 only) Which orientation model was best
%                              1=forward, 2=backward, 3=brownian
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Initialize outputs ---
orientPropagationScheme = [];

%% --- Extract 6D parameters before passing to original ---
% Basic weights
wSpatial = getFieldOrDefault(costMatParam, 'wSpatial', 1.0);
wOrient = getFieldOrDefault(costMatParam, 'wOrient', 1.0);
wOmega = getFieldOrDefault(costMatParam, 'wOmega', 0.0);
maxAngDist = getFieldOrDefault(costMatParam, 'maxAngularDist', pi/2);
useOrient = getFieldOrDefault(costMatParam, 'useOrientation', true);

% Orientation Kalman parameters
if isfield(costMatParam, 'orientKalmanLevel')
    orientKalmanLevel = costMatParam.orientKalmanLevel;
elseif isfield(costMatParam, 'useOrientKalman') && costMatParam.useOrientKalman
    orientKalmanLevel = 2;  % Legacy: useOrientKalman=true maps to Level 2
else
    orientKalmanLevel = 2;  % Default: adaptive
end
D_rot_prior = getFieldOrDefault(costMatParam, 'D_rot_prior', 0.1);  % rad²/s
orientCRLB = getFieldOrDefault(costMatParam, 'orientCRLB', 5 * pi/180);  % 5 degrees default
frameTime = getFieldOrDefault(costMatParam, 'frameTime', 0.05);  % 50 ms default
useSignalWeighting = getFieldOrDefault(costMatParam, 'useSignalWeighting', false);
useOmegaWeighting = getFieldOrDefault(costMatParam, 'useOmegaWeighting', true);
gammaFloor = getFieldOrDefault(costMatParam, 'gammaFloor', 0.15);
minHistoryForAdaptive = getFieldOrDefault(costMatParam, 'minHistoryForAdaptive', 3);

% Diagnostic/save parameters
saveCostMat = getFieldOrDefault(costMatParam, 'saveCostMatrix', false);
costMatSavePath = getFieldOrDefault(costMatParam, 'costMatSavePath', pwd);
costMatSaveFrame = getFieldOrDefault(costMatParam, 'costMatSaveFrame', 5);

%% --- Call original u-track cost function ---
% Remove our custom fields so the original function doesn't choke on them
fieldsToRemove = {'wSpatial', 'wOrient', 'wOmega', 'maxAngularDist', ...
    'useOrientation', 'saveCostMatrix', 'costMatSavePath', 'costMatSaveFrame', ...
    'useOrientKalman', 'orientKalmanLevel', 'D_rot_prior', 'orientCRLB', ...
    'frameTime', 'useSignalWeighting', 'minHistoryForAdaptive'};
origParam = costMatParam;
for i = 1:length(fieldsToRemove)
    if isfield(origParam, fieldsToRemove{i})
        origParam = rmfield(origParam, fieldsToRemove{i});
    end
end

[costMat, propagationScheme, kalmanFilterInfoTmp, nonlinkMarker, errFlag] = ...
    costMatRandomDirectedSwitchingMotionLink(movieInfo, kalmanFilterInfoTmp, ...
    origParam, nnDistFeatures, probDim, prevCostStruct, ...
    featLifetime, trackedFeatureIndx, iFrame);

% Store spatial-only cost matrix immediately after base function call
costMatSpatial = costMat;
origNonlinkMarker = nonlinkMarker;

%% --- Initialize orientation Kalman state if not present ---
if ~isfield(kalmanFilterInfoTmp, 'orient6D')
    kalmanFilterInfoTmp.orient6D = struct();
end

%% --- Check if we should proceed with orientation modification ---
doOrientMod = true;
failReason = '';

% Level 0 means no orientation cost
if orientKalmanLevel == 0
    doOrientMod = false;
    failReason = 'orientKalmanLevel = 0 (disabled)';
    if wSpatial == 0
        warning('costMat6DSMOLMLink:NoCost', ...
            'Both wSpatial=0 and orientKalmanLevel=0. All linking costs will be zero!');
    end
end

% Early exit conditions
if doOrientMod && isempty(costMat)
    doOrientMod = false;
    failReason = 'costMat is empty';
end
if doOrientMod && (~isempty(errFlag) && any(errFlag ~= 0))
    doOrientMod = false;
    failReason = 'errFlag is non-zero';
end
if doOrientMod && (isscalar(costMat) || all(costMat(:) == nonlinkMarker))
    doOrientMod = false;
    failReason = 'costMat is scalar or all nonlinkMarker';
end
if doOrientMod && (~useOrient || wOrient == 0)
    doOrientMod = false;
    failReason = sprintf('useOrient=%d or wOrient=%.2f', useOrient, wOrient);
    if wSpatial == 0
        warning('costMat6DSMOLMLink:NoCost', ...
            'Both wSpatial=0 and wOrient=0. All linking costs will be zero!');
    end
end

%% --- Apply wSpatial weighting even if no orientation modification ---
if ~doOrientMod && ~isempty(costMat) && ~isscalar(costMat) && wSpatial ~= 1.0
    validMask = costMat ~= nonlinkMarker;
    if any(validMask(:))
        costMat(validMask) = wSpatial * costMat(validMask);
        if wSpatial > 0
            newMin = min(costMat(costMat ~= nonlinkMarker));
            if ~isempty(newMin) && isfinite(newMin)
                newNonlinkMarker = min(floor(newMin) - 5, -5);
                costMat(costMat == nonlinkMarker) = newNonlinkMarker;
                nonlinkMarker = newNonlinkMarker;
            end
        end
    end
end

%% --- Get frame info and validate ---
if doOrientMod
    [nRows, nCols] = size(costMat);
    
    if nargin < 9 || isnan(iFrame)
        doOrientMod = false;
        failReason = 'iFrame not provided or NaN';
    else
        iFrame = round(iFrame);
        nextFrame = iFrame + 1;
        if nextFrame > length(movieInfo)
            doOrientMod = false;
            failReason = sprintf('nextFrame %d > length(movieInfo) %d', nextFrame, length(movieInfo));
        end
    end
end

if doOrientMod
    if ~isfield(movieInfo, 'theta') || ~isfield(movieInfo, 'phi')
        doOrientMod = false;
        failReason = 'movieInfo missing theta or phi field';
    elseif isempty(movieInfo(iFrame).theta) || isempty(movieInfo(nextFrame).theta)
        doOrientMod = false;
        failReason = 'theta empty in current or next frame';
    end
end

%% --- Extract orientations ---
if doOrientMod
    nCurrDet = size(movieInfo(iFrame).theta, 1);
    nNextDet = size(movieInfo(nextFrame).theta, 1);
    nOrientRows = min(nRows, nCurrDet);
    
    if nNextDet == 0
        doOrientMod = false;
        failReason = 'nNextDet is 0';
    else
        nOrientCols = min(nCols, nNextDet);
        
        % Current frame orientations
        thetaCurr = movieInfo(iFrame).theta(1:nOrientRows, 1);
        phiCurr = movieInfo(iFrame).phi(1:nOrientRows, 1);
        muCurr = angles2mu(thetaCurr, phiCurr);
        
        % Next frame orientations
        thetaNext = movieInfo(nextFrame).theta(1:nOrientCols, 1);
        phiNext = movieInfo(nextFrame).phi(1:nOrientCols, 1);
        muNext = angles2mu(thetaNext, phiNext);
        
        % Wobble parameters
        if isfield(movieInfo, 'omega') && ~isempty(movieInfo(iFrame).omega)
            omegaCurr = movieInfo(iFrame).omega(1:nOrientRows, 1);
            omegaNext = movieInfo(nextFrame).omega(1:nOrientCols, 1);
        else
            omegaCurr = zeros(nOrientRows, 1);
            omegaNext = zeros(nOrientCols, 1);
        end
    end
end

%% --- Get track history for adaptive D_rot estimation ---
if doOrientMod && orientKalmanLevel >= 2
    % Get learned D_rot for each track from kalmanFilterInfoTmp
    % Or compute from track history using trackedFeatureIndx
    
    D_rot_perTrack = ones(nOrientRows, 1) * D_rot_prior;  % Default to prior
    angVel_perTrack = zeros(nOrientRows, 3);  % Angular velocity for Level 3
    
    if ~isempty(trackedFeatureIndx) && size(trackedFeatureIndx, 2) >= 2
        nFramesCols = size(trackedFeatureIndx, 2);
        for iTrack = 1:min(nOrientRows, size(trackedFeatureIndx, 1))
            % Iterate directly over frame columns to avoid detection index
            % ambiguity. trackedFeatureIndx(iTrack, frameCol) gives the
            % detection index at frame frameCol (0 = not present).
            muHistory = zeros(nFramesCols, 3);
            frameList = zeros(nFramesCols, 1);
            count = 0;

            for frameCol = 1:nFramesCols
                detIdx = trackedFeatureIndx(iTrack, frameCol);
                if detIdx > 0 && frameCol <= length(movieInfo) && ...
                   ~isempty(movieInfo(frameCol).theta) && ...
                   detIdx <= size(movieInfo(frameCol).theta, 1)

                    count = count + 1;
                    th = movieInfo(frameCol).theta(detIdx, 1);
                    ph = movieInfo(frameCol).phi(detIdx, 1);
                    muHistory(count, :) = angles2mu(th, ph);
                    frameList(count) = frameCol;
                end
            end

            % Trim to actual entries
            muHistory = muHistory(1:count, :);
            frameList = frameList(1:count);
            trackLen = count;

            if trackLen >= minHistoryForAdaptive
                % Compute angular displacements between consecutive detections
                angularDisplacements = [];
                for h = 2:trackLen
                    dAlpha = angularDistanceVectors(muHistory(h-1,:)', muHistory(h,:)');
                    angularDisplacements = [angularDisplacements; dAlpha]; %#ok<AGROW>
                end

                % Estimate D_rot from mean squared angular displacements.
                % E[dAngle^2] = 4*D_rot*dt + 2*sigma^2  (two tangent-plane
                % components each with variance 2*D_rot*dt, plus measurement
                % noise from both endpoints).
                % D_rot = (mean(dAngle^2) - 2*sigma^2) / (4*dt)
                if length(angularDisplacements) >= 2
                    meanSqDisp = mean(angularDisplacements.^2);
                    D_rot_est = (meanSqDisp - 2*orientCRLB^2) / (4*frameTime);
                    D_rot_est = max(D_rot_est, 0.001);  % Ensure positive
                    D_rot_est = min(D_rot_est, 10 * D_rot_prior);  % Cap at 10x prior
                    D_rot_perTrack(iTrack) = D_rot_est;

                    if saveCostMat && iFrame == costMatSaveFrame && iTrack <= 3
                        fprintf('  [costMat6DSMOLMLink] Track %d: D_rot_learned = %.4f rad²/s (from %d displacements)\n', ...
                            iTrack, D_rot_est, length(angularDisplacements));
                    end
                end

                % For Level 3: estimate angular velocity from last two frames
                if orientKalmanLevel == 3 && trackLen >= 2
                    mu_prev = muHistory(end-1, :)';
                    mu_curr = muHistory(end, :)';
                    [axis, angle] = estimateAngularVelocity(mu_prev, mu_curr);
                    angVel_perTrack(iTrack, :) = axis' * angle;  % axis-angle representation
                end
            end
        end
    end
end

%% --- Compute orientation cost based on level ---
if doOrientMod
    dt = frameTime;
    
    switch orientKalmanLevel
        case 1
            % =================================================================
            % LEVEL 1: STATIC VARIANCE (original method)
            % =================================================================
            [dAngle, dOmegaVal] = computeAngularDistanceMatrix(muCurr, muNext, omegaCurr, omegaNext);
            
            orientNorm = dAngle.^2 / max(maxAngDist^2, eps);
            omegaNorm = dOmegaVal.^2 / (2*pi)^2;
            
            if saveCostMat && iFrame == costMatSaveFrame
                fprintf('  [costMat6DSMOLMLink] LEVEL 1 (Static): maxAngDist = %.1f°\n', maxAngDist * 180/pi);
            end
            
        case 2
            % =================================================================
            % LEVEL 2: ADAPTIVE KALMAN (Brownian, LEARNS D_rot)
            % =================================================================
            [dAngle, dOmegaVal] = computeAngularDistanceMatrix(muCurr, muNext, omegaCurr, omegaNext);

            % Measurement variance — per-detection, modulated by omega.
            % When useOmegaWeighting is true, each detection's measurement
            % noise is scaled by 1/gamma(omega)^2: wobbling molecules have
            % larger effective angular uncertainty because the camera
            % integrates over many orientations within the wobble cone.
            [R_curr, R_next] = computeMeasurementVariance(orientCRLB, ...
                useSignalWeighting, movieInfo, iFrame, nextFrame, ...
                nOrientRows, nOrientCols, useOmegaWeighting, ...
                omegaCurr, omegaNext, gammaFloor);

            % Process variance using LEARNED D_rot per track.
            % Angular distance measures the total rotation angle on S^2,
            % involving two independent tangent-plane components each with
            % variance 2*D_rot*dt, so E[dAngle^2] = 4*D_rot*dt.
            P_predicted = zeros(nOrientRows, nOrientCols);
            for i = 1:nOrientRows
                P_predicted(i, :) = 4 * D_rot_perTrack(i) * dt;
            end

            % Total variance (innovation covariance).
            % S = P_process + R_curr(i) + R_next(j), where each R depends
            % on that detection's omega via the order parameter gamma.
            S_total = P_predicted + R_curr + R_next;

            % Normalized cost
            orientNorm = dAngle.^2 ./ max(S_total, eps);

            % Wobble cost (use mean D_rot)
            S_wobble = 4 * mean(D_rot_perTrack) * dt + 2 * (0.1)^2;
            omegaNorm = dOmegaVal.^2 ./ max(S_wobble, eps);

            if saveCostMat && iFrame == costMatSaveFrame
                fprintf('  [costMat6DSMOLMLink] LEVEL 2 (Adaptive Brownian):\n');
                fprintf('    D_rot: mean=%.4f, range=[%.4f, %.4f] rad²/s\n', ...
                    mean(D_rot_perTrack), min(D_rot_perTrack), max(D_rot_perTrack));
                fprintf('    S_total: mean=%.6f rad² (sqrt=%.2f°)\n', ...
                    mean(S_total(:)), sqrt(mean(S_total(:)))*180/pi);
                if useOmegaWeighting
                    gammaCurr = omegaToGamma(omegaCurr, gammaFloor);
                    fprintf('    Omega weighting ON: gamma range=[%.3f, %.3f]\n', ...
                        min(gammaCurr), max(gammaCurr));
                end
            end
            
        case 3
            % =================================================================
            % LEVEL 3: FULL ADAPTIVE KALMAN (3 motion models)
            % =================================================================

            % Measurement variance — omega-modulated per detection
            [R_curr, R_next] = computeMeasurementVariance(orientCRLB, ...
                useSignalWeighting, movieInfo, iFrame, nextFrame, ...
                nOrientRows, nOrientCols, useOmegaWeighting, ...
                omegaCurr, omegaNext, gammaFloor);

            % Process variance using LEARNED D_rot per track.
            % Two tangent-plane components: E[dAngle^2] = 4*D_rot*dt.
            P_predicted = zeros(nOrientRows, nOrientCols);
            for i = 1:nOrientRows
                P_predicted(i, :) = 4 * D_rot_perTrack(i) * dt;
            end
            S_total = P_predicted + R_curr + R_next;
            
            % Initialize cost matrices for each scheme
            dAngle_forward = zeros(nOrientRows, nOrientCols);
            dAngle_backward = zeros(nOrientRows, nOrientCols);
            dAngle_brownian = zeros(nOrientRows, nOrientCols);
            
            % For each track, predict orientation under each model
            for i = 1:nOrientRows
                mu_curr_i = muCurr(i, :)';
                angVel_i = angVel_perTrack(i, :)';
                omega_angle = norm(angVel_i);
                
                if omega_angle > 1e-6
                    omega_axis = angVel_i / omega_angle;
                    
                    % Predict under each model
                    mu_forward = rotateVector(mu_curr_i, omega_axis, +omega_angle);
                    mu_backward = rotateVector(mu_curr_i, omega_axis, -omega_angle);
                else
                    % No angular velocity - all models same as Brownian
                    mu_forward = mu_curr_i;
                    mu_backward = mu_curr_i;
                end
                mu_brownian = mu_curr_i;
                
                % Compute angular distance to all next-frame detections
                for j = 1:nOrientCols
                    mu_next_j = muNext(j, :)';
                    dAngle_forward(i,j) = angularDistanceVectors(mu_forward, mu_next_j);
                    dAngle_backward(i,j) = angularDistanceVectors(mu_backward, mu_next_j);
                    dAngle_brownian(i,j) = angularDistanceVectors(mu_brownian, mu_next_j);
                end
            end
            
            % Compute costs for each scheme
            cost_forward = dAngle_forward.^2 ./ max(S_total, eps);
            cost_backward = dAngle_backward.^2 ./ max(S_total, eps);
            cost_brownian = dAngle_brownian.^2 ./ max(S_total, eps);
            
            % Take minimum
            costStack = cat(3, cost_forward, cost_backward, cost_brownian);
            [orientNorm, orientPropagationScheme] = min(costStack, [], 3);
            
            % Wobble cost
            [~, dOmegaVal] = computeAngularDistanceMatrix(muCurr, muNext, omegaCurr, omegaNext);
            S_wobble = 4 * mean(D_rot_perTrack) * dt + 2 * (0.1)^2;
            omegaNorm = dOmegaVal.^2 ./ max(S_wobble, eps);
            dAngle = dAngle_brownian;  % For reporting
            
            if saveCostMat && iFrame == costMatSaveFrame
                fprintf('  [costMat6DSMOLMLink] LEVEL 3 (Full Adaptive 3-model):\n');
                fprintf('    D_rot: mean=%.4f, range=[%.4f, %.4f] rad²/s\n', ...
                    mean(D_rot_perTrack), min(D_rot_perTrack), max(D_rot_perTrack));
                nFwd = sum(orientPropagationScheme(:) == 1);
                nBwd = sum(orientPropagationScheme(:) == 2);
                nBrn = sum(orientPropagationScheme(:) == 3);
                fprintf('    Schemes: Forward=%d, Backward=%d, Brownian=%d\n', nFwd, nBwd, nBrn);
            end
            
        otherwise
            error('Invalid orientKalmanLevel: %d. Must be 0, 1, 2, or 3.', orientKalmanLevel);
    end
    
    %% --- Modify cost matrix ---
    validMask = costMat ~= nonlinkMarker;
    orientBlock = false(nRows, nCols);
    orientBlock(1:nOrientRows, 1:nOrientCols) = true;
    modifyMask = validMask & orientBlock;
    
    if any(modifyMask(:))
        origCosts = costMat(modifyMask);
        
        orientAddFull = zeros(nRows, nCols);
        orientAddBlock = wOrient * orientNorm + wOmega * omegaNorm;
        orientAddFull(1:nOrientRows, 1:nOrientCols) = orientAddBlock;
        
        costMat(modifyMask) = wSpatial * origCosts + orientAddFull(modifyMask);
        
        newMin = min(costMat(costMat ~= nonlinkMarker));
        if ~isempty(newMin)
            newNonlinkMarker = min(floor(newMin) - 5, -5);
            costMat(costMat == nonlinkMarker) = newNonlinkMarker;
            nonlinkMarker = newNonlinkMarker;
        end
    end
    
    %% --- Store orientation state for next iteration ---
    kalmanFilterInfoTmp.orient6D.mu = muCurr;
    kalmanFilterInfoTmp.orient6D.D_rot = D_rot_perTrack;
    if orientKalmanLevel == 3
        kalmanFilterInfoTmp.orient6D.angVel = angVel_perTrack;
    end
    
    %% --- Diagnostic output ---
    if saveCostMat && iFrame == costMatSaveFrame
        fprintf('  [costMat6DSMOLMLink] Frame %d: dAngle range = [%.1f, %.1f]°\n', ...
            iFrame, min(dAngle(:))*180/pi, max(dAngle(:))*180/pi);
        fprintf('  [costMat6DSMOLMLink] Frame %d: orientNorm range = [%.4f, %.4f]\n', ...
            iFrame, min(orientNorm(:)), max(orientNorm(:)));
    end
end

%% --- Save cost matrix for visualization ---
if saveCostMat && ~isempty(costMat) && ~isscalar(costMat)
    if nargin >= 9 && ~isnan(iFrame)
        currentFrame = round(iFrame);
    else
        currentFrame = -1;
    end
    
    if currentFrame == costMatSaveFrame
        levelNames = {'None', 'Static', 'Adaptive', 'Full3Model'};
        if doOrientMod
            fprintf('  [costMat6DSMOLMLink] Frame %d: Orientation cost APPLIED (Level %d: %s)\n', ...
                currentFrame, orientKalmanLevel, levelNames{orientKalmanLevel+1});
        else
            fprintf('  [costMat6DSMOLMLink] Frame %d: Orientation cost NOT applied. Reason: %s\n', ...
                currentFrame, failReason);
        end
        
        % Save data
        costMatFinalViz = costMat;
        costMatFinalViz(costMatFinalViz == nonlinkMarker) = 0;
        costMatSpatialViz = costMatSpatial;
        costMatSpatialViz(costMatSpatialViz == origNonlinkMarker) = 0;
        
        [nRowsSave, nColsSave] = size(costMat);
        nextFrameIdx = currentFrame + 1;
        if nextFrameIdx <= length(movieInfo)
            nDetSave = size(movieInfo(nextFrameIdx).xCoord, 1);
        else
            nDetSave = nColsSave;
        end
        
        costMatDataSave = struct();
        costMatDataSave.costMatFinal = sparse(costMatFinalViz);
        costMatDataSave.costMatSpatial = sparse(costMatSpatialViz);
        costMatDataSave.iFrame = currentFrame;
        costMatDataSave.nTracks = nRowsSave;
        costMatDataSave.nDetections = nDetSave;
        costMatDataSave.parameters = struct(...
            'wSpatial', wSpatial, ...
            'wOrient', wOrient, ...
            'wOmega', wOmega, ...
            'maxAngularDist', maxAngDist, ...
            'orientKalmanLevel', orientKalmanLevel, ...
            'D_rot_prior', D_rot_prior, ...
            'orientCRLB', orientCRLB, ...
            'frameTime', frameTime, ...
            'useSignalWeighting', useSignalWeighting);
        if doOrientMod && orientKalmanLevel >= 2
            costMatDataSave.D_rot_learned = D_rot_perTrack;
        end
        if orientKalmanLevel == 3 && ~isempty(orientPropagationScheme)
            costMatDataSave.orientPropagationScheme = orientPropagationScheme;
        end
        
        saveFile = fullfile(costMatSavePath, 'costMatrix_linking.mat');
        save(saveFile, '-struct', 'costMatDataSave');
        fprintf('  [costMat6DSMOLMLink] Cost matrix saved: %s\n', saveFile);
    end
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function val = getFieldOrDefault(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end

function mu = angles2mu(theta, phi)
%ANGLES2MU Convert spherical angles to unit vectors
    mu = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
end

function [dAngle, dOmega] = computeAngularDistanceMatrix(muCurr, muNext, omegaCurr, omegaNext)
%COMPUTEANGULARDISTANCEMATRIX Compute pairwise angular distances
    nRows = size(muCurr, 1);
    nCols = size(muNext, 1);
    
    dAngle = zeros(nRows, nCols);
    dOmega = zeros(nRows, nCols);
    
    for i = 1:nRows
        for j = 1:nCols
            dotProd = dot(muCurr(i,:), muNext(j,:));
            dAngle(i,j) = acos(min(abs(dotProd), 1));
            dOmega(i,j) = abs(omegaCurr(i) - omegaNext(j));
        end
    end
end

function dAlpha = angularDistanceVectors(mu1, mu2)
%ANGULARDISTANCEVECTORS Angular distance between two unit vectors
    dotProd = dot(mu1, mu2);
    dAlpha = acos(min(abs(dotProd), 1));
end

function [R_curr, R_next] = computeMeasurementVariance(orientCRLB, ...
    useSignalWeighting, movieInfo, iFrame, nextFrame, nRows, nCols, ...
    useOmegaWeighting, omegaCurr, omegaNext, gammaFloor)
%COMPUTEMEASUREMENTVARIANCE Compute per-detection measurement variance.
%
%   Returns separate R_curr [nRows x nCols] and R_next [nRows x nCols]
%   matrices so that the total measurement contribution to innovation
%   variance is S = P_process + R_curr + R_next.
%
%   When useOmegaWeighting is true, each detection's base variance is
%   scaled by 1/gamma(omega)^2, where gamma is the order parameter:
%       gamma = 1 - 3*omega/(4*pi) + omega^2/(8*pi^2)
%   This increases measurement noise for wobbling molecules (low gamma)
%   and keeps it at baseline for fixed dipoles (gamma ~ 1).

    if nargin < 8, useOmegaWeighting = false; end
    if nargin < 9, omegaCurr = []; end
    if nargin < 10, omegaNext = []; end
    if nargin < 11, gammaFloor = 0.15; end

    % --- Base CRLB variance (optionally signal-weighted) ---
    if useSignalWeighting && isfield(movieInfo, 'amp') && ...
       ~isempty(movieInfo(iFrame).amp) && ~isempty(movieInfo(nextFrame).amp)
        ampCurr = movieInfo(iFrame).amp(1:nRows, 1);
        ampNext = movieInfo(nextFrame).amp(1:nCols, 1);

        refSignal = 1000;
        sigma_curr = orientCRLB * sqrt(refSignal ./ max(ampCurr, 1));
        sigma_next = orientCRLB * sqrt(refSignal ./ max(ampNext, 1));
    else
        sigma_curr = repmat(orientCRLB, nRows, 1);
        sigma_next = repmat(orientCRLB, nCols, 1);
    end

    % --- Omega weighting: scale sigma by 1/gamma ---
    % gamma ~ 1: fixed dipole   -> sigma unchanged
    % gamma ~ 0: large wobble   -> sigma inflated (orientation unreliable)
    if useOmegaWeighting && ~isempty(omegaCurr) && ~isempty(omegaNext) && ...
       (any(omegaCurr > 0) || any(omegaNext > 0))
        gammaCurr = omegaToGamma(omegaCurr, gammaFloor);
        gammaNext = omegaToGamma(omegaNext, gammaFloor);
        sigma_curr = sigma_curr ./ gammaCurr;
        sigma_next = sigma_next ./ gammaNext;
    end

    % --- Build [nRows x nCols] variance matrices ---
    R_curr = repmat(sigma_curr.^2, 1, nCols);   % row i contributes R_curr(i)
    R_next = repmat(sigma_next'.^2, nRows, 1);   % col j contributes R_next(j)
end

function gamma = omegaToGamma(omega, gammaFloor)
%OMEGATOGAMMA Convert wobble solid angle to order parameter gamma.
%   gamma = 1 - 3*omega/(4*pi) + omega^2/(8*pi^2)
%   gamma ~ 1: fixed dipole (omega ~ 0)
%   gamma ~ 0: isotropic rotation (omega ~ 2*pi)
%   Clamped to [gammaFloor, 1] to avoid divergent variance.
    if nargin < 2, gammaFloor = 0.15; end
    gamma = 1 - 3*omega./(4*pi) + omega.^2./(8*pi^2);
    gamma = max(gammaFloor, min(1, gamma));
end

function [axis, angle] = estimateAngularVelocity(mu_prev, mu_curr)
%ESTIMATEANGULARVELOCITY Estimate rotation axis and angle from two orientations
    crossProd = cross(mu_prev, mu_curr);
    crossNorm = norm(crossProd);
    
    if crossNorm < 1e-10
        if abs(mu_curr(1)) < 0.9
            axis = cross(mu_curr, [1; 0; 0]);
        else
            axis = cross(mu_curr, [0; 1; 0]);
        end
        axis = axis / norm(axis);
        
        if dot(mu_prev, mu_curr) > 0
            angle = 0;
        else
            angle = pi;
        end
    else
        axis = crossProd / crossNorm;
        dotProd = dot(mu_prev, mu_curr);
        angle = acos(min(max(dotProd, -1), 1));
    end
    
    if angle > pi/2
        angle = pi - angle;
        axis = -axis;
    end
end

function mu_rot = rotateVector(mu, axis, angle)
%ROTATEVECTOR Rotate vector using Rodrigues' formula.
%   Uses the Höfling & Straube (2025) sign convention to match the
%   simulation in generate6DTrajectories.m:
%     v_rot = v*cos(a) - (k x v)*sin(a) + k*(k.v)*(1-cos(a))
    if abs(angle) < 1e-10
        mu_rot = mu;
        return;
    end

    mu_rot = mu * cos(angle) - ...
             cross(axis, mu) * sin(angle) + ...
             axis * dot(axis, mu) * (1 - cos(angle));
    mu_rot = mu_rot / norm(mu_rot);
end
