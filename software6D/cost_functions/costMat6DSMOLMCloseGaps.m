function [costMat, nonlinkMarker, indxMerge, numMerge, indxSplit, numSplit, ...
    errFlag] = costMat6DSMOLMCloseGaps(trackedFeatInfo, trackedFeatIndx, ...
    trackStartTime, trackEndTime, costMatParam, gapCloseParam, ...
    kalmanFilterInfo, nnDistLinkedFeat, probDim, movieInfo)
%COSTMAT6DSMOLMCLOSEGAPS 6D-SMOLM gap closing cost function.
%
%   WRAPPER approach: calls the original costMatRandomDirectedSwitchingMotionCloseGaps
%   to get a properly formatted cost matrix, then modifies costs by adding
%   orientation distance terms for 6D-SMOLM tracking.
%
%   INPUT (same as costMatRandomDirectedSwitchingMotionCloseGaps):
%     trackedFeatInfo  : Positions/amplitudes from linkFeaturesKalman
%     trackedFeatIndx  : Connectivity matrix of features between frames
%     trackStartTime   : Starting time of all tracks
%     trackEndTime     : Ending time of all tracks
%     costMatParam     : Cost matrix parameters (see below for 6D extensions)
%     gapCloseParam    : Gap closing parameters (.timeWindow, .mergeSplit, etc.)
%     kalmanFilterInfo : Kalman filter state info per frame
%     nnDistLinkedFeat : Nearest neighbor distances of linked features
%     probDim          : Problem dimensionality (2 or 3)
%     movieInfo        : movieInfo struct with orientation fields
%
%   ADDITIONAL costMatParam FIELDS (beyond standard u-track):
%     wSpatial       : Weight for spatial cost (default: 1.0).
%     wOrient        : Weight for orientation cost (default: 0.3).
%     wOmega         : Weight for wobble cost (default: 0.0).
%     maxAngularDist : Max angular distance threshold (default: pi/3).
%     useOrientation : Enable orientation cost (default: true).
%     useBoScaling   : Normalize by sqrt(timeGap) (default: true).
%     saveCostMatrix : Save cost matrix for visualization (default: false).
%     costMatSavePath: Path to save cost matrix (default: pwd).
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Extract 6D parameters ---
if isfield(costMatParam, 'wSpatial')
    wSpatial = costMatParam.wSpatial;
else
    wSpatial = 1.0;
end
if isfield(costMatParam, 'wOrient')
    wOrient = costMatParam.wOrient;
else
    wOrient = 0.3;
end
if isfield(costMatParam, 'wOmega')
    wOmega = costMatParam.wOmega;
else
    wOmega = 0.0;
end
if isfield(costMatParam, 'maxAngularDist')
    maxAngDist = costMatParam.maxAngularDist;
else
    maxAngDist = pi/3;
end
if isfield(costMatParam, 'useOrientation')
    useOrient = costMatParam.useOrientation;
else
    useOrient = true;
end
if isfield(costMatParam, 'useBoScaling')
    useBoScaling = costMatParam.useBoScaling;
else
    useBoScaling = true;
end
if isfield(costMatParam, 'saveCostMatrix')
    saveCostMat = costMatParam.saveCostMatrix;
else
    saveCostMat = false;
end
if isfield(costMatParam, 'costMatSavePath')
    costMatSavePath = costMatParam.costMatSavePath;
else
    costMatSavePath = pwd;
end
% Cost-scale normalization (matches costMat6DSMOLMLink.m).
if isfield(costMatParam, 'normalizeSpatial')
    normalizeSpatial = costMatParam.normalizeSpatial;
else
    normalizeSpatial = true;
end
% Quiet diagnostic prints from the base function call by default.
if isfield(costMatParam, 'verboseDiagnostic')
    verboseDiagnostic = costMatParam.verboseDiagnostic;
else
    verboseDiagnostic = false;
end

%% --- Call original u-track gap closing cost function ---
% Remove our custom fields so they don't confuse the base function
fieldsToRemove = {'wSpatial', 'wOrient', 'wOmega', 'maxAngularDist', ...
    'useOrientation', 'useBoScaling', 'saveCostMatrix', 'costMatSavePath', ...
    'normalizeSpatial', 'verboseDiagnostic'};
origParam = costMatParam;
for i = 1:length(fieldsToRemove)
    if isfield(origParam, fieldsToRemove{i})
        origParam = rmfield(origParam, fieldsToRemove{i});
    end
end

% Call with correct 10 arguments
[costMat, nonlinkMarker, indxMerge, numMerge, indxSplit, numSplit, errFlag] = ...
    costMatRandomDirectedSwitchingMotionCloseGaps(trackedFeatInfo, ...
    trackedFeatIndx, trackStartTime, trackEndTime, origParam, ...
    gapCloseParam, kalmanFilterInfo, nnDistLinkedFeat, probDim, movieInfo);

% Save spatial-only cost matrix for comparison
if saveCostMat && ~isempty(costMat)
    costMatSpatial = costMat;
end

% --- Diagnostic output (gated) ---
if verboseDiagnostic
    fprintf('  [costMat6DSMOLMCloseGaps] Base function returned:\n');
    fprintf('    costMat empty: %d, size: [%d x %d]\n', isempty(costMat), size(costMat, 1), size(costMat, 2));
    fprintf('    costMat sparse: %d\n', issparse(costMat));
    fprintf('    errFlag: %s\n', mat2str(errFlag));
    if ~isempty(costMat) && issparse(costMat)
        fprintf('    nonzero entries (potential links): %d\n', nnz(costMat));
    end
    fprintf('    nTracks (ends): %d, (starts): %d\n', length(trackEndTime), length(trackStartTime));
end

%% --- Early exit conditions ---
if isempty(costMat) || (~isempty(errFlag) && any(errFlag ~= 0))
    return;
end

if ~useOrient || wOrient == 0
    return;
end

if ~issparse(costMat)
    return;
end

%% --- Add orientation cost to valid entries ---
% costMat is sparse: nonzero entries are valid link costs.
% Rows = track ends (candidates for gap closing as "sources")
% Cols = track starts (candidates for gap closing as "sinks")
% Plus merge/split columns if enabled.

[nRows, nCols] = size(costMat);
[ii, jj, vv] = find(costMat);

if isempty(ii)
    return;
end

% Check that movieInfo has orientation fields
if ~isfield(movieInfo, 'theta') || ~isfield(movieInfo, 'phi')
    return;
end

% Number of tracks (gap closing operates on track ends -> track starts)
nTracks = length(trackEndTime);

% Normalize spatial cost to median=1 across valid links so wSpatial and
% wOrient operate on comparable scales (see header comment).
if normalizeSpatial
    spatialScale = median(vv(vv > 0));
    if isempty(spatialScale) || ~isfinite(spatialScale) || spatialScale <= eps
        spatialScale = 1;
    end
    vvScaled = vv / spatialScale;
else
    vvScaled = vv;
end

% For each valid link, compute the orientation cost
newVals = vvScaled;

for k = 1:length(ii)
    iEnd = ii(k);    % Track end index (row)
    jStart = jj(k);  % Track start index (column)
    
    % Only process gap-closing links (not merge/split)
    if iEnd > nTracks || jStart > nTracks
        continue;
    end
    
    % Get the last frame of the ending track and first frame of the starting track
    tEnd   = trackEndTime(iEnd);
    tStart = trackStartTime(jStart);
    
    % Validate frame indices
    if tEnd < 1 || tEnd > length(movieInfo) || tStart < 1 || tStart > length(movieInfo)
        continue;
    end
    
    % Time gap
    timeGap = tStart - tEnd;
    if timeGap <= 0
        continue;
    end
    
    % Get orientation at track end
    if isempty(movieInfo(tEnd).theta)
        continue;
    end
    
    % Find the detection index for the end of track iEnd
    % trackedFeatIndx has 1 row per track, columns cycle per frame
    % The detection index at tEnd is stored in the appropriate column
    endDetIdx = findDetectionIndex(trackedFeatIndx, iEnd, tEnd);
    startDetIdx = findDetectionIndex(trackedFeatIndx, jStart, tStart);
    
    if endDetIdx < 1 || endDetIdx > size(movieInfo(tEnd).theta, 1) || ...
       startDetIdx < 1 || startDetIdx > size(movieInfo(tStart).theta, 1)
        continue;
    end
    
    thetaEnd = movieInfo(tEnd).theta(endDetIdx, 1);
    phiEnd   = movieInfo(tEnd).phi(endDetIdx, 1);
    thetaStart = movieInfo(tStart).theta(startDetIdx, 1);
    phiStart   = movieInfo(tStart).phi(startDetIdx, 1);
    
    % Wobble
    omegaEnd = 0;
    omegaStart = 0;
    if isfield(movieInfo, 'omega')
        if ~isempty(movieInfo(tEnd).omega)
            omegaEnd = movieInfo(tEnd).omega(endDetIdx, 1);
        end
        if ~isempty(movieInfo(tStart).omega)
            omegaStart = movieInfo(tStart).omega(startDetIdx, 1);
        end
    end
    
    % Compute angular distance
    [dAngle, dOmega, ~] = angularDistance(thetaEnd, phiEnd, omegaEnd, ...
        thetaStart, phiStart, omegaStart);
    
    % Bo Shuang scaling: normalize by sqrt(timeGap)
    if useBoScaling && timeGap > 1
        dAngle = dAngle / sqrt(timeGap);
        dOmega = dOmega / sqrt(timeGap);
    end
    
    % Add orientation cost
    orientCost = wOrient * (dAngle^2 / max(maxAngDist^2, eps)) + ...
                 wOmega  * (dOmega^2 / (2*pi)^2);  % Wobble range is [0, 2π] steradians
    
    newVals(k) = wSpatial * vvScaled(k) + orientCost;
end

% Rebuild sparse matrix
costMat = sparse(ii, jj, newVals, nRows, nCols);

% Update nonlinkMarker
if nnz(costMat) > 0
    newMin = full(min(nonzeros(costMat)));
    newNonlinkMarker = min(floor(newMin) - 5, nonlinkMarker);
    if newNonlinkMarker ~= nonlinkMarker
        nonlinkMarker = newNonlinkMarker;
    end
end

%% --- Save cost matrix for visualization if requested ---
if saveCostMat
    costMatData = struct();
    costMatData.costMatFinal = costMat;
    costMatData.costMatSpatial = costMatSpatial;
    costMatData.nTracks = nTracks;
    costMatData.trackStartTime = trackStartTime;
    costMatData.trackEndTime = trackEndTime;
    costMatData.parameters = struct('wSpatial', wSpatial, 'wOrient', wOrient, ...
        'wOmega', wOmega, 'maxAngularDist', maxAngDist, 'useBoScaling', useBoScaling);
    
    saveFile = fullfile(costMatSavePath, 'costMatrix_gapClosing.mat');
    save(saveFile, '-struct', 'costMatData');
    fprintf('  [costMat6DSMOLMCloseGaps] Cost matrix saved to: %s\n', saveFile);
end

end

%% === Helper function ===
function detIdx = findDetectionIndex(trackedFeatIndx, trackIdx, frameNum)
%FINDDETECTIONINDEX Get the detection index for a track at a specific frame.
%   trackedFeatIndx has rows=tracks, and columns cycle through frames
%   with 8 columns per frame: [x, y, z, amp, dx, dy, dz, damp].
%   The detection index is the feature number in movieInfo(frameNum).

detIdx = 0;

if trackIdx > size(trackedFeatIndx, 1)
    return;
end

% trackedFeatIndx stores the detection index for each frame
% Column for frameNum is frameNum itself (1-indexed)
nCols = size(trackedFeatIndx, 2);

if frameNum > nCols || frameNum < 1
    return;
end

detIdx = trackedFeatIndx(trackIdx, frameNum);

% Handle 0 or NaN
if isnan(detIdx) || detIdx <= 0
    detIdx = 0;
end

end
