function [trajFiltered, stats] = traj_filt_6D(trajData, varargin)
%TRAJ_FILT_6D Filter and interpolate 6D single-molecule trajectories.
%
%   Extended from traj_filt_3D_new3.m (Landes Research Group) to handle
%   orientation parameters (theta, phi, omega) with proper angular
%   interpolation.
%
%   Operations performed:
%     1. Split trajectories at gaps larger than connectsize.
%     2. Interpolate positions AND orientations for gaps <= connectsize.
%     3. Remove trajectories shorter than minTrajLen.
%
%   SYNOPSIS:
%       [trajFiltered, stats] = traj_filt_6D(trajData)
%       [trajFiltered, stats] = traj_filt_6D(trajData, 'Name', Value, ...)
%
%   INPUT:
%       trajData : [N x 9] matrix or table with columns:
%                  [frame, trackID, x_nm, y_nm, z_nm, theta, phi, omega, intensity]
%                  OR a cell array where each cell is a single trajectory
%                  matrix with the same column format.
%
%   NAME-VALUE PARAMETERS:
%       'ConnectSize'  : Max gap size (frames) to interpolate. Default: 3.
%       'MinTrajLen'   : Min trajectory length (frames) to keep. Default: 5.
%       'MaxSpatialGap': Max spatial jump (nm) allowed in interpolation.
%                        Gaps exceeding this are split instead. Default: Inf.
%       'MaxAngularGap': Max angular jump (rad) allowed. Default: pi/3.
%       'Verbose'      : Print stats. Default: true.
%
%   OUTPUT:
%       trajFiltered : [M x 9] matrix with filtered/interpolated trajectories.
%                      Track IDs are renumbered sequentially.
%       stats        : Struct with filtering statistics:
%                      .nInputTracks, .nOutputTracks, .nSplit, .nRemoved,
%                      .nInterpolatedFrames, .meanTrajLen, .medianTrajLen
%
%   NOTES:
%       - Spatial interpolation: linear in x, y, z.
%       - Theta interpolation: linear (safe for range [0, pi/2]).
%       - Phi interpolation: circular (handles wrapping on [-pi, pi]).
%       - Omega interpolation: linear (safe for range [0, 2*pi]).
%       - Intensity for interpolated frames: set to NaN.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'trajData');
addParameter(p, 'ConnectSize', 3, @isnumeric);
addParameter(p, 'MinTrajLen', 5, @isnumeric);
addParameter(p, 'MaxSpatialGap', Inf, @isnumeric);
addParameter(p, 'MaxAngularGap', pi/3, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
parse(p, trajData, varargin{:});

connectsize = p.Results.ConnectSize;
minTrajLen  = p.Results.MinTrajLen;
maxSpatGap  = p.Results.MaxSpatialGap;
maxAngGap   = p.Results.MaxAngularGap;
verbose     = p.Results.Verbose;

%% --- Convert input to standard format ---
% Expected columns: frame, trackID, x, y, z, theta, phi, omega, intensity

if istable(trajData)
    trajData = table2array(trajData);
end

if iscell(trajData)
    % Cell array of individual trajectories
    allTrajs = trajData;
else
    % Matrix: split by trackID (column 2)
    uniqueIDs = unique(trajData(:, 2));
    allTrajs = cell(length(uniqueIDs), 1);
    for i = 1:length(uniqueIDs)
        allTrajs{i} = trajData(trajData(:, 2) == uniqueIDs(i), :);
    end
end

nInputTracks = length(allTrajs);

%% --- Process each trajectory ---
outputTrajs = {};
nSplit = 0;
nRemoved = 0;
nInterpolated = 0;

for iTraj = 1:nInputTracks
    traj = allTrajs{iTraj};
    
    if isempty(traj) || size(traj, 1) < 2
        nRemoved = nRemoved + 1;
        continue;
    end
    
    % Sort by frame
    traj = sortrows(traj, 1);
    frames = traj(:, 1);
    
    % Find gaps between consecutive frames
    dFrames = diff(frames);
    
    %% --- Step 1: Split at large gaps ---
    splitPoints = find(dFrames > connectsize);
    
    if isempty(splitPoints)
        segments = {traj};
    else
        segments = cell(length(splitPoints) + 1, 1);
        prevIdx = 1;
        for iSplit = 1:length(splitPoints)
            segments{iSplit} = traj(prevIdx:splitPoints(iSplit), :);
            prevIdx = splitPoints(iSplit) + 1;
            nSplit = nSplit + 1;
        end
        segments{end} = traj(prevIdx:end, :);
    end
    
    %% --- Step 2: Interpolate small gaps within each segment ---
    for iSeg = 1:length(segments)
        seg = segments{iSeg};
        
        if isempty(seg) || size(seg, 1) < 1
            continue;
        end
        
        % Find remaining small gaps (1 < dFrame <= connectsize)
        segFrames = seg(:, 1);
        dF = diff(segFrames);
        gapIdx = find(dF > 1 & dF <= connectsize);
        
        if ~isempty(gapIdx)
            % Interpolate gaps (work backwards to preserve indices)
            for iGap = length(gapIdx):-1:1
                idx1 = gapIdx(iGap);      % Index of frame before gap
                idx2 = gapIdx(iGap) + 1;  % Index of frame after gap
                
                f1 = seg(idx1, 1);
                f2 = seg(idx2, 1);
                nGapFrames = f2 - f1 - 1;
                
                % Check spatial and angular jumps before interpolating
                spatJump = norm(seg(idx2, 3:5) - seg(idx1, 3:5));
                
                % Angular jump (with 180° symmetry)
                [angJump, ~, ~] = angularDistance(...
                    seg(idx1, 6), seg(idx1, 7), seg(idx1, 8), ...
                    seg(idx2, 6), seg(idx2, 7), seg(idx2, 8));
                
                if spatJump > maxSpatGap || angJump > maxAngGap
                    % Gap is too large — split here instead
                    nSplit = nSplit + 1;
                    continue;
                end
                
                % Generate interpolated rows
                interpRows = zeros(nGapFrames, size(seg, 2));
                
                for iF = 1:nGapFrames
                    t = iF / (nGapFrames + 1);  % Interpolation parameter [0, 1]
                    interpFrame = f1 + iF;
                    
                    interpRows(iF, 1) = interpFrame;        % frame
                    interpRows(iF, 2) = seg(idx1, 2);       % trackID (same)
                    
                    % Spatial: linear interpolation
                    interpRows(iF, 3) = seg(idx1, 3) + t * (seg(idx2, 3) - seg(idx1, 3));  % x
                    interpRows(iF, 4) = seg(idx1, 4) + t * (seg(idx2, 4) - seg(idx1, 4));  % y
                    interpRows(iF, 5) = seg(idx1, 5) + t * (seg(idx2, 5) - seg(idx1, 5));  % z
                    
                    % Theta: linear (safe for [0, pi/2])
                    interpRows(iF, 6) = seg(idx1, 6) + t * (seg(idx2, 6) - seg(idx1, 6));
                    
                    % Phi: circular interpolation on [-pi, pi]
                    interpRows(iF, 7) = circularInterp(seg(idx1, 7), seg(idx2, 7), t);
                    
                    % Omega: linear (safe for [0, 2*pi])
                    interpRows(iF, 8) = seg(idx1, 8) + t * (seg(idx2, 8) - seg(idx1, 8));
                    
                    % Intensity: NaN for interpolated frames
                    interpRows(iF, 9) = NaN;
                end
                
                nInterpolated = nInterpolated + nGapFrames;
                
                % Insert interpolated rows
                seg = [seg(1:idx1, :); interpRows; seg(idx2:end, :)];
            end
        end
        
        segments{iSeg} = seg;
    end
    
    %% --- Step 3: Filter short segments ---
    for iSeg = 1:length(segments)
        seg = segments{iSeg};
        if size(seg, 1) >= minTrajLen
            outputTrajs{end+1} = seg; %#ok<AGROW>
        else
            nRemoved = nRemoved + 1;
        end
    end
end

%% --- Reassemble and renumber track IDs ---
if isempty(outputTrajs)
    trajFiltered = zeros(0, 9);
else
    % Renumber track IDs sequentially
    for iTraj = 1:length(outputTrajs)
        outputTrajs{iTraj}(:, 2) = iTraj;
    end
    trajFiltered = vertcat(outputTrajs{:});
end

%% --- Statistics ---
nOutputTracks = length(outputTrajs);
if nOutputTracks > 0
    trajLens = cellfun(@(x) size(x, 1), outputTrajs);
    meanLen = mean(trajLens);
    medianLen = median(trajLens);
else
    meanLen = 0;
    medianLen = 0;
end

stats = struct();
stats.nInputTracks       = nInputTracks;
stats.nOutputTracks      = nOutputTracks;
stats.nSplit             = nSplit;
stats.nRemoved           = nRemoved;
stats.nInterpolatedFrames = nInterpolated;
stats.meanTrajLen        = meanLen;
stats.medianTrajLen      = medianLen;

if verbose
    fprintf('--- traj_filt_6D summary ---\n');
    fprintf('  Input tracks:         %d\n', nInputTracks);
    fprintf('  Output tracks:        %d\n', nOutputTracks);
    fprintf('  Splits:               %d\n', nSplit);
    fprintf('  Removed (too short):  %d\n', nRemoved);
    fprintf('  Interpolated frames:  %d\n', nInterpolated);
    fprintf('  Mean traj length:     %.1f frames\n', meanLen);
    fprintf('  Median traj length:   %.1f frames\n', medianLen);
    fprintf('----------------------------\n');
end

end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function phiInterp = circularInterp(phi1, phi2, t)
%CIRCULARINTERP Interpolate angles on a circle, handling wrap-around.
%
%   Interpolates between phi1 and phi2 on the interval [-pi, pi],
%   taking the shortest arc.

    dPhi = phi2 - phi1;
    
    % Wrap difference to [-pi, pi] for shortest path
    dPhi = mod(dPhi + pi, 2*pi) - pi;
    
    phiInterp = phi1 + t * dPhi;
    
    % Wrap result to [-pi, pi]
    phiInterp = mod(phiInterp + pi, 2*pi) - pi;
end
