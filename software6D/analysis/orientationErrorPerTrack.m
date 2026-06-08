function out = orientationErrorPerTrack(trackedTraj, groundTruth, varargin)
%ORIENTATIONERRORPERTRACK Angular error on correctly-linked frames only.
%
%   out = orientationErrorPerTrack(trackedTraj, groundTruth)
%   out = orientationErrorPerTrack(..., 'DistThreshold', 50)
%
%   For each tracked trajectory, finds the GT particle it best matches
%   (most-frequent-spatial-match across frames), then computes the
%   angular distance acos(|mu_tracked . mu_GT|) on the subset of frames
%   where the tracked point is within DistThreshold of the matched GT
%   point. Aggregates across all such correctly-linked frames.
%
%   This is the metric the SMLM-tracking literature is missing: it
%   distinguishes "did 6D help linking" (validateTracking6D / Chenouard)
%   from "did 6D preserve orientation fidelity within the tracks it kept."
%
%   INPUT
%     trackedTraj : [N x 9] tracked matrix
%                   [frame, trackID, x, y, z, theta_rad, phi_rad, omega, intensity]
%     groundTruth : [M x 9] GT matrix, same format
%
%   NAME-VALUE
%     'DistThreshold' : Max spatial distance (nm) for a frame to count
%                       as "correctly linked". Default: 50.
%
%   OUTPUT (struct)
%     .meanErrDeg     : Mean angular error (deg) across correctly-linked frames
%     .medianErrDeg   : Median
%     .p95ErrDeg      : 95th percentile
%     .nLinkedFrames  : Total # frames contributing
%     .nTracks        : # tracks that had at least one correctly-linked frame
%     .perFrameErrDeg : [nLinkedFrames x 1] raw distribution (for histograms)
%     .perTrackMeanDeg: [nTracks x 1] mean error per matched track
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

p = inputParser;
addParameter(p, 'DistThreshold', 50, @isnumeric);
parse(p, varargin{:});
gate = p.Results.DistThreshold;

%% Empty-input guard
if isempty(trackedTraj) || isempty(groundTruth)
    out = emptyResult();
    return;
end

%% Index by trackID for fast lookup
trkIDs = unique(trackedTraj(:,2));
gtIDs  = unique(groundTruth(:,2));

%% For each tracked trajectory, find best-matching GT and accumulate errors
allErrors    = [];
perTrackMean = [];
nTracks      = 0;

for iTrk = 1:length(trkIDs)
    trkID = trkIDs(iTrk);
    trk   = trackedTraj(trackedTraj(:,2) == trkID, :);
    if isempty(trk), continue; end

    % Vote: at each tracked frame, which GT particle is the nearest within gate?
    votes = zeros(length(gtIDs), 1);
    for r = 1:size(trk,1)
        f = trk(r,1);
        gtInFrame = groundTruth(groundTruth(:,1) == f, :);
        if isempty(gtInFrame), continue; end
        d = sqrt(sum((gtInFrame(:,3:5) - trk(r,3:5)).^2, 2));
        [dmin, jmin] = min(d);
        if dmin <= gate
            gtMatchID = gtInFrame(jmin, 2);
            kk = find(gtIDs == gtMatchID, 1);
            if ~isempty(kk), votes(kk) = votes(kk) + 1; end
        end
    end
    if ~any(votes), continue; end
    [~, kBest] = max(votes);
    bestGTID = gtIDs(kBest);

    % Accumulate angular errors on frames matched to that GT particle
    gtOfMatched = groundTruth(groundTruth(:,2) == bestGTID, :);
    trkErrs = [];
    for r = 1:size(trk,1)
        f = trk(r,1);
        gtRow = gtOfMatched(gtOfMatched(:,1) == f, :);
        if isempty(gtRow), continue; end
        d = sqrt(sum((gtRow(1,3:5) - trk(r,3:5)).^2));
        if d > gate, continue; end
        % Both as unit vectors
        mu_t = angles2mu(trk(r,6), trk(r,7));
        mu_g = angles2mu(gtRow(1,6), gtRow(1,7));
        dotProd = dot(mu_t, mu_g);
        cosAng  = min(1, max(0, abs(dotProd)));
        trkErrs(end+1, 1) = acos(cosAng) * 180/pi; %#ok<AGROW>
    end
    if isempty(trkErrs), continue; end
    nTracks = nTracks + 1;
    allErrors    = [allErrors;    trkErrs];        %#ok<AGROW>
    perTrackMean = [perTrackMean; mean(trkErrs)];  %#ok<AGROW>
end

if isempty(allErrors)
    out = emptyResult();
    return;
end

out = struct( ...
    'meanErrDeg',      mean(allErrors), ...
    'medianErrDeg',    median(allErrors), ...
    'p95ErrDeg',       prctile(allErrors, 95), ...
    'nLinkedFrames',   length(allErrors), ...
    'nTracks',         nTracks, ...
    'perFrameErrDeg',  allErrors, ...
    'perTrackMeanDeg', perTrackMean, ...
    'gateDist',        gate);
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function mu = angles2mu(theta, phi)
    mu = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
end

function out = emptyResult()
    out = struct('meanErrDeg', NaN, 'medianErrDeg', NaN, 'p95ErrDeg', NaN, ...
        'nLinkedFrames', 0, 'nTracks', 0, ...
        'perFrameErrDeg', [], 'perTrackMeanDeg', [], 'gateDist', NaN);
end
