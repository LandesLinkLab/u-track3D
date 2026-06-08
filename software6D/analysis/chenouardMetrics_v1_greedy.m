function out = chenouardMetrics(trackedTraj, groundTruth, varargin)
%CHENOUARDMETRICS Particle-tracking metrics from Chenouard et al. 2014.
%
%   out = chenouardMetrics(trackedTraj, groundTruth)
%   out = chenouardMetrics(..., 'GateDist', 100, 'MaxLifetime', [])
%
%   Implements the alpha / beta / JSC metrics from the ISBI Particle
%   Tracking Challenge (Chenouard et al., Nat Methods 11:281, 2014).
%
%   INPUT
%     trackedTraj : [N x 9] tracked matrix, columns
%                   [frame, trackID, x, y, z, theta, phi, omega, intensity]
%     groundTruth : [M x 9] GT matrix, same format. Column 2 = true ID.
%
%   NAME-VALUE
%     'GateDist'    : Max position distance (nm) at which a tracked and GT
%                     point are considered matchable. Default: 100.
%     'MaxLifetime' : Optional cap on track lifetime for per-track scoring.
%                     Default: [] (use the longer of GT/est lifetime).
%
%   OUTPUT (struct)
%     .alpha        : Track-overlap score in [0,1], 1 = perfect.
%                     alpha = sum over matched pairs of (1 - d/gate) / Nmax
%     .beta         : alpha with spurious / missed track penalty in [0,1].
%                     beta = (alpha_sum - n_unmatched_est) / Nmax  (clipped 0)
%     .JSC          : Jaccard on per-frame positions across all tracks
%                     = |matched| / (|GT detections| + |est detections| - |matched|)
%     .nMatched     : Number of track pairs successfully matched
%     .nMissed      : GT tracks with no acceptable est track
%     .nSpurious    : Est tracks with no acceptable GT track
%     .perPairAlpha : [nMatched x 1] alpha contribution per matched pair
%
%   NOTES
%     - Pure-position metric; orientation columns are not used.
%     - Track-to-track matching is bipartite, solved greedily on the
%       length-normalized alpha score (suitable for hundreds of tracks).
%     - alpha and beta normalize by the sum of GT track lengths, so
%       reported values are bounded in [0, 1].
%
%   REFERENCES
%     Chenouard et al., Nature Methods 11:281, 2014.
%     "Objective comparison of particle tracking methods."
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

p = inputParser;
addParameter(p, 'GateDist',    100, @isnumeric);
addParameter(p, 'MaxLifetime', [],  @isnumeric);
parse(p, varargin{:});
gateDist = p.Results.GateDist;

%% --- Sanity ---
if isempty(trackedTraj)
    out = emptyResult();
    return;
end

trkIDs = unique(trackedTraj(:,2));
gtIDs  = unique(groundTruth(:,2));
nTrk   = length(trkIDs);
nGT    = length(gtIDs);

%% --- Build per-track frame -> position maps ---
trkMap = buildMap(trackedTraj, trkIDs);
gtMap  = buildMap(groundTruth, gtIDs);

%% --- Score every (estimated, GT) pair on length-normalized alpha ---
% alpha(i,j) = sum over shared frames of (1 - d_ij(t)/gate)+, divided by
% the LONGER of the two track lifetimes.
scoreMat = zeros(nTrk, nGT);
matchMat = zeros(nTrk, nGT);  % number of matched frames per pair

for i = 1:nTrk
    fi = trkMap(i).frames;
    pi = trkMap(i).pos;
    for j = 1:nGT
        fj = gtMap(j).frames;
        pj = gtMap(j).pos;
        [common, ia, ib] = intersect(fi, fj);
        if isempty(common), continue; end
        d = sqrt(sum((pi(ia,:) - pj(ib,:)).^2, 2));
        contrib = max(0, 1 - d/gateDist);  % zero where d > gate
        Lmax = max(length(fi), length(fj));
        scoreMat(i,j)  = sum(contrib) / max(Lmax, 1);
        matchMat(i,j)  = sum(d <= gateDist);
    end
end

%% --- Greedy bipartite matching on score (descending) ---
[srt, idx] = sort(scoreMat(:), 'descend');
estUsed = false(nTrk, 1);
gtUsed  = false(nGT, 1);
pairs   = zeros(0, 2);
pairScore = zeros(0,1);
for k = 1:length(srt)
    if srt(k) <= 0, break; end
    [iE, iG] = ind2sub([nTrk, nGT], idx(k));
    if estUsed(iE) || gtUsed(iG), continue; end
    estUsed(iE) = true;
    gtUsed(iG)  = true;
    pairs(end+1, :) = [iE, iG]; %#ok<AGROW>
    pairScore(end+1) = srt(k);  %#ok<AGROW>
end

nMatched  = size(pairs, 1);
nMissed   = nGT - nMatched;
nSpurious = nTrk - nMatched;

%% --- alpha and beta ---
% alpha is the mean per-track score weighted by GT track length so longer
% GT tracks count more. Normalized to [0,1] by dividing by sum of GT
% lengths. Per Chenouard's formulation alpha is a similarity, beta is
% alpha minus the spurious-track penalty.
gtLengths = cellfun(@length, {gtMap.frames})';
trkLengths = cellfun(@length, {trkMap.frames})';

% sum-of-contributions over all matched pairs (already length-weighted)
% = sum over pairs of pairScore * Lmax_pair
pairLmax = zeros(nMatched, 1);
for k = 1:nMatched
    iE = pairs(k,1); iG = pairs(k,2);
    pairLmax(k) = max(length(trkMap(iE).frames), length(gtMap(iG).frames));
end
alphaSum  = sum(pairScore(:) .* pairLmax);
normLen   = max(sum(gtLengths), 1);
alpha     = alphaSum / normLen;

% Spurious-track penalty: total length of un-matched estimated tracks,
% normalized by sum of GT lengths. Beta clipped at 0.
spuriousLen = sum(trkLengths(~estUsed));
beta = max(0, (alphaSum - spuriousLen) / normLen);

%% --- JSC: per-frame position Jaccard across all tracks ---
% Treat each (frame, GT_id) as a GT point; each (frame, est_id) as an est
% point. A point pair is "matched" if (a) the est track is matched to that
% GT id, and (b) the per-frame distance <= gate.
nMatchedFrames = 0;
for k = 1:nMatched
    iE = pairs(k,1); iG = pairs(k,2);
    nMatchedFrames = nMatchedFrames + matchMat(iE, iG);
end
nGTDetections  = sum(gtLengths);
nEstDetections = sum(trkLengths);
denom = nGTDetections + nEstDetections - nMatchedFrames;
JSC = nMatchedFrames / max(denom, 1);

%% --- Pack output ---
out = struct( ...
    'alpha',        alpha, ...
    'beta',         beta, ...
    'JSC',          JSC, ...
    'nMatched',     nMatched, ...
    'nMissed',      nMissed, ...
    'nSpurious',    nSpurious, ...
    'perPairAlpha', pairScore(:), ...
    'gateDist',     gateDist);
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function m = buildMap(traj, ids)
    n = length(ids);
    m(n) = struct('id', [], 'frames', [], 'pos', []);
    for k = 1:n
        rows = traj(traj(:,2) == ids(k), :);
        rows = sortrows(rows, 1);   % by frame
        m(k).id     = ids(k);
        m(k).frames = rows(:,1);
        m(k).pos    = rows(:, 3:5);
    end
end

function out = emptyResult()
    out = struct('alpha',0,'beta',0,'JSC',0,'nMatched',0,'nMissed',0, ...
        'nSpurious',0,'perPairAlpha',[],'gateDist',NaN);
end
