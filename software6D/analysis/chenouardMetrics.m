function out = chenouardMetrics(trackedTraj, groundTruth, varargin)
%CHENOUARDMETRICS  Chenouard 2014 alpha/beta/JSC, faithful to the paper.
%
%   out = chenouardMetrics(trackedTraj, groundTruth)
%   out = chenouardMetrics(..., 'GateDist', 100)
%
%   Implements the alpha / beta / JSC metrics as defined in the main text of
%   Chenouard et al., Nature Methods 11:281-289 (2014), page 282, numbered
%   items (1)-(3) of "Quantitative performance measures":
%
%       alpha(X, Y) = 1 - d(X, Y) / d(X, empty)
%       beta(X, Y)  = (d(X, empty) - d(X, Y)) / (d(X, empty) + d(Y, empty))
%       JSC         = TP / (TP + FN + FP)
%
%   where:
%     - X = ground-truth tracks, Y = estimated tracks
%     - Y is extended with dummy tracks so every kX gets paired
%     - Munkres optimal assignment (MATLAB's matchpairs) pairs kX with kZ
%       drawn from the dummy-extended Y; minimal total distance
%     - d(kX, kZ) = sum over t in union(frames(kX), frames(kZ)) of the
%       gated Euclidean distance |kX(t) - kZ(t)|_{2,eps} = min(|.|_2, eps)
%       with a dummy point contributing the full gate at any t where one
%       of the tracks is missing a point
%     - d(X, empty) = sum over kX of |frames(kX)| * gate (maximum possible)
%     - d(Y, empty) = sum over spurious (unpaired-to-real-GT) est tracks
%       of |frames(kZ)| * gate
%     - TP, FN, FP are detection-level counts using the gated pairing:
%         TP = matched detection pairs within the Munkres pairing whose
%              per-frame distance <= gate
%         FN = GT detections not matched
%         FP = estimated detections not matched
%
%   Differences vs the earlier (deprecated) greedy implementation:
%     - Bipartite matching is Munkres-optimal (was greedy by score)
%     - Y is extended with explicit dummies (was implicit via length norm)
%     - Distance d(kX, kZ) follows the paper's truncated-distance form
%       (was rewritten as a per-frame similarity 1 - d/gate)
%
%   The deprecated greedy version is kept in chenouardMetrics_v1_greedy.m
%   for reproducibility of any earlier results.
%
%   REFERENCES
%     Chenouard, N. et al. Objective Comparison of Particle Tracking
%     Methods. Nat. Methods 11(3), 281-289 (2014). DOI 10.1038/nmeth.2808.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

p = inputParser;
addParameter(p, 'GateDist', 100, @isnumeric);
parse(p, varargin{:});
gateDist = p.Results.GateDist;

if isempty(trackedTraj) || isempty(groundTruth)
    out = emptyResult();
    return;
end

trkIDs = unique(trackedTraj(:, 2));
gtIDs  = unique(groundTruth(:, 2));
nTrk   = length(trkIDs);
nGT    = length(gtIDs);

if nGT == 0
    out = emptyResult();
    return;
end

trkMap = buildMap(trackedTraj, trkIDs);
gtMap  = buildMap(groundTruth,  gtIDs);

%% -- Per-pair distances d(kX, kZ) for kX in GT, kZ in real estimated --
% Uses union of frame indices; dummy fill contributes gate per missing frame.
realDistMat = zeros(nGT, nTrk);
realTPmat   = zeros(nGT, nTrk);     % # detection pairs <= gate (for JSC)

for i = 1:nGT
    fi = gtMap(i).frames;
    pi = gtMap(i).pos;
    for j = 1:nTrk
        fj = trkMap(j).frames;
        pj = trkMap(j).pos;
        [framesUnion, ~, ~] = union(fi, fj);
        d   = 0;
        tp  = 0;
        for k = 1:length(framesUnion)
            t   = framesUnion(k);
            iA  = find(fi == t, 1);
            iB  = find(fj == t, 1);
            if ~isempty(iA) && ~isempty(iB)
                dt = norm(pi(iA, :) - pj(iB, :));
                if dt <= gateDist
                    d  = d  + dt;
                    tp = tp + 1;
                else
                    d  = d  + gateDist;       % truncated to gate
                end
            else
                d = d + gateDist;             % dummy fill
            end
        end
        realDistMat(i, j) = d;
        realTPmat(i,   j) = tp;
    end
end

%% -- Dummy column: d(kX, dummy) = |frames(kX)| * gate -----------------
gtLengths  = arrayfun(@(s) length(s.frames), gtMap);
gtLengths  = gtLengths(:);                     % column
trkLengths = arrayfun(@(s) length(s.frames), trkMap);
trkLengths = trkLengths(:);
dummyCol   = gtLengths * gateDist;

% Cost matrix for Munkres: nGT rows x (nTrk + nGT) cols.
% Cols 1..nTrk are real estimated tracks; cols nTrk+1..end are nGT dummies.
% Each dummy column has cost dummyCol(i) only in row i; off-diagonal dummies
% have cost Inf (so they cannot be used by another GT).
dummyMat = inf(nGT, nGT);
for i = 1:nGT
    dummyMat(i, i) = dummyCol(i);
end
costMat = [realDistMat, dummyMat];

%% -- Munkres optimal assignment (MATLAB built-in) ---------------------
% costUnmatched large enough that every row gets paired (real or dummy).
M = matchpairs(costMat, sum(dummyCol) + 1);  % [rowIdx, colIdx] pairs
M = sortrows(M, 1);

% Tally d(X, Y), record real pairings
dXY            = 0;
matchedTrkIdx  = [];
realPairs      = zeros(0, 2);    % [iGT, iTrk] for real pairings
for k = 1:size(M, 1)
    iX = M(k, 1);
    iZ = M(k, 2);
    if iZ <= nTrk
        dXY            = dXY + realDistMat(iX, iZ);
        matchedTrkIdx  = [matchedTrkIdx; iZ]; %#ok<AGROW>
        realPairs(end+1, :) = [iX, iZ];        %#ok<AGROW>
    else
        % paired with this GT track's own dummy
        dXY = dXY + dummyCol(iX);
    end
end

%% -- d(X, empty), d(Y, empty), alpha, beta ----------------------------
dXempty   = sum(dummyCol);
% Spurious tracks: estimated tracks NOT paired with any GT
spuriousIdx = setdiff(1:nTrk, matchedTrkIdx(:)');
dYempty   = sum(trkLengths(spuriousIdx)) * gateDist;

alpha = 1 - dXY / max(dXempty, eps);
beta  = (dXempty - dXY) / max(dXempty + dYempty, eps);
beta  = max(0, beta);   % clip to >= 0 per paper convention

%% -- JSC: TP / (TP + FN + FP) on detection points ---------------------
TP = 0;
for k = 1:size(realPairs, 1)
    TP = TP + realTPmat(realPairs(k, 1), realPairs(k, 2));
end
% Total GT detection points
nGTDet  = sum(gtLengths);
nTrkDet = sum(trkLengths);
% FN: GT detections NOT counted as TP
FN = nGTDet - TP;
% FP: estimated detections NOT counted as TP
FP = nTrkDet - TP;
JSC = TP / max(TP + FN + FP, eps);

%% -- Pack ---------------------------------------------------------------
out = struct( ...
    'alpha',     alpha, ...
    'beta',      beta, ...
    'JSC',       JSC, ...
    'nMatched',  size(realPairs, 1), ...
    'nMissed',   nGT - size(realPairs, 1), ...
    'nSpurious', length(spuriousIdx), ...
    'nTPdet',    TP, ...
    'nFNdet',    FN, ...
    'nFPdet',    FP, ...
    'gateDist',  gateDist);
end


%% ========================================================================
%  HELPERS
%  ========================================================================
function m = buildMap(traj, ids)
n = length(ids);
m(n) = struct('id', [], 'frames', [], 'pos', []);
for k = 1:n
    rows = traj(traj(:, 2) == ids(k), :);
    rows = sortrows(rows, 1);
    m(k).id     = ids(k);
    m(k).frames = rows(:, 1);
    m(k).pos    = rows(:, 3:5);
end
end


function out = emptyResult()
out = struct('alpha', 0, 'beta', 0, 'JSC', 0, 'nMatched', 0, ...
    'nMissed', 0, 'nSpurious', 0, 'nTPdet', 0, 'nFNdet', 0, ...
    'nFPdet', 0, 'gateDist', NaN);
end
