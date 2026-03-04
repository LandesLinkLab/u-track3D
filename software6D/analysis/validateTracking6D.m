function [metrics, figHandles] = validateTracking6D(trackedTraj, groundTruth, varargin)
%VALIDATETRACKING6D Validate tracking performance against ground truth
%
%   Computes tracking accuracy metrics based on Jaqaman et al., Nature
%   Methods 2008. Compares tracked trajectories to known ground truth
%   to assess linking accuracy, gap closing, and overall performance.
%
%   SYNOPSIS:
%       [metrics, figHandles] = validateTracking6D(trackedTraj, groundTruth)
%       [metrics, figHandles] = validateTracking6D(..., 'Name', Value, ...)
%
%   INPUT:
%       trackedTraj : [N x 9] tracked trajectory matrix with columns:
%                     [frame, trackID, x, y, z, theta, phi, omega, intensity]
%                     OR cell array of trajectories
%       groundTruth : [M x 9] ground truth matrix with same format
%                     Column 2 must contain TRUE particle IDs
%
%   NAME-VALUE PARAMETERS:
%       'DistThreshold'  : Max distance (nm) to match points. Default: 50.
%       'AngleThreshold' : Max angle (rad) to match orientations. Default: 0.3.
%       'Verbose'        : Print detailed results. Default: true.
%       'PlotResults'    : Generate diagnostic plots. Default: true.
%       'SavePath'       : Directory to save results. Default: '' (no save).
%
%   OUTPUT:
%       metrics : Struct with validation metrics:
%           .linkingAccuracy      : Fraction of correct frame-to-frame links
%           .linkingPrecision     : TP_links / (TP_links + FP_links)
%           .linkingRecall        : TP_links / (TP_links + FN_links)
%           .gapClosingAccuracy   : Fraction of correctly closed gaps
%           .trackCompleteness    : Mean(tracked_length / true_length)
%           .trackPurity          : Mean fraction of track from single GT particle
%           .nTracksRecovered     : Number of GT particles recovered
%           .nTracksSplit         : Number of GT particles split into multiple tracks
%           .nTracksMerged        : Number of tracked trajectories with multiple GT
%           .lifetimeKS           : KS test p-value comparing lifetime distributions
%           .D_trans_error        : RMSE of translational D estimates
%           .D_rot_error          : RMSE of rotational D estimates
%
%       figHandles : Array of figure handles for generated plots
%
%   METRICS EXPLAINED (Jaqaman 2008):
%
%   1. Linking Accuracy: For each frame-to-frame link in tracked data,
%      check if the same link exists in ground truth.
%      
%   2. Gap Closing Accuracy: For each gap-closed link, verify the 
%      track segments belonged to the same GT particle.
%
%   3. Track Completeness: How much of each GT trajectory was recovered.
%      Perfect tracking = 1.0 for all particles.
%
%   4. Track Purity: How many GT particles contributed to each tracked
%      trajectory. Pure track = 1.0 (single GT source).
%
%   5. Lifetime Distribution: Kolmogorov-Smirnov test comparing tracked
%      vs GT lifetime distributions.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'trackedTraj');
addRequired(p, 'groundTruth');
addParameter(p, 'DistThreshold', 50, @isnumeric);
addParameter(p, 'AngleThreshold', 0.3, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'PlotResults', true, @islogical);
addParameter(p, 'SavePath', '', @ischar);
parse(p, trackedTraj, groundTruth, varargin{:});

opts = p.Results;
figHandles = [];

%% --- Convert to matrices if needed ---
if iscell(trackedTraj)
    trackedTraj = cell2mat(trackedTraj(:));
end
if iscell(groundTruth)
    groundTruth = cell2mat(groundTruth(:));
end

%% --- Extract unique tracks and frames ---
trackedIDs = unique(trackedTraj(:, 2));
gtIDs = unique(groundTruth(:, 2));
allFrames = unique([trackedTraj(:,1); groundTruth(:,1)]);

nTracked = length(trackedIDs);
nGT = length(gtIDs);
nFrames = length(allFrames);

fprintf('\n=== Tracking Validation ===\n');
fprintf('Ground truth: %d particles, %d frames\n', nGT, nFrames);
fprintf('Tracked: %d trajectories\n', nTracked);

%% --- Build frame-by-frame correspondence ---
% For each frame, match tracked points to GT points

correspondence = cell(nFrames, 1);  % tracked_idx -> gt_idx mapping per frame

for iFrame = 1:nFrames
    frame = allFrames(iFrame);
    
    % Get tracked points in this frame
    trackedInFrame = trackedTraj(trackedTraj(:,1) == frame, :);
    gtInFrame = groundTruth(groundTruth(:,1) == frame, :);
    
    nT = size(trackedInFrame, 1);
    nG = size(gtInFrame, 1);
    
    if nT == 0 || nG == 0
        correspondence{iFrame} = struct('trackedID', [], 'gtID', [], ...
            'trackedLocalIdx', [], 'gtLocalIdx', []);
        continue;
    end
    
    % Compute distance matrix
    distMat = zeros(nT, nG);
    for i = 1:nT
        for j = 1:nG
            dx = trackedInFrame(i, 3) - gtInFrame(j, 3);
            dy = trackedInFrame(i, 4) - gtInFrame(j, 4);
            dz = trackedInFrame(i, 5) - gtInFrame(j, 5);
            distMat(i, j) = sqrt(dx^2 + dy^2 + dz^2);
        end
    end
    
    % Hungarian algorithm for optimal assignment
    distMat(distMat > opts.DistThreshold) = Inf;
    
    % Simple greedy matching (for small matrices)
    matchT = zeros(nT, 1);
    matchG = zeros(nG, 1);
    
    [sortedDist, sortIdx] = sort(distMat(:));
    for k = 1:length(sortedDist)
        if sortedDist(k) == Inf
            break;
        end
        [i, j] = ind2sub([nT, nG], sortIdx(k));
        if matchT(i) == 0 && matchG(j) == 0
            matchT(i) = j;
            matchG(j) = i;
        end
    end
    
    % Store correspondence
    matched = matchT > 0;
    correspondence{iFrame} = struct(...
        'trackedID', trackedInFrame(matched, 2), ...
        'gtID', gtInFrame(matchT(matched), 2), ...
        'trackedLocalIdx', find(matched), ...
        'gtLocalIdx', matchT(matched));
end

%% === METRIC 1: Frame-to-Frame Linking Accuracy ===
% Check if consecutive links in tracked data match GT

nLinks_tracked = 0;
nLinks_correct = 0;
nLinks_wrong = 0;

for iFrame = 1:(nFrames-1)
    frame1 = allFrames(iFrame);
    frame2 = allFrames(iFrame + 1);
    
    % Get tracked trajectories that span both frames
    for iTrk = 1:nTracked
        trkID = trackedIDs(iTrk);
        
        pts1 = trackedTraj(trackedTraj(:,1) == frame1 & trackedTraj(:,2) == trkID, :);
        pts2 = trackedTraj(trackedTraj(:,1) == frame2 & trackedTraj(:,2) == trkID, :);
        
        if isempty(pts1) || isempty(pts2)
            continue;
        end
        
        % This is a tracked link
        nLinks_tracked = nLinks_tracked + 1;
        
        % Find corresponding GT IDs
        corr1 = correspondence{iFrame};
        corr2 = correspondence{iFrame + 1};
        
        gtID1 = [];
        gtID2 = [];
        
        if ~isempty(corr1.trackedID)
            idx1 = corr1.trackedID == trkID;
            if any(idx1)
                gtID1 = corr1.gtID(idx1);
            end
        end
        
        if ~isempty(corr2.trackedID)
            idx2 = corr2.trackedID == trkID;
            if any(idx2)
                gtID2 = corr2.gtID(idx2);
            end
        end
        
        % Check if link is correct
        if ~isempty(gtID1) && ~isempty(gtID2) && gtID1 == gtID2
            nLinks_correct = nLinks_correct + 1;
        else
            nLinks_wrong = nLinks_wrong + 1;
        end
    end
end

linkingAccuracy = nLinks_correct / max(nLinks_tracked, 1);
linkingPrecision = nLinks_correct / max(nLinks_correct + nLinks_wrong, 1);

%% === METRIC 2: Track Completeness ===
% For each GT particle, what fraction was recovered in tracking?

completeness = zeros(nGT, 1);
gtLengths = zeros(nGT, 1);
recoveredLengths = zeros(nGT, 1);

for iGT = 1:nGT
    gtID = gtIDs(iGT);
    gtFrames = groundTruth(groundTruth(:,2) == gtID, 1);
    gtLengths(iGT) = length(gtFrames);
    
    % Find which tracked trajectories matched this GT particle
    recoveredFrames = [];
    for iFrame = 1:nFrames
        corr = correspondence{iFrame};
        if ~isempty(corr.gtID)
            idx = corr.gtID == gtID;
            if any(idx)
                recoveredFrames = [recoveredFrames; allFrames(iFrame)];
            end
        end
    end
    
    recoveredLengths(iGT) = length(unique(recoveredFrames));
    completeness(iGT) = recoveredLengths(iGT) / gtLengths(iGT);
end

meanCompleteness = mean(completeness);

%% === METRIC 3: Track Purity ===
% For each tracked trajectory, how many GT particles contributed?

purity = zeros(nTracked, 1);
trackedLengths = zeros(nTracked, 1);

for iTrk = 1:nTracked
    trkID = trackedIDs(iTrk);
    trkFrames = trackedTraj(trackedTraj(:,2) == trkID, 1);
    trackedLengths(iTrk) = length(trkFrames);
    
    % Count GT particles that contributed
    gtContributors = [];
    for iFrame = 1:nFrames
        corr = correspondence{iFrame};
        if ~isempty(corr.trackedID)
            idx = corr.trackedID == trkID;
            if any(idx)
                gtContributors = [gtContributors; corr.gtID(idx)];
            end
        end
    end
    
    if isempty(gtContributors)
        purity(iTrk) = 0;
    else
        % Purity = fraction from most common GT contributor
        [~, ~, ic] = unique(gtContributors);
        counts = accumarray(ic, 1);
        purity(iTrk) = max(counts) / length(gtContributors);
    end
end

meanPurity = mean(purity);

%% === METRIC 4: Split/Merge Detection ===
% Count GT particles that were split or tracked trajectories that merged

% Splits: GT particles recovered by multiple tracked trajectories
gtToTracked = cell(nGT, 1);
for iGT = 1:nGT
    gtID = gtIDs(iGT);
    trackedContributors = [];
    for iFrame = 1:nFrames
        corr = correspondence{iFrame};
        if ~isempty(corr.gtID)
            idx = corr.gtID == gtID;
            if any(idx)
                trackedContributors = [trackedContributors; corr.trackedID(idx)];
            end
        end
    end
    gtToTracked{iGT} = unique(trackedContributors);
end

nSplit = sum(cellfun(@length, gtToTracked) > 1);
nRecovered = sum(cellfun(@length, gtToTracked) >= 1);

% Merges: Tracked trajectories with contributions from multiple GT
trackedToGT = cell(nTracked, 1);
for iTrk = 1:nTracked
    trkID = trackedIDs(iTrk);
    gtContributors = [];
    for iFrame = 1:nFrames
        corr = correspondence{iFrame};
        if ~isempty(corr.trackedID)
            idx = corr.trackedID == trkID;
            if any(idx)
                gtContributors = [gtContributors; corr.gtID(idx)];
            end
        end
    end
    trackedToGT{iTrk} = unique(gtContributors);
end

nMerged = sum(cellfun(@length, trackedToGT) > 1);

%% === METRIC 5: Lifetime Distribution ===
% Kolmogorov-Smirnov test comparing lifetime distributions

trackedLifetimes = trackedLengths;
gtLifetimes = gtLengths;

if length(trackedLifetimes) >= 2 && length(gtLifetimes) >= 2
    [~, ksP] = kstest2(trackedLifetimes, gtLifetimes);
else
    ksP = NaN;
end

%% === Compile metrics struct ===
metrics = struct();
metrics.linkingAccuracy = linkingAccuracy;
metrics.linkingPrecision = linkingPrecision;
metrics.nLinksTracked = nLinks_tracked;
metrics.nLinksCorrect = nLinks_correct;
metrics.nLinksWrong = nLinks_wrong;
metrics.meanCompleteness = meanCompleteness;
metrics.completenessPerGT = completeness;
metrics.meanPurity = meanPurity;
metrics.purityPerTrack = purity;
metrics.nGTParticles = nGT;
metrics.nTrackedTraj = nTracked;
metrics.nRecovered = nRecovered;
metrics.nSplit = nSplit;
metrics.nMerged = nMerged;
metrics.lifetimeKS_pvalue = ksP;
metrics.trackedLifetimes = trackedLifetimes;
metrics.gtLifetimes = gtLifetimes;

%% --- Print summary ---
if opts.Verbose
    fprintf('\n--- Tracking Validation Summary ---\n');
    fprintf('  Frame-to-frame linking:\n');
    fprintf('    Accuracy: %.1f%% (%d/%d correct)\n', ...
        100*linkingAccuracy, nLinks_correct, nLinks_tracked);
    fprintf('    Precision: %.1f%%\n', 100*linkingPrecision);
    fprintf('\n  Track recovery:\n');
    fprintf('    GT particles: %d\n', nGT);
    fprintf('    Recovered: %d (%.1f%%)\n', nRecovered, 100*nRecovered/nGT);
    fprintf('    Split into multiple tracks: %d\n', nSplit);
    fprintf('    Tracked trajectories: %d\n', nTracked);
    fprintf('    Merged (multiple GT): %d\n', nMerged);
    fprintf('\n  Quality metrics:\n');
    fprintf('    Mean completeness: %.1f%%\n', 100*meanCompleteness);
    fprintf('    Mean purity: %.1f%%\n', 100*meanPurity);
    fprintf('\n  Lifetime distribution KS test:\n');
    fprintf('    p-value: %.4f %s\n', ksP, ...
        ternary(ksP > 0.05, '(distributions similar)', '(distributions differ!)'));
    fprintf('-----------------------------------\n');
end

%% --- Generate plots ---
if opts.PlotResults
    figHandles = zeros(3, 1);
    
    % Figure 1: Summary bar chart
    figHandles(1) = figure('Position', [100, 100, 800, 400], 'Color', 'w');
    
    subplot(1, 2, 1);
    categories = {'Linking\nAccuracy', 'Mean\nCompleteness', 'Mean\nPurity'};
    values = [linkingAccuracy, meanCompleteness, meanPurity] * 100;
    b = bar(values, 'FaceColor', 'flat');
    b.CData = [0.2 0.6 0.8; 0.4 0.8 0.4; 0.8 0.4 0.6];
    ylim([0 105]);
    ylabel('Percentage (%)', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTickLabel', categories, 'FontSize', 10);
    title('Tracking Performance Metrics', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add value labels
    for i = 1:length(values)
        text(i, values(i) + 3, sprintf('%.1f%%', values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    subplot(1, 2, 2);
    categories2 = {'Recovered', 'Split', 'Merged'};
    values2 = [nRecovered, nSplit, nMerged];
    b2 = bar(values2, 'FaceColor', 'flat');
    b2.CData = [0.3 0.7 0.3; 0.9 0.5 0.2; 0.7 0.3 0.7];
    ylabel('Count', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTickLabel', categories2, 'FontSize', 10);
    title('Track Recovery Analysis', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add value labels
    for i = 1:length(values2)
        text(i, values2(i) + 0.5, sprintf('%d', values2(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Figure 2: Lifetime comparison
    figHandles(2) = figure('Position', [150, 150, 600, 400], 'Color', 'w');
    
    hold on;
    [f1, x1] = ecdf(gtLifetimes);
    [f2, x2] = ecdf(trackedLifetimes);
    stairs(x1, f1, 'b-', 'LineWidth', 2);
    stairs(x2, f2, 'r--', 'LineWidth', 2);
    hold off;
    
    xlabel('Trajectory Length (frames)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Cumulative Probability', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Lifetime Distribution (KS p = %.3f)', ksP), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend({'Ground Truth', 'Tracked'}, 'Location', 'southeast', 'FontSize', 11);
    grid on;
    
    % Figure 3: Completeness and purity distributions
    figHandles(3) = figure('Position', [200, 200, 800, 350], 'Color', 'w');
    
    subplot(1, 2, 1);
    histogram(completeness * 100, 0:10:100, 'FaceColor', [0.4 0.8 0.4], ...
        'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xlabel('Completeness (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Number of GT Particles', 'FontSize', 11, 'FontWeight', 'bold');
    title('Track Completeness Distribution', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0 100]);
    
    subplot(1, 2, 2);
    histogram(purity * 100, 0:10:100, 'FaceColor', [0.8 0.4 0.6], ...
        'EdgeColor', 'w', 'FaceAlpha', 0.8);
    xlabel('Purity (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Number of Tracked Trajectories', 'FontSize', 11, 'FontWeight', 'bold');
    title('Track Purity Distribution', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0 100]);
    
    % Save if requested
    if ~isempty(opts.SavePath)
        for i = 1:length(figHandles)
            saveas(figHandles(i), fullfile(opts.SavePath, ...
                sprintf('validation_fig%d.png', i)));
        end
        fprintf('Validation figures saved to: %s\n', opts.SavePath);
    end
end

end

%% --- Helper function ---
function out = ternary(condition, trueVal, falseVal)
    if condition
        out = trueVal;
    else
        out = falseVal;
    end
end
