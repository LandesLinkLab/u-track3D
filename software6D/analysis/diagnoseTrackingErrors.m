function [errorReport, figHandle] = diagnoseTrackingErrors(trajTracked, groundTruth, varargin)
%DIAGNOSETRACKINGERRORS Compare tracked trajectories to ground truth and identify errors.
%
%   Analyzes tracking results to find:
%   - Track fragmentation (one true particle split into multiple tracks)
%   - Track merging (multiple true particles merged into one track)
%   - ID swaps (track switches from one particle to another mid-trajectory)
%   - Spatial jumps (sudden large displacements suggesting linking errors)
%   - Orientation jumps (sudden large angular changes)
%
%   USAGE:
%       [errorReport, figHandle] = diagnoseTrackingErrors(trajTracked, groundTruth)
%       [errorReport, figHandle] = diagnoseTrackingErrors(..., 'Name', Value)
%
%   INPUTS:
%       trajTracked : [N x 9] tracked trajectory matrix
%                     Columns: [frame, trackID, x, y, z, theta, phi, omega, signal]
%       groundTruth : struct from generate6DTrajectories containing:
%                     .trajectories - cell array of true trajectories
%                     .NumParticles - number of simulated particles
%
%   OPTIONAL PARAMETERS:
%       'SpatialJumpThreshold' : Distance (nm) to flag as spatial jump (default: 200)
%       'AngularJumpThreshold' : Angle (degrees) to flag as orientation jump (default: 30)
%       'MatchRadius'          : Max distance (nm) to match detection to GT (default: 100)
%       'Verbose'              : Print detailed output (default: true)
%       'PlotResults'          : Generate diagnostic figures (default: true)
%
%   OUTPUTS:
%       errorReport : struct containing:
%           .trackAssignments - which GT particle each track point belongs to
%           .fragmentedTracks - tracks split from same particle
%           .mergedTracks     - tracks containing multiple particles
%           .idSwaps          - frames where track switches particles
%           .spatialJumps     - frames with large spatial jumps
%           .orientJumps      - frames with large orientation jumps
%           .summary          - text summary of errors
%       figHandle   : handle to diagnostic figure (if PlotResults=true)
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% Parse inputs
p = inputParser;
addRequired(p, 'trajTracked', @isnumeric);
addRequired(p, 'groundTruth', @isstruct);
addParameter(p, 'SpatialJumpThreshold', 200, @isnumeric);  % nm
addParameter(p, 'AngularJumpThreshold', 30, @isnumeric);   % degrees
addParameter(p, 'MatchRadius', 100, @isnumeric);           % nm
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'PlotResults', true, @islogical);
parse(p, trajTracked, groundTruth, varargin{:});

spatialThresh = p.Results.SpatialJumpThreshold;
angularThresh = p.Results.AngularJumpThreshold * pi/180;  % Convert to radians
matchRadius = p.Results.MatchRadius;
verbose = p.Results.Verbose;
doPlot = p.Results.PlotResults;

%% Initialize output
errorReport = struct();
errorReport.trackAssignments = [];
errorReport.fragmentedTracks = {};
errorReport.mergedTracks = {};
errorReport.idSwaps = [];
errorReport.spatialJumps = [];
errorReport.orientJumps = [];
errorReport.summary = '';
figHandle = [];

%% Build ground truth lookup table
% Combine all GT trajectories into one matrix for easy lookup
if verbose
    fprintf('\n=== TRACKING ERROR DIAGNOSIS ===\n\n');
    fprintf('Building ground truth lookup table...\n');
end

% Get number of particles from groundTruth structure
% Handle both direct field and params substructure
if isfield(groundTruth, 'NumParticles')
    nParticles = groundTruth.NumParticles;
elseif isfield(groundTruth, 'params') && isfield(groundTruth.params, 'NumParticles')
    nParticles = groundTruth.params.NumParticles;
elseif isfield(groundTruth, 'trajectories')
    nParticles = length(groundTruth.trajectories);
else
    error('Cannot determine number of particles from groundTruth structure');
end

gtAll = [];
for i = 1:nParticles
    gtTraj = groundTruth.trajectories{i};
    % GT columns: [frame, signal, x, y, z, theta_deg, phi_deg, gamma, trackID, isInterp]
    % Add particle ID as extra column
    gtWithID = [gtTraj, repmat(i, size(gtTraj, 1), 1)];
    gtAll = [gtAll; gtWithID];
end

% Create frame-indexed lookup
nFrames = max(gtAll(:, 1));
gtByFrame = cell(nFrames, 1);
for f = 1:nFrames
    gtByFrame{f} = gtAll(gtAll(:, 1) == f, :);
end

if verbose
    fprintf('  Ground truth: %d particles, %d frames, %d total points\n', ...
        nParticles, nFrames, size(gtAll, 1));
end

%% Match each tracked point to ground truth
if verbose
    fprintf('Matching tracked points to ground truth...\n');
end

nTrackedPts = size(trajTracked, 1);
% Columns: [frame, trackID, x, y, z, theta, phi, omega, signal]
trackAssignments = zeros(nTrackedPts, 1);  % Which GT particle each point matches
matchDistances = zeros(nTrackedPts, 1);    % Distance to matched GT point

for i = 1:nTrackedPts
    frame = trajTracked(i, 1);
    x = trajTracked(i, 3);
    y = trajTracked(i, 4);
    z = trajTracked(i, 5);
    
    if frame < 1 || frame > nFrames || isempty(gtByFrame{frame})
        trackAssignments(i) = 0;  % No match
        matchDistances(i) = Inf;
        continue;
    end
    
    % Find closest GT point in this frame
    gtFrame = gtByFrame{frame};
    gtX = gtFrame(:, 3);
    gtY = gtFrame(:, 4);
    gtZ = gtFrame(:, 5);
    
    dists = sqrt((gtX - x).^2 + (gtY - y).^2 + (gtZ - z).^2);
    [minDist, minIdx] = min(dists);
    
    if minDist <= matchRadius
        trackAssignments(i) = gtFrame(minIdx, end);  % GT particle ID
        matchDistances(i) = minDist;
    else
        trackAssignments(i) = 0;  % No match within radius
        matchDistances(i) = minDist;
    end
end

errorReport.trackAssignments = trackAssignments;

nMatched = sum(trackAssignments > 0);
nUnmatched = sum(trackAssignments == 0);
if verbose
    fprintf('  Matched: %d points (%.1f%%)\n', nMatched, 100*nMatched/nTrackedPts);
    fprintf('  Unmatched: %d points (dist > %d nm)\n', nUnmatched, matchRadius);
end

%% Analyze track composition (fragmentation and merging)
if verbose
    fprintf('Analyzing track composition...\n');
end

trackIDs = unique(trajTracked(:, 2));
nTracks = length(trackIDs);

% For each track, count how many GT particles it contains
trackComposition = cell(nTracks, 1);
for t = 1:nTracks
    tid = trackIDs(t);
    mask = trajTracked(:, 2) == tid;
    gtIDs = trackAssignments(mask);
    gtIDs = gtIDs(gtIDs > 0);  % Remove unmatched
    
    if isempty(gtIDs)
        trackComposition{t} = struct('trackID', tid, 'gtParticles', [], 'counts', []);
    else
        [uniqueGT, ~, ic] = unique(gtIDs);
        counts = accumarray(ic, 1);
        trackComposition{t} = struct('trackID', tid, 'gtParticles', uniqueGT', 'counts', counts');
    end
end

% Find merged tracks (contain multiple GT particles)
mergedTracks = {};
for t = 1:nTracks
    if length(trackComposition{t}.gtParticles) > 1
        mergedTracks{end+1} = trackComposition{t};
    end
end
errorReport.mergedTracks = mergedTracks;

% Find fragmented particles (GT particle split across multiple tracks)
particleToTracks = cell(nParticles, 1);
for p = 1:nParticles
    tracksWithParticle = [];
    for t = 1:nTracks
        if ismember(p, trackComposition{t}.gtParticles)
            tracksWithParticle = [tracksWithParticle, trackComposition{t}.trackID];
        end
    end
    particleToTracks{p} = tracksWithParticle;
end

fragmentedTracks = {};
for p = 1:nParticles
    if length(particleToTracks{p}) > 1
        fragmentedTracks{end+1} = struct('gtParticle', p, 'tracks', particleToTracks{p});
    end
end
errorReport.fragmentedTracks = fragmentedTracks;

if verbose
    fprintf('  Tracks containing multiple GT particles (merging): %d\n', length(mergedTracks));
    fprintf('  GT particles split across multiple tracks (fragmentation): %d\n', length(fragmentedTracks));
end

%% Find ID swaps within tracks
if verbose
    fprintf('Finding ID swaps within tracks...\n');
end

idSwaps = [];
for t = 1:nTracks
    tid = trackIDs(t);
    mask = trajTracked(:, 2) == tid;
    trackData = trajTracked(mask, :);
    trackGT = trackAssignments(mask);
    
    % Sort by frame
    [~, sortIdx] = sort(trackData(:, 1));
    trackData = trackData(sortIdx, :);
    trackGT = trackGT(sortIdx);
    
    % Find where GT assignment changes
    for i = 2:length(trackGT)
        if trackGT(i) > 0 && trackGT(i-1) > 0 && trackGT(i) ~= trackGT(i-1)
            swapInfo = struct();
            swapInfo.trackID = tid;
            swapInfo.frame = trackData(i, 1);
            swapInfo.fromGT = trackGT(i-1);
            swapInfo.toGT = trackGT(i);
            swapInfo.position = trackData(i, 3:5);
            idSwaps = [idSwaps; swapInfo];
        end
    end
end
errorReport.idSwaps = idSwaps;

if verbose
    fprintf('  ID swaps detected: %d\n', length(idSwaps));
end

%% Find spatial jumps
if verbose
    fprintf('Finding spatial jumps (threshold: %d nm)...\n', spatialThresh);
end

spatialJumps = [];
for t = 1:nTracks
    tid = trackIDs(t);
    mask = trajTracked(:, 2) == tid;
    trackData = trajTracked(mask, :);
    
    % Sort by frame
    [~, sortIdx] = sort(trackData(:, 1));
    trackData = trackData(sortIdx, :);
    
    for i = 2:size(trackData, 1)
        dx = trackData(i, 3) - trackData(i-1, 3);
        dy = trackData(i, 4) - trackData(i-1, 4);
        dz = trackData(i, 5) - trackData(i-1, 5);
        dist = sqrt(dx^2 + dy^2 + dz^2);
        
        % Scale by frame gap (allow larger jumps for gaps)
        frameGap = trackData(i, 1) - trackData(i-1, 1);
        scaledThresh = spatialThresh * sqrt(frameGap);
        
        if dist > scaledThresh
            jumpInfo = struct();
            jumpInfo.trackID = tid;
            jumpInfo.frame = trackData(i, 1);
            jumpInfo.frameGap = frameGap;
            jumpInfo.distance = dist;
            jumpInfo.threshold = scaledThresh;
            jumpInfo.positionBefore = trackData(i-1, 3:5);
            jumpInfo.positionAfter = trackData(i, 3:5);
            spatialJumps = [spatialJumps; jumpInfo];
        end
    end
end
errorReport.spatialJumps = spatialJumps;

if verbose
    fprintf('  Spatial jumps detected: %d\n', length(spatialJumps));
end

%% Find orientation jumps
if verbose
    fprintf('Finding orientation jumps (threshold: %.0f°)...\n', angularThresh * 180/pi);
end

orientJumps = [];
for t = 1:nTracks
    tid = trackIDs(t);
    mask = trajTracked(:, 2) == tid;
    trackData = trajTracked(mask, :);
    
    % Sort by frame
    [~, sortIdx] = sort(trackData(:, 1));
    trackData = trackData(sortIdx, :);
    
    for i = 2:size(trackData, 1)
        theta1 = trackData(i-1, 6);
        phi1 = trackData(i-1, 7);
        theta2 = trackData(i, 6);
        phi2 = trackData(i, 7);
        
        % Compute angular distance (with head-tail symmetry)
        mu1 = [sin(theta1)*cos(phi1), sin(theta1)*sin(phi1), cos(theta1)];
        mu2 = [sin(theta2)*cos(phi2), sin(theta2)*sin(phi2), cos(theta2)];
        dotProd = abs(dot(mu1, mu2));
        dAngle = acos(min(dotProd, 1));
        
        % Scale by frame gap
        frameGap = trackData(i, 1) - trackData(i-1, 1);
        scaledThresh = angularThresh * sqrt(frameGap);
        
        if dAngle > scaledThresh
            jumpInfo = struct();
            jumpInfo.trackID = tid;
            jumpInfo.frame = trackData(i, 1);
            jumpInfo.frameGap = frameGap;
            jumpInfo.angularDistance = dAngle * 180/pi;  % degrees
            jumpInfo.threshold = scaledThresh * 180/pi;
            jumpInfo.orientBefore = [theta1, phi1] * 180/pi;
            jumpInfo.orientAfter = [theta2, phi2] * 180/pi;
            orientJumps = [orientJumps; jumpInfo];
        end
    end
end
errorReport.orientJumps = orientJumps;

if verbose
    fprintf('  Orientation jumps detected: %d\n', length(orientJumps));
end

%% Generate summary
summary = sprintf('=== TRACKING ERROR SUMMARY ===\n');
summary = [summary, sprintf('Ground truth: %d particles\n', nParticles)];
summary = [summary, sprintf('Tracked: %d tracks, %d points\n', nTracks, nTrackedPts)];
summary = [summary, sprintf('Matched to GT: %d points (%.1f%%)\n\n', nMatched, 100*nMatched/nTrackedPts)];

summary = [summary, sprintf('FRAGMENTATION (1 particle → multiple tracks):\n')];
if isempty(fragmentedTracks)
    summary = [summary, sprintf('  None\n')];
else
    for i = 1:length(fragmentedTracks)
        ft = fragmentedTracks{i};
        summary = [summary, sprintf('  Particle %d → Tracks [%s]\n', ...
            ft.gtParticle, num2str(ft.tracks))];
    end
end

summary = [summary, sprintf('\nMERGING (multiple particles → 1 track):\n')];
if isempty(mergedTracks)
    summary = [summary, sprintf('  None\n')];
else
    for i = 1:length(mergedTracks)
        mt = mergedTracks{i};
        summary = [summary, sprintf('  Track %d contains particles [%s]\n', ...
            mt.trackID, num2str(mt.gtParticles))];
    end
end

summary = [summary, sprintf('\nID SWAPS (track switches particles mid-trajectory):\n')];
if isempty(idSwaps)
    summary = [summary, sprintf('  None\n')];
else
    for i = 1:length(idSwaps)
        sw = idSwaps(i);
        summary = [summary, sprintf('  Track %d @ frame %d: particle %d → %d\n', ...
            sw.trackID, sw.frame, sw.fromGT, sw.toGT)];
    end
end

summary = [summary, sprintf('\nSPATIAL JUMPS (> %d nm scaled):\n', spatialThresh)];
if isempty(spatialJumps)
    summary = [summary, sprintf('  None\n')];
else
    for i = 1:min(10, length(spatialJumps))  % Show first 10
        sj = spatialJumps(i);
        summary = [summary, sprintf('  Track %d @ frame %d: %.0f nm (gap=%d frames)\n', ...
            sj.trackID, sj.frame, sj.distance, sj.frameGap)];
    end
    if length(spatialJumps) > 10
        summary = [summary, sprintf('  ... and %d more\n', length(spatialJumps) - 10)];
    end
end

summary = [summary, sprintf('\nORIENTATION JUMPS (> %.0f° scaled):\n', angularThresh * 180/pi)];
if isempty(orientJumps)
    summary = [summary, sprintf('  None\n')];
else
    for i = 1:min(10, length(orientJumps))  % Show first 10
        oj = orientJumps(i);
        summary = [summary, sprintf('  Track %d @ frame %d: %.1f° (gap=%d frames)\n', ...
            oj.trackID, oj.frame, oj.angularDistance, oj.frameGap)];
    end
    if length(orientJumps) > 10
        summary = [summary, sprintf('  ... and %d more\n', length(orientJumps) - 10)];
    end
end

errorReport.summary = summary;

if verbose
    fprintf('\n%s\n', summary);
end

%% Plot results
if doPlot
    figHandle = figure('Position', [50, 50, 1600, 900], 'Color', 'w', ...
        'Name', 'Tracking Error Diagnosis');
    
    % === STANDARD FORMATTING CONSTANTS ===
    FMT_FONTSIZE_TITLE = 14;
    FMT_FONTSIZE_LABEL = 13;
    FMT_FONTSIZE_TICK = 12;
    FMT_LINEWIDTH_AXES = 2.0;
    FMT_FONTWEIGHT = 'bold';
    FMT_EXPORT_DPI = 600;
    
    % Helper function for axes formatting (includes grid)
    applyFmt = @(ax) set(ax, 'FontSize', FMT_FONTSIZE_TICK, ...
        'FontWeight', FMT_FONTWEIGHT, 'LineWidth', FMT_LINEWIDTH_AXES, ...
        'TickLength', [0.02 0.02], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
    
    % Color map for GT particles
    gtColors = lines(nParticles);
    
    % --- Panel 1: XY trajectories colored by GT assignment ---
    ax1 = subplot(2, 3, 1);
    hold on;
    for p = 1:nParticles
        mask = trackAssignments == p;
        if any(mask)
            scatter(trajTracked(mask, 3), trajTracked(mask, 4), 15, ...
                gtColors(p, :), 'filled', 'MarkerFaceAlpha', 0.6);
        end
    end
    % Unmatched in gray
    mask = trackAssignments == 0;
    if any(mask)
        scatter(trajTracked(mask, 3), trajTracked(mask, 4), 10, ...
            [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.3);
    end
    % Mark ID swaps with red X
    for i = 1:length(idSwaps)
        plot(idSwaps(i).position(1), idSwaps(i).position(2), 'rx', ...
            'MarkerSize', 12, 'LineWidth', 2);
    end
    hold off;
    xlabel('X (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    ylabel('Y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('XY Colored by GT Particle', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    axis equal; grid on;
    applyFmt(gca);
    
    % Legend
    legendStr = {};
    for p = 1:nParticles
        legendStr{end+1} = sprintf('Particle %d', p);
    end
    legendStr{end+1} = 'Unmatched';
    if ~isempty(idSwaps)
        legendStr{end+1} = 'ID Swap';
    end
    legend(legendStr, 'Location', 'best', 'FontSize', FMT_FONTSIZE_TICK-1);
    
    % --- Panel 2: XY trajectories colored by track ID ---
    % Use GT colors if track maps to single particle, otherwise use track colors
    ax2 = subplot(2, 3, 2);
    hold on;
    trackColors = lines(nTracks);  % Fallback colors
    
    % Determine dominant GT particle for each track
    trackToGT = zeros(nTracks, 1);
    for t = 1:nTracks
        tid = trackIDs(t);
        mask = trajTracked(:, 2) == tid;
        gtAssign = trackAssignments(mask);
        gtAssign = gtAssign(gtAssign > 0);  % Remove unmatched
        if ~isempty(gtAssign)
            % Get most common GT assignment
            trackToGT(t) = mode(gtAssign);
        end
    end
    
    for t = 1:nTracks
        tid = trackIDs(t);
        mask = trajTracked(:, 2) == tid;
        
        % Use GT color if track maps to a GT particle
        if trackToGT(t) > 0
            plotColor = gtColors(trackToGT(t), :);
        else
            plotColor = [0.5 0.5 0.5];  % Gray for unmatched
        end
        
        scatter(trajTracked(mask, 3), trajTracked(mask, 4), 15, ...
            plotColor, 'filled', 'MarkerFaceAlpha', 0.6);
    end
    % Mark spatial jumps with black circles
    for i = 1:length(spatialJumps)
        plot(spatialJumps(i).positionAfter(1), spatialJumps(i).positionAfter(2), 'ko', ...
            'MarkerSize', 15, 'LineWidth', 2);
    end
    hold off;
    xlabel('X (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    ylabel('Y (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title(sprintf('XY Colored by Track ID (%d tracks)', nTracks), 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    axis equal; grid on;
    applyFmt(gca);
    
    % --- Panel 3: Time series of GT assignment per track ---
    ax3 = subplot(2, 3, 3);
    hold on;
    for t = 1:nTracks
        tid = trackIDs(t);
        mask = trajTracked(:, 2) == tid;
        frames = trajTracked(mask, 1);
        gtAssign = trackAssignments(mask);
        
        % Plot GT assignment (no offset - align with y-ticks)
        plot(frames, gtAssign, '-', 'Color', trackColors(t, :), 'LineWidth', 2);
    end
    hold off;
    xlabel('Frame', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    ylabel('GT Particle Assignment', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
    title('GT Assignment Over Time', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    ylim([0.5, nParticles + 0.5]);
    yticks(1:nParticles);
    grid on;
    applyFmt(gca);
    
    % --- Panel 4: Spatial displacement histogram ---
    ax4 = subplot(2, 3, 4);
    allDisplacements = [];
    for t = 1:nTracks
        tid = trackIDs(t);
        mask = trajTracked(:, 2) == tid;
        trackData = trajTracked(mask, :);
        [~, sortIdx] = sort(trackData(:, 1));
        trackData = trackData(sortIdx, :);
        
        for i = 2:size(trackData, 1)
            frameGap = trackData(i, 1) - trackData(i-1, 1);
            if frameGap == 1  % Only consecutive frames
                dx = trackData(i, 3) - trackData(i-1, 3);
                dy = trackData(i, 4) - trackData(i-1, 4);
                dz = trackData(i, 5) - trackData(i-1, 5);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                allDisplacements = [allDisplacements; dist];
            end
        end
    end
    if ~isempty(allDisplacements) && any(isfinite(allDisplacements))
        maxDisp = max(allDisplacements(isfinite(allDisplacements)));
        % Set x-axis limit based on data, with minimum range
        xMax = max(maxDisp * 1.2, 1);  % At least 1 nm range
        
        histogram(allDisplacements, 50, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w');
        if xMax > 0 && isfinite(xMax)
            xlim([0, xMax]);  % Set limit based on data
        end
        hold on;
        % Only show threshold line if it's within the visible range
        if spatialThresh <= xMax && spatialThresh > 0
            xline(spatialThresh, 'r--', 'LineWidth', 2, 'Label', sprintf('Thresh=%.0f', spatialThresh));
        end
        hold off;
        xlabel('Frame-to-frame displacement (nm)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Spatial Displacements (mean=%.0f nm, n=%d)', mean(allDisplacements), length(allDisplacements)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyFmt(gca);
    else
        text(0.5, 0.5, 'No consecutive frame data', 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
            'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Spatial Displacements: No data', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    grid on;
    
    % --- Panel 5: Angular displacement histogram ---
    ax5 = subplot(2, 3, 5);
    allAngular = [];
    for t = 1:nTracks
        tid = trackIDs(t);
        mask = trajTracked(:, 2) == tid;
        trackData = trajTracked(mask, :);
        [~, sortIdx] = sort(trackData(:, 1));
        trackData = trackData(sortIdx, :);
        
        for i = 2:size(trackData, 1)
            frameGap = trackData(i, 1) - trackData(i-1, 1);
            if frameGap == 1  % Only consecutive frames
                theta1 = trackData(i-1, 6);
                phi1 = trackData(i-1, 7);
                theta2 = trackData(i, 6);
                phi2 = trackData(i, 7);
                
                mu1 = [sin(theta1)*cos(phi1), sin(theta1)*sin(phi1), cos(theta1)];
                mu2 = [sin(theta2)*cos(phi2), sin(theta2)*sin(phi2), cos(theta2)];
                dotProd = abs(dot(mu1, mu2));
                dAngle = acos(min(dotProd, 1)) * 180/pi;
                allAngular = [allAngular; dAngle];
            end
        end
    end
    angularThreshDeg = angularThresh * 180/pi;
    if ~isempty(allAngular) && any(isfinite(allAngular))
        maxAngular = max(allAngular(isfinite(allAngular)));
        % Set x-axis limit based on data, with minimum range
        xMax = max(maxAngular * 1.2, 1);  % At least 1 degree range
        
        histogram(allAngular, 50, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'w');
        if xMax > 0 && isfinite(xMax)
            xlim([0, xMax]);  % Set limit based on data
        end
        hold on;
        % Only show threshold line if it's within the visible range
        if angularThreshDeg <= xMax && angularThreshDeg > 0
            xline(angularThreshDeg, 'r--', 'LineWidth', 2, 'Label', sprintf('Thresh=%.0f°', angularThreshDeg));
        end
        hold off;
        xlabel('Frame-to-frame angular change (°)', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        ylabel('Count', 'FontSize', FMT_FONTSIZE_LABEL, 'FontWeight', FMT_FONTWEIGHT);
        title(sprintf('Angular Displacements (mean=%.1f°, n=%d)', mean(allAngular), length(allAngular)), ...
            'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
        applyFmt(gca);
    else
        text(0.5, 0.5, 'No consecutive frame data', 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
            'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT);
        title('Angular Displacements: No data', 'FontSize', FMT_FONTSIZE_TITLE, 'FontWeight', FMT_FONTWEIGHT);
    end
    grid on;
    
    % --- Panel 6: Error details with visual representation ---
    ax6 = subplot(2, 3, 6);
    axis off;
    
    % Build detailed error text
    yPos = 0.95;
    lineHeight = 0.065;
    
    text(0.02, yPos, 'ERROR SUMMARY', 'FontWeight', 'bold', 'FontSize', FMT_FONTSIZE_TITLE, ...
        'VerticalAlignment', 'top');
    yPos = yPos - lineHeight * 1.2;
    
    % Basic stats
    text(0.02, yPos, sprintf('GT particles: %d → Tracked: %d tracks', nParticles, nTracks), ...
        'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT, 'VerticalAlignment', 'top');
    yPos = yPos - lineHeight;
    
    text(0.02, yPos, sprintf('Matched to GT: %.1f%%', 100*nMatched/nTrackedPts), ...
        'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT, 'VerticalAlignment', 'top');
    yPos = yPos - lineHeight * 1.5;
    
    % Fragmentation details
    if ~isempty(fragmentedTracks)
        text(0.02, yPos, sprintf('FRAGMENTATION (%d particles split):', length(fragmentedTracks)), ...
            'FontWeight', 'bold', 'FontSize', 9, 'Color', [0.8 0.4 0], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
        for i = 1:min(5, length(fragmentedTracks))
            ft = fragmentedTracks{i};
            trackStr = sprintf('%d', ft.tracks(1));
            for j = 2:length(ft.tracks)
                trackStr = [trackStr, sprintf(',%d', ft.tracks(j))];
            end
            text(0.04, yPos, sprintf('P%d → T[%s]', ft.gtParticle, trackStr), ...
                'FontSize', 8, 'FontName', 'FixedWidth', 'VerticalAlignment', 'top');
            yPos = yPos - lineHeight * 0.8;
        end
        if length(fragmentedTracks) > 5
            text(0.04, yPos, sprintf('... +%d more', length(fragmentedTracks)-5), ...
                'FontSize', 8, 'VerticalAlignment', 'top');
            yPos = yPos - lineHeight * 0.8;
        end
    else
        text(0.02, yPos, 'Fragmentation: None', ...
            'FontSize', 9, 'Color', [0 0.6 0], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
    end
    yPos = yPos - lineHeight * 0.5;
    
    % Merging details
    if ~isempty(mergedTracks)
        text(0.02, yPos, sprintf('MERGING (%d tracks merged):', length(mergedTracks)), ...
            'FontWeight', 'bold', 'FontSize', 9, 'Color', [0.8 0 0], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
        for i = 1:min(3, length(mergedTracks))
            mt = mergedTracks{i};
            particleStr = sprintf('%d', mt.gtParticles(1));
            for j = 2:length(mt.gtParticles)
                particleStr = [particleStr, sprintf(',%d', mt.gtParticles(j))];
            end
            text(0.04, yPos, sprintf('T%d ← P[%s]', mt.trackID, particleStr), ...
                'FontSize', 8, 'FontName', 'FixedWidth', 'VerticalAlignment', 'top');
            yPos = yPos - lineHeight * 0.8;
        end
    else
        text(0.02, yPos, 'Merging: None', ...
            'FontSize', 9, 'Color', [0 0.6 0], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
    end
    yPos = yPos - lineHeight * 0.5;
    
    % ID swaps
    if ~isempty(idSwaps)
        text(0.02, yPos, sprintf('ID SWAPS (%d):', length(idSwaps)), ...
            'FontWeight', 'bold', 'FontSize', 9, 'Color', [0.6 0 0.6], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
        for i = 1:min(3, length(idSwaps))
            sw = idSwaps(i);
            text(0.04, yPos, sprintf('T%d@f%d: P%d→P%d', sw.trackID, sw.frame, sw.fromGT, sw.toGT), ...
                'FontSize', 8, 'FontName', 'FixedWidth', 'VerticalAlignment', 'top');
            yPos = yPos - lineHeight * 0.8;
        end
    else
        text(0.02, yPos, 'ID swaps: None', ...
            'FontSize', 9, 'Color', [0 0.6 0], 'VerticalAlignment', 'top');
        yPos = yPos - lineHeight;
    end
    yPos = yPos - lineHeight * 0.5;
    
    % Jumps summary
    jumpText = sprintf('Jumps: %d spatial, %d orient', length(spatialJumps), length(orientJumps));
    if isempty(spatialJumps) && isempty(orientJumps)
        text(0.02, yPos, [jumpText ''], 'FontSize', 9, 'Color', [0 0.6 0], 'VerticalAlignment', 'top');
    else
        text(0.02, yPos, jumpText, 'FontSize', FMT_FONTSIZE_TICK, 'FontWeight', FMT_FONTWEIGHT, 'Color', [0.8 0.4 0], 'VerticalAlignment', 'top');
    end
    
    sgtitle('Tracking Error Diagnosis', 'FontWeight', 'bold', 'FontSize', 16);
    
    % Store axes handles for potential panel export
    errorReport.figAxes = {ax1, ax2, ax3, ax4, ax5, ax6};
    errorReport.figPanelNames = {'XY_by_GT_Particle', 'XY_by_Track_ID', 'GT_Assignment_Time', ...
                                  'Spatial_Displacement_Hist', 'Angular_Displacement_Hist', 'Error_Summary'};
end

end
