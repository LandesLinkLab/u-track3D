function [figHandle, costData] = visualizeCostMatrix(costMat, varargin)
%VISUALIZECOSTMATRIX Visualize tracking cost matrix as heatmap
%
%   Creates publication-quality visualizations of the LAP cost matrix
%   used in particle tracking. Useful for debugging and understanding
%   tracking decisions.
%
%   Based on Jaqaman et al., Nature Methods 2008 - Figure 1.
%
%   SYNOPSIS:
%       [figHandle, costData] = visualizeCostMatrix(costMat)
%       [figHandle, costData] = visualizeCostMatrix(costMat, 'Name', Value, ...)
%
%   INPUT:
%       costMat : Sparse or full cost matrix from tracking
%                 Can be frame-to-frame linking or gap closing matrix
%
%   NAME-VALUE PARAMETERS:
%       'Type'           : 'linking' or 'gapclosing'. Default: auto-detect
%       'TrackLabels'    : Cell array of track segment labels. Default: numeric
%       'Title'          : Figure title. Default: 'Cost Matrix'
%       'Colormap'       : Colormap name. Default: 'hot'
%       'LogScale'       : Use log scale for costs. Default: true
%       'ShowValues'     : Show cost values in cells. Default: false (auto for small matrices)
%       'NonlinkMarker'  : Value indicating impossible link. Default: -1
%       'MaxDisplay'     : Max matrix size to display. Default: 100
%       'SavePath'       : Path to save figure. Default: '' (no save)
%
%   OUTPUT:
%       figHandle : Handle to the figure
%       costData  : Struct with cost statistics
%           .minCost      : Minimum valid cost
%           .maxCost      : Maximum valid cost
%           .meanCost     : Mean valid cost
%           .medianCost   : Median valid cost
%           .nValidLinks  : Number of valid (non-infinite) links
%           .sparsity     : Fraction of impossible links
%           .costDistribution : Histogram of costs
%
%   COST MATRIX STRUCTURE (Jaqaman 2008):
%
%   Frame-to-frame linking:
%       ┌─────────────────┬─────────────────┐
%       │    λ_ij         │      d          │
%       │  (link costs)   │  (death costs)  │
%       ├─────────────────┼─────────────────┤
%       │      b          │    auxiliary    │
%       │ (birth costs)   │     block       │
%       └─────────────────┴─────────────────┘
%
%   Gap closing, merging, splitting:
%       ┌─────────┬─────────┬─────────┬─────────┐
%       │  g_IJ   │   ×     │    d    │    ×    │
%       │ (gaps)  │         │         │         │
%       ├─────────┼─────────┼─────────┼─────────┤
%       │  m_IJ   │   ×     │    ×    │   d'    │
%       │(merges) │         │         │         │
%       ├─────────┼─────────┼─────────┼─────────┤
%       │    b    │   ×     │  aux    │    ×    │
%       │         │         │         │         │
%       ├─────────┼─────────┼─────────┼─────────┤
%       │    ×    │   b'    │    ×    │  aux    │
%       │         │         │         │         │
%       └─────────┴─────────┴─────────┴─────────┘
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'costMat');
addParameter(p, 'Type', 'auto', @ischar);
addParameter(p, 'TrackLabels', {}, @iscell);
addParameter(p, 'Title', 'Cost Matrix', @ischar);
addParameter(p, 'Colormap', 'hot', @ischar);
addParameter(p, 'LogScale', true, @islogical);
addParameter(p, 'ShowValues', [], @islogical);
addParameter(p, 'NonlinkMarker', -1, @isnumeric);
addParameter(p, 'MaxDisplay', 100, @isnumeric);
addParameter(p, 'SavePath', '', @ischar);
addParameter(p, 'OrientationCosts', [], @isnumeric);  % For 6D visualization
addParameter(p, 'SpatialCosts', [], @isnumeric);      % For 6D visualization
parse(p, costMat, varargin{:});

opts = p.Results;

%% --- Convert sparse to full if needed ---
if issparse(costMat)
    costMat = full(costMat);
end

[nRows, nCols] = size(costMat);

%% --- Identify valid links ---
% Valid links have finite, positive costs (not nonlinkMarker)
validMask = isfinite(costMat) & costMat > 0 & costMat ~= opts.NonlinkMarker;
validCosts = costMat(validMask);

%% --- Compute statistics ---
costData = struct();
costData.matrixSize = [nRows, nCols];
costData.nValidLinks = numel(validCosts);
costData.sparsity = 1 - costData.nValidLinks / numel(costMat);

if ~isempty(validCosts)
    costData.minCost = min(validCosts);
    costData.maxCost = max(validCosts);
    costData.meanCost = mean(validCosts);
    costData.medianCost = median(validCosts);
    costData.stdCost = std(validCosts);
    
    % Histogram of costs
    if opts.LogScale && costData.minCost > 0
        edges = logspace(log10(costData.minCost), log10(costData.maxCost), 50);
        [counts, edges] = histcounts(validCosts, edges);
    else
        [counts, edges] = histcounts(validCosts, 50);
    end
    costData.costDistribution.counts = counts;
    costData.costDistribution.edges = edges;
else
    costData.minCost = NaN;
    costData.maxCost = NaN;
    costData.meanCost = NaN;
    costData.medianCost = NaN;
    costData.stdCost = NaN;
    costData.costDistribution = [];
end

%% --- Prepare display matrix ---
displayMat = costMat;

% Replace invalid entries with NaN for visualization
displayMat(~validMask) = NaN;

% Downsample if too large
if nRows > opts.MaxDisplay || nCols > opts.MaxDisplay
    stepR = ceil(nRows / opts.MaxDisplay);
    stepC = ceil(nCols / opts.MaxDisplay);
    displayMat = displayMat(1:stepR:end, 1:stepC:end);
    warning('Matrix downsampled from [%d x %d] to [%d x %d] for display', ...
        nRows, nCols, size(displayMat, 1), size(displayMat, 2));
end

% Log transform for better visualization
if opts.LogScale && ~isempty(validCosts) && costData.minCost > 0
    displayMat = log10(displayMat);
    colorLabel = 'log_{10}(Cost)';
else
    colorLabel = 'Cost';
end

%% --- Determine whether to show values ---
if isempty(opts.ShowValues)
    % Auto-determine based on matrix size
    opts.ShowValues = numel(displayMat) <= 400;  % 20x20 or smaller
end

%% --- Create figure ---
figHandle = figure('Position', [100, 100, 900, 700], 'Color', 'w');

% Main heatmap
ax1 = subplot(2, 2, [1, 3]);
imagesc(displayMat, 'AlphaData', ~isnan(displayMat));
colormap(ax1, opts.Colormap);
cb = colorbar;
cb.Label.String = colorLabel;
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

% Labels
xlabel('Target (Track Starts / Particles at t+1)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Source (Track Ends / Particles at t)', 'FontSize', 12, 'FontWeight', 'bold');
title(opts.Title, 'FontSize', 14, 'FontWeight', 'bold');

% Set background color for NaN values
set(ax1, 'Color', [0.9 0.9 0.9]);

% Add grid for small matrices
if numel(displayMat) <= 400
    ax1.XTick = 0.5:1:size(displayMat, 2)+0.5;
    ax1.YTick = 0.5:1:size(displayMat, 1)+0.5;
    ax1.XTickLabel = {};
    ax1.YTickLabel = {};
    grid(ax1, 'on');
    ax1.GridColor = [0.5 0.5 0.5];
    ax1.GridAlpha = 0.3;
end

% Show values if requested
if opts.ShowValues
    [nR, nC] = size(displayMat);
    for i = 1:nR
        for j = 1:nC
            if ~isnan(displayMat(i, j))
                val = displayMat(i, j);
                if abs(val) < 1
                    txt = sprintf('%.2f', val);
                else
                    txt = sprintf('%.1f', val);
                end
                text(j, i, txt, 'HorizontalAlignment', 'center', ...
                    'FontSize', 8, 'Color', 'k');
            end
        end
    end
end

%% --- Cost distribution histogram ---
ax2 = subplot(2, 2, 2);
if ~isempty(validCosts)
    if opts.LogScale && costData.minCost > 0
        histogram(log10(validCosts), 30, 'FaceColor', [0.2 0.6 0.8], ...
            'EdgeColor', 'w', 'FaceAlpha', 0.8);
        xlabel('log_{10}(Cost)', 'FontSize', 11, 'FontWeight', 'bold');
    else
        histogram(validCosts, 30, 'FaceColor', [0.2 0.6 0.8], ...
            'EdgeColor', 'w', 'FaceAlpha', 0.8);
        xlabel('Cost', 'FontSize', 11, 'FontWeight', 'bold');
    end
    ylabel('Count', 'FontSize', 11, 'FontWeight', 'bold');
    title('Cost Distribution', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add statistics text
    statText = sprintf('N = %d\nMedian = %.2f\nMean = %.2f', ...
        costData.nValidLinks, costData.medianCost, costData.meanCost);
    text(0.95, 0.95, statText, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k');
end

%% --- Statistics summary panel ---
ax3 = subplot(2, 2, 4);
axis off;

% Create text summary
summaryText = {
    sprintf('\\bf{Cost Matrix Statistics}')
    ''
    sprintf('Matrix size: %d × %d', nRows, nCols)
    sprintf('Valid links: %d (%.1f%%)', costData.nValidLinks, 100*(1-costData.sparsity))
    sprintf('Sparsity: %.1f%%', 100*costData.sparsity)
    ''
    sprintf('Cost range: [%.2e, %.2e]', costData.minCost, costData.maxCost)
    sprintf('Mean cost: %.2e', costData.meanCost)
    sprintf('Median cost: %.2e', costData.medianCost)
    sprintf('Std dev: %.2e', costData.stdCost)
};

text(0.1, 0.9, summaryText, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, ...
    'FontName', 'FixedWidth');

%% --- Add 6D-specific visualization if provided ---
if ~isempty(opts.OrientationCosts) && ~isempty(opts.SpatialCosts)
    % Create additional figure for 6D breakdown
    figure('Position', [150, 150, 800, 400], 'Color', 'w');
    
    subplot(1, 2, 1);
    scatter(opts.SpatialCosts(:), opts.OrientationCosts(:), 20, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Spatial Cost', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Orientation Cost', 'FontSize', 11, 'FontWeight', 'bold');
    title('Spatial vs Orientation Costs', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    subplot(1, 2, 2);
    totalCost = opts.SpatialCosts(:) + opts.OrientationCosts(:);
    spatialFrac = opts.SpatialCosts(:) ./ totalCost;
    histogram(spatialFrac, 20, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'w');
    xlabel('Spatial Fraction of Total Cost', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Count', 'FontSize', 11, 'FontWeight', 'bold');
    title('Cost Breakdown', 'FontSize', 12, 'FontWeight', 'bold');
end

%% --- Save if requested ---
if ~isempty(opts.SavePath)
    saveas(figHandle, opts.SavePath);
    fprintf('Cost matrix visualization saved to: %s\n', opts.SavePath);
end

end
