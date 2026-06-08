function plot_fig3_draft_v7()
%PLOT_FIG3_DRAFT_V7  Same as v6 but relabel the alpha metric as
%  "Chenouard alpha" (track-overlap score from Chenouard et al. 2014)
%  not "MSD alpha" (which would be the diffusion-exponent power law).
%
%  v6 swaps tiledlayout for direct axes('Position', ...) calls so the four
%  plot boxes are exactly the same size (square), evenly spaced, and the
%  whitespace between rows is small. Panel labels live in the figure
%  margin via annotation(), so they never collide with tick labels or
%  axis labels (this was the v5 failure mode).
%
%  Layout (figure 7.0" wide x 6.5" tall):
%    - Square plot box on every panel: 2.55" x 2.55"
%    - Colorbar slot to the right of each plot: 0.20"
%    - Inter-column gap: 0.50" (for colorbar + label + breathing room)
%    - Inter-row gap:    ~0.45"
%
%  All text: Calibri Bold. Axes/ticks/colorbar 14 pt; panel labels 24 pt.
%
%  v6 also restores BlinkOffRate, GammaMean, orientCRLB to Panel D (v5
%  had dropped them; they actually have |rho| >= 0.25 on at least one
%  metric so they belong in the figure, not the SI).
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
csvPath = fullfile(thisDir, 'benchmark_output', 'mc', 'mc_results.csv');
T = readtable(csvPath);
T3 = T(strcmp(T.mode, '3D'), :);
T6 = T(strcmp(T.mode, '6D'), :);
[~, i3, i6] = intersect(T3.sample, T6.sample);
T3 = T3(i3, :); T6 = T6(i6, :);
N  = height(T3);
fprintf('Paired samples loaded: %d\n', N);

% metrics: name, sign convention (+1 higher better, -1 lower better), short label
metrics = { ...
    'linkAcc',          +1, 'Linking acc.';
    'meanPurity',       +1, 'Purity';
    'nMerged',          -1, 'False merges';
    'oriMedianDeg',     -1, 'Orient. err';
    'meanCompleteness', +1, 'Completeness';
    'JSC',              +1, 'JSC';
    'alpha',            +1, 'Chenouard \alpha';
};
nM = size(metrics, 1);
deltas = zeros(N, nM);
for k = 1:nM
    deltas(:, k) = metrics{k, 2} * (T6.(metrics{k, 1}) - T3.(metrics{k, 1}));
end

% ---------------- figure ----------------
figW = 7.0; figH = 6.5;
fig = figure('Units', 'inches', 'Position', [0.5 0.5 figW figH], ...
             'Color', 'w', 'Visible', 'off');

% Manual layout (figure-normalized). Margins were tuned by checking that:
%   - "Completeness" (longest Panel C y-tick label) fits in left margin
%   - A's colorbar label ("sigma_xy (nm)") + B's y-label both fit in mid gap
%   - B's and D's colorbar labels fit in right margin
plotIn = 1.95;                       % plot side, inches (square)
plotW  = plotIn/figW;
plotH  = plotIn/figH;
xL1   = 1.10/figW;                   % left margin
xL2   = (1.10 + 1.95 + 1.25)/figW;   % col1 + col-gap of 1.25"
% row 1 (top) bottom-edge: figure_top - top_margin - plot_height
yB2   = (figH - 0.45 - 2.55)/figH;
% row 2 (bottom) bottom-edge: 0.55" from bottom for x-axis label
yB1   = 0.55/figH;

posA = [xL1, yB2, plotW, plotH];
posB = [xL2, yB2, plotW, plotH];
posC = [xL1, yB1, plotW, plotH];
posD = [xL2, yB1, plotW, plotH];


% common defaults
set(fig, ...
    'DefaultAxesFontName',      'Calibri', ...
    'DefaultAxesFontWeight',    'bold', ...
    'DefaultAxesFontSize',      14, ...
    'DefaultAxesLineWidth',     1.5, ...
    'DefaultAxesTickDir',       'in', ...
    'DefaultAxesTickLength',    [0.020 0.020], ...
    'DefaultTextFontName',      'Calibri', ...
    'DefaultTextFontWeight',    'bold', ...
    'DefaultTextFontSize',      14, ...
    'DefaultColorbarFontName',  'Calibri', ...
    'DefaultColorbarFontWeight','bold', ...
    'DefaultColorbarFontSize',  14, ...
    'DefaultLineLineWidth',     1.8);

% ---------------- Panel A: linkAcc paired scatter ----------------
ax = axes(fig, 'Units', 'normalized', 'Position', posA);
hold(ax, 'on');
scatter(ax, 100*T3.linkAcc, 100*T6.linkAcc, 22, T6.NoiseXY, 'filled', ...
        'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.25);
ax.ColorScale = 'log';
xyMax = ceil(max(100*[T3.linkAcc; T6.linkAcc])/10)*10;
plot(ax, [0 xyMax], [0 xyMax], 'k--', 'LineWidth', 1.2);
xlim(ax, [0 xyMax]); ylim(ax, [0 xyMax]);
set(ax, 'XTick', [0 round(xyMax/2) xyMax], 'YTick', [0 round(xyMax/2) xyMax]);
xlabel(ax, '3D linking acc. (%)');
ylabel(ax, '6D linking acc. (%)');
colormap(ax, parula);
cb = colorbar(ax, 'eastoutside');
cb.Label.String = '\sigma_{xy} (nm)';
cb.Ticks = [5 10 20 50 100];
cb.TickLabels = {'5','10','20','50','100'};
box(ax, 'on');
ax.Position = posA;                                  % restore axes size
cb.Position = [posA(1)+posA(3)+0.005, posA(2), 0.015, posA(4)];
addPanelLabel(fig, posA, 'A');

% ---------------- Panel B: nMerged paired scatter ----------------
ax = axes(fig, 'Units', 'normalized', 'Position', posB);
hold(ax, 'on');
x = T3.nMerged + 0.5; y = T6.nMerged + 0.5;
scatter(ax, x, y, 22, T6.NumParticles, 'filled', ...
        'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.25);
xyMax = 10^ceil(log10(max([x; y])));
plot(ax, [0.5 xyMax], [0.5 xyMax], 'k--', 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [0.5 xyMax]); ylim(ax, [0.5 xyMax]);
tickVals = [1 10 100]; tickVals = tickVals(tickVals <= xyMax);
set(ax, 'XTick', tickVals, 'YTick', tickVals);
xlabel(ax, '3D false merges');
ylabel(ax, '6D false merges');
colormap(ax, parula);
cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'Particles/frame';
box(ax, 'on');
ax.Position = posB;
cb.Position = [posB(1)+posB(3)+0.005, posB(2), 0.015, posB(4)];
addPanelLabel(fig, posB, 'B');

% ---------------- Panel C: per-metric 6D-better fraction ----------------
ax = axes(fig, 'Units', 'normalized', 'Position', posC);
hold(ax, 'on');
nBoot = 2000;
rng(20260522);
fracs   = zeros(nM, 1);
ciLo    = zeros(nM, 1);
ciHi    = zeros(nM, 1);
for k = 1:nM
    d = deltas(:, k);
    fracs(k) = mean(d > 0);
    bs = zeros(nBoot, 1);
    n = length(d);
    for b = 1:nBoot
        bs(b) = mean(d(randi(n, n, 1)) > 0);
    end
    ciLo(k) = prctile(bs, 2.5);
    ciHi(k) = prctile(bs, 97.5);
end
[~, ord] = sort(fracs, 'descend');
fracs  = fracs(ord); ciLo = ciLo(ord); ciHi = ciHi(ord);
labels = metrics(ord, 3);

yPos = (1:nM)';
greenC = [0.20 0.58 0.27];
redC   = [0.80 0.27 0.27];
for k = 1:nM
    col = greenC; if fracs(k) < 0.5, col = redC; end
    barh(ax, yPos(k), 100*fracs(k), 'FaceColor', col, ...
        'EdgeColor', 'none', 'BarWidth', 0.7);
    errorbar(ax, 100*fracs(k), yPos(k), ...
             100*(fracs(k) - ciLo(k)), 100*(ciHi(k) - fracs(k)), ...
             'horizontal', 'k', 'LineStyle', 'none', ...
             'LineWidth', 1.3, 'CapSize', 6);
    if fracs(k) >= 0.20
        % LEFT-aligned inside bar so the bar-end errorbar cap is unobstructed
        text(ax, 2, yPos(k), sprintf('%.0f%%', 100*fracs(k)), ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
             'FontName', 'Calibri');
    else
        text(ax, 100*ciHi(k) + 2, yPos(k), sprintf('%.0f%%', 100*fracs(k)), ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', ...
             'FontName', 'Calibri');
    end
end
xline(ax, 50, 'k:', 'LineWidth', 1.3);
set(ax, 'YTick', yPos, 'YTickLabel', labels, 'YDir', 'reverse', ...
        'XTick', [0 50 100]);
xlim(ax, [0 100]);
ylim(ax, [0.4 nM + 0.6]);
xlabel(ax, '% MC samples where 6D > 3D');
box(ax, 'on');
ax.Position = posC;
addPanelLabel(fig, posC, 'C');

% ---------------- Panel D: Spearman driver heatmap ----------------
ax = axes(fig, 'Units', 'normalized', 'Position', posD);
% All 8 inputs that show |rho| >= 0.25 on at least one metric.
inputs = {'NumParticles', 'GammaMean', 'BlinkOffRate', 'NoiseXY', 'NoiseTheta', ...
          'wOrientSampled', 'maxSearchRadius', 'orientCRLB'};
prettyInputs = {'\rho_p', '\gamma', 'k_{off}', '\sigma_{xy}', '\sigma_\theta', ...
                'w_O', 'R', '\sigma_O'};
nI = length(inputs);
rhoMat = zeros(nM, nI);
for k = 1:nM
    for j = 1:nI
        rhoMat(k, j) = corr(T6.(inputs{j}), deltas(:, ord(k)), 'Type', 'Spearman');
    end
end
imagesc(ax, rhoMat);
colormap(ax, divergingColormap(64));
clim(ax, [-0.7 0.7]);
set(ax, 'XTick', 1:nI, 'XTickLabel', prettyInputs, ...
        'XTickLabelRotation', 0, 'YTick', 1:nM, 'YTickLabel', labels);
cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'Spearman \rho';
cb.Ticks = [-0.7 -0.35 0 0.35 0.7];
box(ax, 'on');
ax.Position = posD;
cb.Position = [posD(1)+posD(3)+0.005, posD(2), 0.015, posD(4)];
addPanelLabel(fig, posD, 'D');

% ---------------- save ----------------
outDir = ['C:\Users\emilg\OneDrive - University of Illinois - Urbana\' ...
          'Manuscripts\6D Tracking Paper_REAL\Figure Outline'];
if ~exist(outDir, 'dir'), mkdir(outDir); end
outPath = fullfile(outDir, 'fig3_draft_v7.png');
exportgraphics(fig, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
close(fig);
end


function addPanelLabel(fig, pos, letter)
% Place panel label ABOVE the plot box (not aligned with axes top),
% so it never overlaps the topmost tick value.
labelW = 0.05; labelH = 0.05;
labelX = pos(1) - 0.075;
labelY = pos(2) + pos(4) + 0.005;   % start just above axes top
annotation(fig, 'textbox', [labelX, labelY, labelW, labelH], ...
    'String', letter, ...
    'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Calibri', ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end


function cmap = divergingColormap(n)
if nargin < 1, n = 64; end
half = floor(n/2);
t = linspace(0, 1, half)';
top    = [t, t, ones(half, 1)];
bottom = [ones(half,1), 1-t, 1-t];
if mod(n, 2) == 1
    cmap = [top; [1 1 1]; bottom];
else
    cmap = [top; bottom];
end
cmap = cmap(1:n, :);
end
