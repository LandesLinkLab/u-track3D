function plot_fig3_draft_v4()
%PLOT_FIG3_DRAFT_V4  Publication-size Figure 3. v4 = v3 cleanup:
%   - Panel labels A/B/C/D moved OUTSIDE each tile (top-left, above axes)
%     so they no longer cover data.
%   - Panel D: x-tick labels rotated 45 deg; cell-value annotations dropped
%     (color carries the signal; the strongest cells are listed in the
%     caption / text). Far cleaner read.
%   - Panel B: explicit log ticks at 1, 10.
%   - Panel A: y-axis tick at xyMax shown explicitly.
%   - Figure height raised to 6.5" (still under the 6.56" hard cap).
%
%  Target: double-column (7.0" wide x 6.4" tall, well under the 6.56"
%  hard ceiling Emil set with the black reference rectangles in slide 8 of
%  the figure-outline pptx).
%
%  Typography (Emil's spec):
%   - Calibri Bold for all text in all panels.
%   - Panel labels A/B/C/D at 24 pt, top-left corner of each subplot.
%   - Axes, tick, and colorbar labels at 14 pt.
%   - No per-panel titles (frees vertical space and reduces clutter).
%
%  Panel C: med-Delta annotations removed from the figure (they go in the
%  caption); only frac6Dbetter +/- bootstrap CI is shown, with the % value
%  printed inside each bar. Y-tick labels shortened.
%
%  Panel D: cell numbers shown where |rho| >= 0.3 (slightly stricter than
%  v2 to reduce visual clutter). Shorter input-parameter abbreviations.
%
%  Output: Manuscripts\6D Tracking Paper_REAL\Figure Outline\fig3_draft_v3.png
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

% Metric, sign (+1 higher better / -1 lower better), short label
metrics = { ...
    'linkAcc',          +1, 'Linking acc.';
    'meanPurity',       +1, 'Purity';
    'nMerged',          -1, 'False merges';
    'oriMedianDeg',     -1, 'Orient. err';
    'meanCompleteness', +1, 'Completeness';
    'JSC',              +1, 'JSC';
    'alpha',            +1, 'MSD \alpha';
};
nM = size(metrics, 1);
deltas = zeros(N, nM);
for k = 1:nM
    deltas(:, k) = metrics{k, 2} * (T6.(metrics{k, 1}) - T3.(metrics{k, 1}));
end

% ---------------- figure ----------------
% Target 7.0" x 6.4" final, which becomes the embedded image size.
fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.0 6.5], ...
             'Color', 'w', 'Visible', 'off');
set(fig, ...
    'DefaultAxesFontName',     'Calibri', ...
    'DefaultAxesFontWeight',   'bold', ...
    'DefaultAxesFontSize',     14, ...
    'DefaultAxesLineWidth',    1.5, ...
    'DefaultAxesTickDir',      'in', ...
    'DefaultAxesTickLength',   [0.02 0.02], ...
    'DefaultTextFontName',     'Calibri', ...
    'DefaultTextFontWeight',   'bold', ...
    'DefaultTextFontSize',     14, ...
    'DefaultColorbarFontName', 'Calibri', ...
    'DefaultColorbarFontWeight','bold', ...
    'DefaultColorbarFontSize', 14, ...
    'DefaultLineLineWidth',    1.8);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

% ---------------- Panel A: linkAcc paired scatter ----------------
ax = nexttile(tl);
hold(ax, 'on');
scatter(ax, 100*T3.linkAcc, 100*T6.linkAcc, 22, T6.NoiseXY, 'filled', ...
        'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.25);
ax.ColorScale = 'log';
xyMax = ceil(max(100*[T3.linkAcc; T6.linkAcc])/10)*10;
plot(ax, [0 xyMax], [0 xyMax], 'k--', 'LineWidth', 1.2);
xlim(ax, [0 xyMax]); ylim(ax, [0 xyMax]);
set(ax, 'XTick', 0:20:xyMax, 'YTick', 0:20:xyMax);
axis(ax, 'square'); box(ax, 'on');
xlabel(ax, '3D linking acc. (%)');
ylabel(ax, '6D linking acc. (%)');
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = '\sigma_{xy} (nm)';
cb.Ticks = [5 10 20 50 100];
cb.TickLabels = {'5','10','20','50','100'};
addPanelLabel(ax, 'A');

% ---------------- Panel B: nMerged paired scatter ----------------
ax = nexttile(tl);
hold(ax, 'on');
x = T3.nMerged + 0.5; y = T6.nMerged + 0.5;
scatter(ax, x, y, 22, T6.NumParticles, 'filled', ...
        'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.25);
xyMax = 10^ceil(log10(max([x; y])));
plot(ax, [0.5 xyMax], [0.5 xyMax], 'k--', 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [0.5 xyMax]); ylim(ax, [0.5 xyMax]);
tickVals = [1 10 100];
tickVals = tickVals(tickVals <= xyMax);
set(ax, 'XTick', tickVals, 'YTick', tickVals);
axis(ax, 'square'); box(ax, 'on');
xlabel(ax, '3D false merges');
ylabel(ax, '6D false merges');
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = 'Particles/frame';
addPanelLabel(ax, 'B');

% ---------------- Panel C: per-metric 6D-better fraction ----------------
ax = nexttile(tl);
hold(ax, 'on');
nBoot = 2000;
rng(20260522);
fracs   = zeros(nM, 1);
ciLo    = zeros(nM, 1);
ciHi    = zeros(nM, 1);
medDelt = zeros(nM, 1);
for k = 1:nM
    d = deltas(:, k);
    fracs(k)   = mean(d > 0);
    medDelt(k) = median(d);
    bs = zeros(nBoot, 1);
    n = length(d);
    for b = 1:nBoot
        bs(b) = mean(d(randi(n, n, 1)) > 0);
    end
    ciLo(k) = prctile(bs, 2.5);
    ciHi(k) = prctile(bs, 97.5);
end

[~, ord] = sort(fracs, 'descend');
fracs   = fracs(ord);
ciLo    = ciLo(ord);
ciHi    = ciHi(ord);
medDelt = medDelt(ord);
labels  = metrics(ord, 3);

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
    % % value inside bar (white) for long bars, just past CI (black) for short
    if fracs(k) >= 0.20
        text(ax, 100*fracs(k) - 3, yPos(k), sprintf('%.0f%%', 100*fracs(k)), ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
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
set(ax, 'YTick', yPos, 'YTickLabel', labels, 'YDir', 'reverse');
xlim(ax, [0 100]);
ylim(ax, [0.4 nM + 0.6]);
xlabel(ax, '% MC samples where 6D > 3D');
box(ax, 'on');
addPanelLabel(ax, 'C');

% ---------------- Panel D: Spearman driver heatmap ----------------
ax = nexttile(tl);
% Keep only the 6 parameters with at least one |rho| >= 0.25 across metrics.
% Excluded (all weak across the 7 metrics): D_T, D_R, gamma, k_off, k_on, S, sigma_O.
% This is disclosed in the caption.
inputs = {'NumParticles', 'NoiseXY', 'NoiseTheta', ...
          'wOrientSampled', 'maxSearchRadius', 'timeWindow'};
prettyInputs = {'\rho_p', '\sigma_{xy}', '\sigma_\theta', ...
                'w_O', 'R', '\Delta t'};
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
        'XTickLabelRotation', 45, 'YTick', 1:nM, 'YTickLabel', labels);
% With only 6 columns we can keep tick labels at the default 14 pt.
% no cell-value annotations: color is the channel; keep ax uncluttered.
% drop xlabel: prettyInputs are already labeled.
cb = colorbar(ax);
cb.Label.String = 'Spearman \rho';
cb.Ticks = [-0.7 -0.35 0 0.35 0.7];
box(ax, 'on');
addPanelLabel(ax, 'D');

% ---------------- save ----------------
outDir = ['C:\Users\emilg\OneDrive - University of Illinois - Urbana\' ...
          'Manuscripts\6D Tracking Paper_REAL\Figure Outline'];
if ~exist(outDir, 'dir'), mkdir(outDir); end
outPath = fullfile(outDir, 'fig3_draft_v4.png');
% Resolution 600 gives 7" x 6.4" -> 4200 x 3840 px
exportgraphics(fig, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
close(fig);
end


function addPanelLabel(ax, letter)
% Panel label OUTSIDE axes (above-left), so it never overlaps data.
text(ax, -0.20, 1.08, letter, 'Units', 'normalized', ...
     'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Calibri', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'Clipping', 'off');
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
