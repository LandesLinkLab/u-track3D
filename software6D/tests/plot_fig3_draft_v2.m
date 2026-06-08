function plot_fig3_draft_v2()
%PLOT_FIG3_DRAFT_V2  Figure 3 algorithm benchmarking, draft v2.
%
%  Refinements vs v1:
%   - Panel C xlim tightened to [0, 100]; med-delta annotations placed
%     INSIDE each bar (white text) when the bar is long enough to hold them,
%     otherwise just outside the CI cap (black text). No more wasted whitespace.
%   - Per-metric labels left untouched; definitions live in the slide-side text
%     boxes and the figure caption.
%
%  Output: Manuscripts\6D Tracking Paper_REAL\Figure Outline\fig3_draft_v2.png
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

% metric, sign convention (+1: higher better, -1: lower better), label
metrics = { ...
    'linkAcc',          +1, 'Linking accuracy';
    'meanPurity',       +1, 'Track purity';
    'nMerged',          -1, 'Wrong-particle fusions';
    'oriMedianDeg',     -1, 'Orientation error';
    'meanCompleteness', +1, 'Track completeness';
    'JSC',              +1, 'Jaccard similarity';
    'alpha',            +1, 'MSD exponent \alpha'; ...
};
nM = size(metrics, 1);
deltas = zeros(N, nM);
for k = 1:nM
    deltas(:, k) = metrics{k, 2} * (T6.(metrics{k, 1}) - T3.(metrics{k, 1}));
end

% ---------------- figure ----------------
fig = figure('Position', [40 40 1600 1180], 'Color', 'w', 'Visible', 'off');
set(fig, 'DefaultAxesFontName', 'Calibri', 'DefaultAxesFontWeight', 'bold', ...
         'DefaultAxesFontSize', 11, ...
         'DefaultAxesLineWidth', 1.5, 'DefaultAxesTickDir', 'in', ...
         'DefaultLineLineWidth', 2);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------------- Panel A: linkAcc paired scatter ----------------
ax = nexttile(tl);
hold(ax, 'on');
scatter(ax, 100*T3.linkAcc, 100*T6.linkAcc, 36, T6.NoiseXY, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.3);
ax.ColorScale = 'log';
xyMax = max(100*[T3.linkAcc; T6.linkAcc]) * 1.05;
plot(ax, [0 xyMax], [0 xyMax], 'k--', 'LineWidth', 1.3);
xlim(ax, [0 xyMax]); ylim(ax, [0 xyMax]);
axis(ax, 'square'); box(ax, 'on');
xlabel(ax, '3D linking accuracy (%)');
ylabel(ax, '6D linking accuracy (%)');
fracW = 100*mean(T6.linkAcc > T3.linkAcc);
title(ax, sprintf('A. Linking accuracy  (6D > 3D on %.0f%% of samples)', fracW));
cb = colorbar(ax); cb.Label.String = 'NoiseXY (nm)';
cb.Label.FontWeight = 'bold';
colormap(ax, parula);

% ---------------- Panel B: nMerged paired scatter ----------------
ax = nexttile(tl);
hold(ax, 'on');
x = T3.nMerged + 0.5; y = T6.nMerged + 0.5;
scatter(ax, x, y, 36, T6.NumParticles, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.3);
xyMax = max([x; y]) * 1.1;
plot(ax, [0.5 xyMax], [0.5 xyMax], 'k--', 'LineWidth', 1.3);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [0.5 xyMax]); ylim(ax, [0.5 xyMax]);
axis(ax, 'square'); box(ax, 'on');
xlabel(ax, '3D wrong-particle fusions');
ylabel(ax, '6D wrong-particle fusions');
fracB = 100*mean(T6.nMerged < T3.nMerged);
title(ax, sprintf('B. Wrong-particle fusions  (6D < 3D on %.0f%% of samples)', fracB));
cb = colorbar(ax); cb.Label.String = 'Particle density';
cb.Label.FontWeight = 'bold';
colormap(ax, parula);

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

% sort metrics by frac so the figure reads top-down "where 6D wins"
[~, ord] = sort(fracs, 'descend');
fracs    = fracs(ord);
ciLo     = ciLo(ord);
ciHi     = ciHi(ord);
medDelt  = medDelt(ord);
labels   = metrics(ord, 3);

yPos = (1:nM)';
for k = 1:nM
    if fracs(k) >= 0.5
        col = [0.20 0.55 0.20];   % green = 6D-better majority
    else
        col = [0.80 0.30 0.30];   % red   = 3D-better majority
    end
    barh(ax, yPos(k), 100*fracs(k), 'FaceColor', col, ...
        'EdgeColor', 'none', 'BarWidth', 0.62);
    errorbar(ax, 100*fracs(k), yPos(k), ...
             100*(fracs(k) - ciLo(k)), 100*(ciHi(k) - fracs(k)), ...
             'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 1.4, ...
             'CapSize', 8);
end
xline(ax, 50, 'k:', 'LineWidth', 1.5);

% smart med-delta placement: inside bar (white) for long bars, else outside CI cap (black)
for k = 1:nM
    msg = sprintf('med \\Delta = %s', fmtDelta(medDelt(k), labels{k}));
    if fracs(k) >= 0.30
        text(ax, 2, yPos(k), msg, ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10, ...
             'FontName', 'Calibri');
    else
        text(ax, 100*ciHi(k) + 2, yPos(k), msg, ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'k', 'FontWeight', 'bold', 'FontSize', 10, ...
             'FontName', 'Calibri');
    end
end

set(ax, 'YTick', yPos, 'YTickLabel', labels, 'YDir', 'reverse');
xlim(ax, [0 100]);
ylim(ax, [0.4 nM + 0.6]);
xlabel(ax, '% of MC samples where 6D > 3D (bootstrap 95% CI)');
title(ax, 'C. Where does the orientation channel help?');
box(ax, 'on');

% ---------------- Panel D: Spearman driver heatmap ----------------
ax = nexttile(tl);
inputs = {'NumParticles', 'Dtrans', 'Drot', 'GammaMean', 'BlinkOffRate', ...
          'BlinkOnRate', 'Signal', 'NoiseXY', 'NoiseTheta', ...
          'wOrientSampled', 'maxSearchRadius', 'timeWindow', 'orientCRLB'};
prettyInputs = {'density', 'D_{trans}', 'D_{rot}', '\gamma', 'blink off', ...
                'blink on', 'signal', 'noise XY', 'noise \theta', ...
                'w_{Orient}', 'search R', 't window', 'orient CRLB'};
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
title(ax, 'D. What drives the 6D advantage  (Spearman \rho)');
cb = colorbar(ax); cb.Label.String = '\rho  (input vs sign-aligned \Delta)';
cb.Label.FontWeight = 'bold';
for k = 1:nM
    for j = 1:nI
        if abs(rhoMat(k, j)) >= 0.25
            txtCol = [1 1 1]; if abs(rhoMat(k,j)) < 0.4, txtCol = [0.1 0.1 0.1]; end
            text(ax, j, k, sprintf('%.2f', rhoMat(k, j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8.5, ...
                'FontWeight', 'bold', 'Color', txtCol);
        end
    end
end
box(ax, 'on');
axis(ax, 'tight');

% ---------------- title + save ----------------
title(tl, sprintf(['Figure 3 (draft v2). Orientation-aware 6D tracking: ' ...
    'where it helps and where it costs.  ' ...
    'N = %d paired Monte-Carlo samples (15-dim parameter space).'], N), ...
    'FontName', 'Calibri', 'FontWeight', 'bold', 'FontSize', 13);

outDir = ['C:\Users\emilg\OneDrive - University of Illinois - Urbana\' ...
          'Manuscripts\6D Tracking Paper_REAL\Figure Outline'];
if ~exist(outDir, 'dir'), mkdir(outDir); end
outPath = fullfile(outDir, 'fig3_draft_v2.png');
exportgraphics(fig, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
close(fig);
end


function s = fmtDelta(v, label)
%FMTDELTA  Format Delta with sensible precision and units.
if contains(label, 'accuracy') || contains(label, 'purity') || ...
   contains(label, 'completeness') || contains(label, 'Jaccard') || ...
   contains(label, 'alpha')
    s = sprintf('%+.2f pp', 100*v);
elseif contains(label, 'Orientation')
    s = sprintf('%+.2f\\circ', v);
elseif contains(label, 'fusions')
    s = sprintf('%+.1f', v);
else
    s = sprintf('%+.3g', v);
end
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
