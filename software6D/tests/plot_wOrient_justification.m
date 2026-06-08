function plot_wOrient_justification()
%PLOT_WORIENT_JUSTIFICATION  Publication figure: why wOrient = 0.3 is the
%recommended default.
%
%  Pulls the wOrient sweep rows from sweeps_v3/sweeps_results.mat and
%  builds a 2x3 panel layout:
%    Top row:    Linking accuracy, # merges, Median orientation error
%    Bottom row: Mean completeness, Chenouard alpha, Combined score
%
%  A shaded "sweet spot" band marks 0.1 <= wOrient <= 0.5, with a bold
%  vertical line at wOrient=0.3 (recommended default). w=0 (3D baseline)
%  is shown as a star on each panel for reference.
%
%  Output:
%    tests/benchmark_output/sweeps_v3/wOrient_justification.png  (600 DPI)
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(genpath(fullfile(thisDir, '..')));

dataPath = fullfile(thisDir, 'benchmark_output', 'sweeps_v3', 'sweeps_results.mat');
S = load(dataPath, 'T', 'B', 'seeds');
T = S.T; B = S.B; seeds = S.seeds;

%% Extract wOrient sweep rows
mask = strcmp(T.sweep, 'wOrient');
Tw = T(mask, :);
wVals = unique(Tw.axisValue);

[mAcc, sAcc] = aggregate(Tw, 'linkAcc',          wVals);
[mMrg, sMrg] = aggregate(Tw, 'nMerged',          wVals);
[mCmp, sCmp] = aggregate(Tw, 'meanCompleteness', wVals);
[mAlf, sAlf] = aggregate(Tw, 'alpha',            wVals);
[mOri, sOri] = aggregate(Tw, 'oriMedianDeg',     wVals);

%% Plot setup
fig = figure('Position', [60 60 1500 820], 'Color', 'w', 'Visible', 'off');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

sweetSpot = [0.1, 0.5];
recW      = 0.3;
posWMin   = max(0.05, min(wVals(wVals > 0)));
xLog      = [posWMin * 0.7, max(wVals) * 1.5];

% For log scale, replace w=0 with a small positive value so it can be plotted
wPlot = wVals; wPlot(wPlot == 0) = xLog(1) * 1.4;
i0 = find(wVals == 0, 1);

%% Panel 1: Linking accuracy (WIN)
nexttile;
errorbar(wPlot, 100*mAcc, 100*sAcc, '-o', 'Color', [0.20 0.55 0.20], ...
    'LineWidth', 2.2, 'MarkerFaceColor', [0.20 0.55 0.20], 'MarkerSize', 9);
markBaseline(wPlot, 100*mAcc, i0, -1.5);
ylim([min(100*mAcc) - 2, max(100*mAcc) + 2]);
decorate(xLog, sweetSpot, recW, 'Linking accuracy (%)', 'Linking accuracy');

%% Panel 2: # merges (WIN — lower is better)
nexttile;
errorbar(wPlot, mMrg, sMrg, '-s', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 2.2, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 9);
markBaseline(wPlot, mMrg, i0, +1.5);
ylim([min(mMrg) * 0.85, max(mMrg) * 1.15]);
decorate(xLog, sweetSpot, recW, '# merges (lower better)', 'Wrong-particle fusions');

%% Panel 3: Median orientation error (NEUTRAL — should NOT change)
nexttile;
errorbar(wPlot, mOri, sOri, '-d', 'Color', [0.40 0.40 0.85], ...
    'LineWidth', 2.2, 'MarkerFaceColor', [0.40 0.40 0.85], 'MarkerSize', 9);
markBaseline(wPlot, mOri, i0, -0.10);
yPad = 0.5;
ylim([min(mOri) - yPad, max(mOri) + yPad]);
decorate(xLog, sweetSpot, recW, 'Median orient. error (deg)', 'Orientation fidelity');

%% Panel 4: Mean completeness (COST — drops at high w)
nexttile;
errorbar(wPlot, 100*mCmp, 100*sCmp, '-^', 'Color', [0.65 0.20 0.65], ...
    'LineWidth', 2.2, 'MarkerFaceColor', [0.65 0.20 0.65], 'MarkerSize', 9);
markBaseline(wPlot, 100*mCmp, i0, +2);
ylim([min(100*mCmp) - 3, max(100*mCmp) + 3]);
decorate(xLog, sweetSpot, recW, 'Mean completeness (%)', 'Track completeness');

%% Panel 5: Chenouard alpha (COST — drops at high w)
nexttile;
errorbar(wPlot, mAlf, sAlf, '-v', 'Color', [0.85 0.50 0.10], ...
    'LineWidth', 2.2, 'MarkerFaceColor', [0.85 0.50 0.10], 'MarkerSize', 9);
markBaseline(wPlot, mAlf, i0, -0.025);
ylim([0, max(mAlf) * 1.15]);
decorate(xLog, sweetSpot, recW, 'Chenouard \alpha', 'Chenouard \alpha (overall score)');

%% Panel 6: Combined improvement score vs 3D baseline
nexttile;
ref = i0; if isempty(ref), ref = 1; end
gainAcc = (mAcc - mAcc(ref)) / max(mAcc(ref), eps);    % up = good
gainMrg = (mMrg(ref) - mMrg) / max(mMrg(ref), eps);    % down = good
gainCmp = (mCmp - mCmp(ref)) / max(mCmp(ref), eps);    % up = good
gainAlf = (mAlf - mAlf(ref)) / max(mAlf(ref), eps);    % up = good
combined = (gainAcc + gainMrg + gainCmp + gainAlf) / 4;

bar(wPlot, combined, 'FaceColor', [0.30 0.55 0.85], 'EdgeColor', 'k', 'LineWidth', 1.3);
ylim([min(combined) - 0.05, max(combined) + 0.05]);

[~, iRec] = min(abs(wVals - recW));
text(wPlot(iRec), combined(iRec) + 0.015, 'best balance', ...
    'HorizontalAlignment','center', 'FontWeight','bold', ...
    'FontName','Calibri', 'FontSize', 10, 'Color', [0.10 0.35 0.10]);

decorate(xLog, sweetSpot, recW, 'Mean rel. improvement vs 3D', ...
    'Combined improvement score');

%% Title
title(tl, sprintf( ...
    ['Choosing wOrient: gains saturate by w \\approx 0.3, costs rise above w \\approx 1   ' ...
     '(density=%d, %d seeds)'], B.density, length(seeds)), ...
    'FontName','Calibri', 'FontWeight','bold', 'FontSize', 15);

%% Save
outPng = fullfile(thisDir, 'benchmark_output', 'sweeps_v3', 'wOrient_justification.png');
exportgraphics(fig, outPng, 'Resolution', 600);
close(fig);
fprintf('Saved: %s\n', outPng);
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function [m, s] = aggregate(Tw, col, wVals)
    m = zeros(size(wVals));
    s = zeros(size(wVals));
    for i = 1:length(wVals)
        v = Tw.(col)(Tw.axisValue == wVals(i));
        m(i) = mean(v);
        s(i) = std(v) / sqrt(max(length(v), 1));
    end
end

function decorate(xLog, sweetSpot, recW, ylab, ttl)
    yl = ylim;
    % Sweet-spot shaded band
    patch([sweetSpot(1) sweetSpot(2) sweetSpot(2) sweetSpot(1)], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.96 0.96 0.70], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    % Bring data on top
    ch = get(gca, 'Children');
    set(gca, 'Children', flipud(ch));
    % Vertical line at recommended w
    xline(recW, 'k--', 'LineWidth', 2, 'Alpha', 0.7, ...
        'Label', ' recommended w=0.3', 'LabelOrientation','horizontal', ...
        'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','top', ...
        'FontWeight','bold', 'FontName','Calibri', 'FontSize', 10, ...
        'HandleVisibility', 'off');
    xlim(xLog); ylim(yl);
    set(gca, 'XScale', 'log', 'FontName','Calibri', 'FontWeight','bold', ...
        'FontSize', 12, 'LineWidth', 2, 'TickDir','in', 'Box', 'on');
    xlabel('wOrient', 'FontWeight','bold','FontSize', 13);
    ylabel(ylab, 'FontWeight','bold','FontSize', 13);
    title(ttl, 'FontWeight','bold','FontSize', 14);
    grid on; box on;
end

function markBaseline(wPlot, y, i0, txtOffset)
    if isempty(i0), return; end
    hold on;
    plot(wPlot(i0), y(i0), 'p', 'MarkerSize', 16, ...
        'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(wPlot(i0), y(i0) + txtOffset, '3D', 'HorizontalAlignment','center', ...
        'FontWeight','bold', 'FontName','Calibri', 'FontSize', 10);
end
