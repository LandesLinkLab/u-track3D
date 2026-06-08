%% replot_benchmark6D.m
%  Re-plot benchmark6D results (no re-run). Adds #merges panel.

clear; clc;
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
outDir = fullfile(thisDir, 'benchmark_output');
load(fullfile(outDir, 'benchmark6D_results.mat'), 'T', 'densities', 'wOrientVals', 'seeds');

%% Build the 5-panel comparison: link accuracy / completeness / purity / splits / merges
fig = figure('Position', [50 50 1600 360], 'Color', 'w', 'Visible', 'off');

metricCols  = {'linkingAccuracy','meanCompleteness','meanPurity','nSplit','nMerged'};
metricLabel = {'Linking accuracy','Mean completeness','Mean purity','# splits','# merges'};
metricIsPct = [true,  true,  true,  false, false];

colors  = lines(length(wOrientVals));
markers = {'o','s','^','d','v','p'};

for iM = 1:5
    subplot(1, 5, iM); hold on; grid on; box on;
    colName = metricCols{iM};
    scale = 100 * metricIsPct(iM) + 1 * (~metricIsPct(iM));
    for iW = 1:length(wOrientVals)
        wO = wOrientVals(iW);
        m = zeros(size(densities)); s = zeros(size(densities));
        for iD = 1:length(densities)
            v = T.(colName)(T.density==densities(iD) & T.wOrient==wO);
            m(iD) = mean(v) * scale; s(iD) = std(v)/sqrt(length(v)) * scale;
        end
        mk = markers{1 + mod(iW-1, length(markers))};
        if wO == 0
            dispName = '3D (wOrient=0)';
        else
            dispName = sprintf('6D wOrient=%g', wO);
        end
        errorbar(densities, m, s, ['-' mk], 'Color', colors(iW,:), ...
            'LineWidth', 2, 'MarkerFaceColor', colors(iW,:), 'MarkerSize', 7, ...
            'DisplayName', dispName);
    end
    xlabel('# particles', 'FontWeight','bold', 'FontSize', 12);
    if metricIsPct(iM)
        ylabel([metricLabel{iM} ' (%)'], 'FontWeight','bold','FontSize', 12);
        ylim([0 105]);
    else
        ylabel(metricLabel{iM}, 'FontWeight','bold','FontSize', 12);
    end
    title(metricLabel{iM}, 'FontSize', 13, 'FontWeight','bold');
    set(gca, 'XTick', densities, 'FontSize', 11, 'LineWidth', 1.5, 'TickDir','in');
    if iM == 1
        legend('Location','southwest','FontSize',9);
    end
end

sgtitle(sprintf('3D vs 6D tracking — heavy blinking, Level 2 Kalman (%d seeds/cond)', ...
    length(seeds)), 'FontSize', 14, 'FontWeight', 'bold');

pngPath = fullfile(outDir, 'benchmark6D_comparison.png');
exportgraphics(fig, pngPath, 'Resolution', 300);
fprintf('Plot saved: %s\n', pngPath);
close(fig);
