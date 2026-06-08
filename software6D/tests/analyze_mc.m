function analyze_mc()
%ANALYZE_MC  Post-hoc analysis of the LHS Monte Carlo benchmark.
%
%  Loads tests/benchmark_output/mc/mc_results.mat and produces:
%    1. mc_sensitivity.png — rank-correlation sensitivity per metric per
%       input variable (which inputs predict the 6D advantage?)
%    2. mc_partial_dependence.png — small-multiples of metric Δ (6D−3D)
%       vs each input variable, smoothed
%    3. mc_pair_scatter.png — 2D scatter of (wOrient, NumParticles) etc.
%       colored by 6D advantage
%    4. mc_summary.csv — table of overall and per-quartile 6D advantage
%    5. console summary printout
%
%  Emil Gillett, Landes Research Group, UIUC, 2026.

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
outDir = fullfile(thisDir, 'benchmark_output', 'mc');
S = load(fullfile(outDir, 'mc_results.mat'), 'T', 'P', 'F', 'Tdesign', 'nSamples');
T = S.T; P = S.P;

%% Split into 3D and 6D and merge on sample
is3 = strcmp(T.mode, '3D');
is6 = strcmp(T.mode, '6D');
T3 = T(is3, :);
T6 = T(is6, :);
[~, i3, i6] = intersect(T3.sample, T6.sample);
T3 = T3(i3, :); T6 = T6(i6, :);
fprintf('Loaded %d paired samples (3D vs 6D) from %d total.\n', height(T3), S.nSamples);

%% Inputs (sampled vars) and outputs (delta metrics)
inputCols = {'NumParticles','Dtrans','Drot','GammaMean','BlinkOffRate', ...
             'BlinkOnRate','Signal','NoiseXY','NoiseTheta','wOrientSampled', ...
             'maxSearchRadius','brownStdMult','timeWindow', ...
             'orientKalmanLevel','orientCRLB'};
inputDisp = {'density','Dtrans','Drot','GammaMean','BlinkOffRate', ...
             'BlinkOnRate','Signal','NoiseXY','NoiseTheta','wOrient', ...
             'maxSearchRadius','brownStdMult','timeWindow', ...
             'KalmanLevel','orientCRLB'};

% Output metrics (delta = 6D - 3D; for # merges and orient err, NEGATIVE delta = 6D better)
outputCols  = {'linkAcc','meanCompleteness','meanPurity','nMerged','alpha','JSC','oriMedianDeg'};
outputDisp  = {'\Delta linkAcc','\Delta completeness','\Delta purity', ...
               '\Delta # merges','\Delta \alpha','\Delta JSC','\Delta oriErr (deg)'};
outputDir   = [+1, +1, +1, -1, +1, +1, -1];  % +1: 6D better when delta positive; -1: 6D better when delta negative

% Build Δ table
nP = height(T3);
X = nan(nP, length(inputCols));
Y = nan(nP, length(outputCols));
for v = 1:length(inputCols)
    X(:, v) = T6.(inputCols{v});  % design inputs are identical between modes
end
for k = 1:length(outputCols)
    Y(:, k) = T6.(outputCols{k}) - T3.(outputCols{k});
end

%% ========================================================================
%  PLOT 1: Sensitivity heatmap (Spearman rank correlation between input and Δmetric)
%  ========================================================================
S_rho = nan(length(inputCols), length(outputCols));
S_p   = nan(length(inputCols), length(outputCols));
for v = 1:length(inputCols)
    for k = 1:length(outputCols)
        ok = ~isnan(X(:,v)) & ~isnan(Y(:,k));
        if sum(ok) >= 5
            [r, p] = corr(X(ok, v), Y(ok, k), 'Type', 'Spearman');
        else
            r = NaN; p = NaN;
        end
        S_rho(v, k) = r;
        S_p(v, k)   = p;
    end
end

% Re-orient so that LARGE positive value = "input increases 6D advantage"
% For metrics where negative Δ is "6D better", flip sign.
S_signed = S_rho .* repmat(outputDir, length(inputCols), 1);

fig1 = figure('Position', [50 50 1100 700], 'Color', 'w', 'Visible', 'off');
imagesc(S_signed);
colormap(redblueMap()); cb = colorbar; cb.Label.String = 'Spearman \rho (signed: + = 6D advantage grows)';
cb.Label.FontWeight = 'bold';
caxis([-0.6 0.6]);

set(gca, 'XTick', 1:length(outputCols), 'XTickLabel', outputDisp, ...
    'YTick', 1:length(inputCols), 'YTickLabel', inputDisp, ...
    'FontName','Calibri', 'FontWeight','bold','FontSize', 11, 'TickDir','out', 'LineWidth', 1.5);
xtickangle(30);

% Annotate significant cells
for v = 1:length(inputCols)
    for k = 1:length(outputCols)
        if ~isnan(S_p(v, k)) && S_p(v, k) < 0.01 && abs(S_signed(v, k)) > 0.10
            txtCol = 'w'; if abs(S_signed(v, k)) < 0.3, txtCol = 'k'; end
            text(k, v, sprintf('%.2f', S_signed(v, k)), ...
                'HorizontalAlignment','center', 'FontWeight','bold', ...
                'FontName','Calibri', 'FontSize', 10, 'Color', txtCol);
        end
    end
end

title({'MC sensitivity: rank-correlation of each input with each metric''s 6D advantage', ...
       sprintf('Red = increasing this input STRENGTHENS the 6D advantage  (N=%d)', nP)}, ...
    'FontWeight','bold', 'FontName','Calibri', 'FontSize', 13);
xlabel('Δ Metric (6D − 3D, signed for "6D better")', 'FontWeight','bold');
ylabel('Input variable', 'FontWeight','bold');

png1 = fullfile(outDir, 'mc_sensitivity.png');
exportgraphics(fig1, png1, 'Resolution', 600); close(fig1);
fprintf('Saved: %s\n', png1);

%% ========================================================================
%  PLOT 2: Partial dependence — for each (input, metric) bin the input
%  values into 8 quantiles and plot the mean Δmetric per bin
%  ========================================================================
nBins = 8;
fig2 = figure('Position', [40 40 1900 1100], 'Color', 'w', 'Visible', 'off');
tl2 = tiledlayout(length(outputCols), length(inputCols), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:length(outputCols)
    for v = 1:length(inputCols)
        nexttile;
        xv = X(:, v); yv = Y(:, k) * outputDir(k);  % flip so + = 6D better
        ok = ~isnan(xv) & ~isnan(yv);
        if sum(ok) < 10
            text(0.5, 0.5, 'insufficient', 'HorizontalAlignment','center', ...
                'Units','normalized', 'FontName','Calibri');
            axis off; continue;
        end
        xv = xv(ok); yv = yv(ok);
        % Quantile bins
        edges = quantile(xv, linspace(0, 1, nBins + 1));
        edges(1) = min(xv) - eps; edges(end) = max(xv) + eps;
        binIdx = discretize(xv, edges);
        m = nan(nBins, 1); s = nan(nBins, 1); xb = nan(nBins, 1);
        for b = 1:nBins
            mask = binIdx == b;
            if any(mask)
                m(b)  = mean(yv(mask));
                s(b)  = std(yv(mask)) / sqrt(sum(mask));
                xb(b) = mean(xv(mask));
            end
        end

        hold on;
        yline(0, 'k-', 'Alpha', 0.4, 'LineWidth', 1);
        % light raw scatter
        scatter(xv, yv, 6, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.25);
        % binned mean with SE
        errorbar(xb, m, s, '-o', 'Color', [0.85 0.33 0.10], ...
            'LineWidth', 2, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 6);
        box on; grid on;
        set(gca, 'FontName','Calibri','FontSize', 9, 'LineWidth', 1.0, 'TickDir', 'in');

        % log x where appropriate
        pv = P{v};
        if strcmp(pv{4}, 'log')
            set(gca, 'XScale', 'log');
        end

        if k == length(outputCols)
            xlabel(inputDisp{v}, 'FontWeight','bold','FontSize', 9);
        end
        if v == 1
            ylabel(outputDisp{k}, 'FontWeight','bold','FontSize', 9);
        end
    end
end
title(tl2, sprintf( ...
    'Partial dependence: mean (6D − 3D) per input bin (positive = 6D better; N=%d)', nP), ...
    'FontName','Calibri', 'FontWeight','bold', 'FontSize', 14);

png2 = fullfile(outDir, 'mc_partial_dependence.png');
exportgraphics(fig2, png2, 'Resolution', 600); close(fig2);
fprintf('Saved: %s\n', png2);

%% ========================================================================
%  PLOT 3: 2D scatter of (wOrient, density) colored by 6D advantage on linkAcc
%  ========================================================================
fig3 = figure('Position', [50 50 1450 450], 'Color', 'w', 'Visible', 'off');
tl3 = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

panelInputs  = {'wOrientSampled', 'BlinkOffRate', 'maxSearchRadius'};
panelLabels  = {'wOrient', 'BlinkOffRate', 'maxSearchRadius'};
panelLogX    = [true, false, true];

% Color = 6D advantage on linkAcc (positive = 6D wins)
cAdv = (T6.linkAcc - T3.linkAcc);
density = T6.NumParticles;
for p = 1:3
    nexttile;
    xv = T6.(panelInputs{p});
    sc = scatter(xv, density, 50, cAdv, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    colormap(gca, redblueMap()); cb = colorbar; cb.Label.String = '\Delta linkAcc (6D−3D)';
    cb.Label.FontWeight = 'bold';
    caxis([-0.1 0.1]);
    if panelLogX(p), set(gca, 'XScale', 'log'); end
    xlabel(panelLabels{p}, 'FontWeight','bold','FontSize', 12);
    ylabel('density (# particles)', 'FontWeight','bold','FontSize', 12);
    set(gca, 'FontName','Calibri','FontSize', 11, 'LineWidth', 1.5, 'TickDir', 'in', 'Box', 'on');
    grid on;
    title(sprintf('density vs %s', panelLabels{p}), 'FontWeight','bold','FontSize', 13);
end
title(tl3, sprintf('Where is the 6D advantage on linking accuracy largest? (N=%d)', nP), ...
    'FontName','Calibri', 'FontWeight','bold', 'FontSize', 14);

png3 = fullfile(outDir, 'mc_pair_scatter.png');
exportgraphics(fig3, png3, 'Resolution', 600); close(fig3);
fprintf('Saved: %s\n', png3);

%% ========================================================================
%  CONSOLE + CSV SUMMARY
%  ========================================================================
fprintf('\n=== Overall 6D vs 3D delta means (across N=%d samples) ===\n', nP);
fprintf('  %-22s | %10s | %10s\n', 'metric', 'mean Δ', 'fraction 6D-better');
fprintf('  %s\n', repmat('-', 1, 55));
sumRows = {};
for k = 1:length(outputCols)
    d = Y(:, k);
    fracBetter = mean(d * outputDir(k) > 0);
    fprintf('  %-22s | %+10.4f | %10.1f%%\n', ...
        outputCols{k}, mean(d), 100*fracBetter);
    sumRows(end+1, :) = {outputCols{k}, mean(d), median(d), fracBetter, ...
        std(d) / sqrt(length(d))}; %#ok<AGROW>
end
Tsum = cell2table(sumRows, 'VariableNames', ...
    {'metric','meanDelta','medianDelta','fraction6Dbetter','SE'});
writetable(Tsum, fullfile(outDir, 'mc_summary.csv'));
fprintf('\nSaved summary CSV: %s\n', fullfile(outDir, 'mc_summary.csv'));

%% Top-3 drivers per metric (largest |signed Spearman rho|)
fprintf('\n=== Top-3 input drivers of each metric''s 6D advantage ===\n');
for k = 1:length(outputCols)
    [~, idx] = sort(abs(S_signed(:, k)), 'descend');
    fprintf('  %-22s : ', outputCols{k});
    for j = 1:min(3, length(idx))
        if isnan(S_signed(idx(j), k)), continue; end
        fprintf('%s(\\rho=%+.2f)  ', inputDisp{idx(j)}, S_signed(idx(j), k));
    end
    fprintf('\n');
end

fprintf('\nDone.\n');
end

%% Diverging red-white-blue colormap
function cm = redblueMap()
    n = 128;
    r = [linspace(0.05, 1, n)'; linspace(1, 0.85, n)'];
    g = [linspace(0.20, 1, n)'; linspace(1, 0.15, n)'];
    b = [linspace(0.55, 1, n)'; linspace(1, 0.10, n)'];
    cm = [r g b];
end
