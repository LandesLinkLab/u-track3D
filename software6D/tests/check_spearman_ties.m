function check_spearman_ties()
%CHECK_SPEARMAN_TIES  Does our MC data have ties that would invalidate
%the simplified Spearman formula rho = 1 - 6*sum(d^2)/(n(n^2-1))?
%
%  Compares the simplified formula against MATLAB's tie-aware Spearman
%  on (input, sign-aligned-delta) pairs. If they agree to ~1e-12, the
%  simplified form is OK; if they disagree, the manuscript should cite
%  a tie-aware reference.

thisDir = fileparts(mfilename('fullpath'));
csvPath = fullfile(thisDir, 'benchmark_output', 'mc', 'mc_results.csv');
T = readtable(csvPath);
T3 = T(strcmp(T.mode, '3D'), :);
T6 = T(strcmp(T.mode, '6D'), :);
[~, i3, i6] = intersect(T3.sample, T6.sample);
T3 = T3(i3, :); T6 = T6(i6, :);

inputs = {'NumParticles', 'Dtrans', 'Drot', 'GammaMean', 'BlinkOffRate', ...
          'BlinkOnRate', 'Signal', 'NoiseXY', 'NoiseTheta', ...
          'wOrientSampled', 'maxSearchRadius', 'timeWindow', 'orientCRLB'};
metrics = {'linkAcc', 'meanCompleteness', 'meanPurity', 'nMerged', ...
           'alpha', 'JSC', 'oriMedianDeg'};

fprintf('%-18s | %-22s | %12s | %12s | %12s\n', ...
    'Input', 'Metric', 'rho_simple', 'rho_matlab', 'diff');
fprintf(repmat('-', 1, 90)); fprintf('\n');

maxDiff = 0;
nTies   = 0;
for ii = 1:length(inputs)
    x = T6.(inputs{ii});
    for jj = 1:length(metrics)
        y = T6.(metrics{jj}) - T3.(metrics{jj});
        % Tie counts in raw values
        nTiesX = length(x) - length(unique(x));
        nTiesY = length(y) - length(unique(y));
        if nTiesX > 0 || nTiesY > 0, nTies = nTies + 1; end

        % Simplified Spearman (no-ties assumption)
        Rx = tiedrank(x);
        Ry = tiedrank(y);
        n  = length(x);
        d  = Rx - Ry;
        rhoSimple = 1 - 6 * sum(d.^2) / (n * (n^2 - 1));

        % MATLAB built-in (tie-aware Pearson on tied ranks)
        rhoMatlab = corr(x, y, 'Type', 'Spearman');

        diff = abs(rhoSimple - rhoMatlab);
        if diff > maxDiff, maxDiff = diff; end
        if diff > 1e-6
            fprintf('%-18s | %-22s | %+.6f | %+.6f | %.2e\n', ...
                inputs{ii}, metrics{jj}, rhoSimple, rhoMatlab, diff);
        end
    end
end
fprintf('\nMax |rho_simple - rho_matlab| across all pairings: %.4e\n', maxDiff);
fprintf('Pairings with ties in raw values: %d / %d\n', ...
    nTies, length(inputs) * length(metrics));
if maxDiff < 1e-6
    fprintf('Simplified Spearman formula is exact for this data.\n');
else
    fprintf('Simplified formula differs - cite tie-aware reference.\n');
end
end
