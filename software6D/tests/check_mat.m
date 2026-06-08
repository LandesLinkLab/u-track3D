fprintf('=== sweeps_results.mat contents ===\n');
S = load(fullfile(fileparts(mfilename('fullpath')), 'benchmark_output', 'sweeps_v3', 'sweeps_results.mat'));
disp(fieldnames(S));
if isfield(S,'T')
    fprintf('T table columns:\n');
    disp(S.T.Properties.VariableNames);
    fprintf('T rows = %d\n', height(S.T));
end

fprintf('\n=== density_blink.mat contents ===\n');
S = load(fullfile(fileparts(mfilename('fullpath')), 'benchmark_output', 'sweeps_v3', 'density_blink.mat'));
disp(fieldnames(S));
if isfield(S,'tracksGT_3D') || isfield(S,'tracksGT'), fprintf('HAS raw tracks!\n'); end

fprintf('\n=== mc_results.mat contents ===\n');
S = load(fullfile(fileparts(mfilename('fullpath')), 'benchmark_output', 'mc', 'mc_results.mat'));
disp(fieldnames(S));
