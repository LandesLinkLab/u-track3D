function movieInfo = convert6DSMOLMtoMovieInfo(inputData, varargin)
%CONVERT6DSMOLMTOMOVIEINFO Convert 6D-SMOLM fusion output to u-track movieInfo format.
%
%   Converts the output of the mDeepSTORM3D + mDeep-SMOLM neural network
%   fusion pipeline into the movieInfo structure required by u-track3D's
%   trackCloseGapsKalmanSparse function.
%
%   The extra orientation fields (theta, phi, omega) are stored as
%   additional fields in the movieInfo struct. The core u-track3D tracker
%   ignores unknown fields; our custom cost functions read them.
%
%   SYNOPSIS:
%       movieInfo = convert6DSMOLMtoMovieInfo(inputData)
%       movieInfo = convert6DSMOLMtoMovieInfo(inputData, 'Name', Value, ...)
%
%   INPUT:
%       inputData : Can be one of:
%           (a) Path to a CSV file with columns:
%               frame, x_nm, y_nm, z_nm, theta, phi, omega, intensity
%               (Column order can be specified via 'ColumnMap'.)
%           (b) Path to a MAT file containing a table or matrix variable.
%           (c) A numeric matrix [N x 8] with the columns above.
%           (d) A MATLAB table with the corresponding column names.
%
%   NAME-VALUE PARAMETERS:
%       'ColumnMap'    : Struct mapping field names to column indices.
%                        Default (Fusion output format):
%                        struct('frame',1, 'intensity',2, 'x_nm',3,
%                               'y_nm',4, 'z_nm',5, 'theta',6, 'phi',7, 'gamma',8)
%                        Use 'gamma' key for gamma input or 'omega' for omega input.
%       'AngleUnits'   : 'degrees' (default) or 'radians'. Fusion outputs degrees.
%       'WobbleInput'  : 'gamma' (default) or 'omega'. Fusion outputs gamma [0,1].
%                        If gamma: converts via omega = pi*(3 - sqrt(1 + 8*gamma)).
%       'UncertaintyX' : Default x uncertainty in nm (if not in data). Default: 20.
%       'UncertaintyY' : Default y uncertainty in nm. Default: 20.
%       'UncertaintyZ' : Default z uncertainty in nm. Default: 50.
%       'UncertaintyTheta' : Default theta uncertainty in rad. Default: 0.1.
%       'UncertaintyPhi'   : Default phi uncertainty in rad. Default: 0.2.
%       'UncertaintyOmega' : Default omega uncertainty in sr. Default: 0.3.
%       'HasUncertainties' : Logical. If true, expects paired columns. Default: false.
%       'VariableName'     : For MAT files, name of variable to load. Default: auto-detect.
%       'Verbose'          : Logical. Print progress info. Default: true.
%
%   OUTPUT:
%       movieInfo : [numFrames x 1] struct array with fields:
%           .xCoord  - [P x 2] x-positions and uncertainties (nm)
%           .yCoord  - [P x 2] y-positions and uncertainties (nm)
%           .zCoord  - [P x 2] z-positions and uncertainties (nm)
%           .amp     - [P x 2] intensities and uncertainties
%           .theta   - [P x 2] polar angles and uncertainties (rad)
%           .phi     - [P x 2] azimuthal angles and uncertainties (rad)
%           .omega   - [P x 2] wobble solid angles and uncertainties (sr)
%
%   NOTES:
%       - Frame indices in the input should be 1-based integers.
%       - Empty frames are allowed (the struct will have empty arrays).
%       - The function validates angle ranges and warns about out-of-range values.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'inputData');
addParameter(p, 'ColumnMap', struct('frame',1, 'intensity',2, 'x_nm',3, ...
    'y_nm',4, 'z_nm',5, 'theta',6, 'phi',7, 'gamma',8), @isstruct);
addParameter(p, 'AngleUnits', 'degrees', @ischar);  % 'degrees' or 'radians'
addParameter(p, 'WobbleInput', 'gamma', @ischar);   % 'gamma' or 'omega'
addParameter(p, 'UncertaintyX', 20, @isnumeric);
addParameter(p, 'UncertaintyY', 20, @isnumeric);
addParameter(p, 'UncertaintyZ', 50, @isnumeric);
addParameter(p, 'UncertaintyTheta', 0.1, @isnumeric);  % in rad (internal)
addParameter(p, 'UncertaintyPhi', 0.2, @isnumeric);    % in rad (internal)
addParameter(p, 'UncertaintyOmega', 0.3, @isnumeric);  % in sr (internal)
addParameter(p, 'HasUncertainties', false, @islogical);
addParameter(p, 'VariableName', '', @ischar);
addParameter(p, 'Verbose', true, @islogical);
parse(p, inputData, varargin{:});

opts = p.Results;
colMap = opts.ColumnMap;

%% --- Load data ---
if ischar(inputData) || isstring(inputData)
    inputData = char(inputData);
    [~, ~, ext] = fileparts(inputData);
    
    switch lower(ext)
        case '.csv'
            if opts.Verbose
                fprintf('Loading CSV file: %s\n', inputData);
            end
            data = readmatrix(inputData);
            
        case '.mat'
            if opts.Verbose
                fprintf('Loading MAT file: %s\n', inputData);
            end
            S = load(inputData);
            fnames = fieldnames(S);
            if ~isempty(opts.VariableName)
                data = S.(opts.VariableName);
            elseif length(fnames) == 1
                data = S.(fnames{1});
            else
                % Try to find a table or large matrix
                found = false;
                for i = 1:length(fnames)
                    if istable(S.(fnames{i})) || (isnumeric(S.(fnames{i})) && size(S.(fnames{i}),2) >= 6)
                        data = S.(fnames{i});
                        found = true;
                        if opts.Verbose
                            fprintf('  Auto-selected variable: %s\n', fnames{i});
                        end
                        break;
                    end
                end
                if ~found
                    error('convert6DSMOLMtoMovieInfo:ambiguousMAT', ...
                        'MAT file has multiple variables. Specify ''VariableName''.');
                end
            end
            
        otherwise
            error('convert6DSMOLMtoMovieInfo:unsupportedFormat', ...
                'Unsupported file format: %s. Use .csv or .mat.', ext);
    end
elseif istable(inputData)
    data = inputData;
elseif isnumeric(inputData)
    data = inputData;
else
    error('convert6DSMOLMtoMovieInfo:invalidInput', ...
        'Input must be a file path, table, or numeric matrix.');
end

%% --- Convert table to matrix if needed ---
if istable(data)
    % Try to match column names to expected fields
    colNames = data.Properties.VariableNames;
    expectedFields = {'frame','x_nm','y_nm','z_nm','theta','phi','omega','intensity'};
    
    % Check if table has named columns matching our expectations
    hasNamedCols = all(ismember(expectedFields, colNames));
    
    if hasNamedCols
        dataMatrix = [data.frame, data.x_nm, data.y_nm, data.z_nm, ...
            data.theta, data.phi, data.omega, data.intensity];
        % Reset column map to standard order since we extracted in order
        colMap = struct('frame',1, 'x_nm',2, 'y_nm',3, 'z_nm',4, ...
            'theta',5, 'phi',6, 'omega',7, 'intensity',8);
    else
        dataMatrix = table2array(data);
    end
    data = dataMatrix;
end

%% --- Extract columns ---
frames    = data(:, colMap.frame);
x_nm      = data(:, colMap.x_nm);
y_nm      = data(:, colMap.y_nm);
z_nm      = data(:, colMap.z_nm);
intensity = data(:, colMap.intensity);

% Extract angles — handle field name 'gamma' or 'omega' in column map
theta_raw = data(:, colMap.theta);
phi_raw   = data(:, colMap.phi);

if isfield(colMap, 'gamma')
    gamma_raw = data(:, colMap.gamma);
    hasGamma = true;
elseif isfield(colMap, 'omega')
    omega_raw = data(:, colMap.omega);
    hasGamma = false;
else
    error('convert6DSMOLMtoMovieInfo:noWobble', ...
        'ColumnMap must have either ''gamma'' or ''omega'' field.');
end

%% --- Unit conversions ---

% Convert angles to radians (internal convention)
if strcmpi(opts.AngleUnits, 'degrees')
    theta = theta_raw * (pi / 180);  % degrees -> radians
    phi   = phi_raw * (pi / 180);
    if opts.Verbose
        fprintf('  Converting angles: degrees -> radians\n');
    end
elseif strcmpi(opts.AngleUnits, 'radians')
    theta = theta_raw;
    phi   = phi_raw;
else
    error('convert6DSMOLMtoMovieInfo:badUnits', ...
        'AngleUnits must be ''degrees'' or ''radians''. Got: %s', opts.AngleUnits);
end

% Convert wobble: gamma -> omega (solid angle in steradians)
% gamma = 1 - 3*omega/(4*pi) + omega^2/(8*pi^2)
% Inverse: omega = pi * (3 - sqrt(1 + 8*gamma))    [smaller root]
if hasGamma || strcmpi(opts.WobbleInput, 'gamma')
    if ~exist('gamma_raw', 'var')
        gamma_raw = omega_raw;  % In case user set WobbleInput override
    end
    gamma_raw = max(0, min(1, gamma_raw));  % Clamp to [0, 1]
    omega = pi * (3 - sqrt(1 + 8 * gamma_raw));
    omega = real(omega);  % Safety: avoid complex from numerical edge cases
    omega = max(0, omega);
    if opts.Verbose
        fprintf('  Converting wobble: gamma [0,1] -> omega [0, 2pi] sr\n');
        fprintf('    gamma range: [%.3f, %.3f] -> omega range: [%.3f, %.3f] sr\n', ...
            min(gamma_raw), max(gamma_raw), min(omega), max(omega));
    end
else
    omega = omega_raw;
end

%% --- Validate data ---
numLoc = size(data, 1);

% Ensure frames are integers
if any(frames ~= round(frames))
    warning('convert6DSMOLMtoMovieInfo:nonIntegerFrames', ...
        'Frame indices are not integers. Rounding.');
    frames = round(frames);
end

% Check for frame 0 (convert to 1-based if needed)
if any(frames == 0)
    warning('convert6DSMOLMtoMovieInfo:zeroBased', ...
        'Frame indices start at 0. Converting to 1-based indexing.');
    frames = frames + 1;
end

minFrame = min(frames);
maxFrame = max(frames);

if minFrame < 1
    error('convert6DSMOLMtoMovieInfo:negativeFrames', ...
        'Frame indices contain negative values.');
end

% Validate angle ranges
nBadTheta = sum(theta < 0 | theta > pi/2);
if nBadTheta > 0
    warning('convert6DSMOLMtoMovieInfo:thetaRange', ...
        '%d localizations have theta outside [0, pi/2]. Clamping.', nBadTheta);
    theta = max(0, min(pi/2, theta));
end

nBadPhi = sum(phi < -pi | phi > pi);
if nBadPhi > 0
    warning('convert6DSMOLMtoMovieInfo:phiRange', ...
        '%d localizations have phi outside [-pi, pi]. Wrapping.', nBadPhi);
    phi = mod(phi + pi, 2*pi) - pi;  % Wrap to [-pi, pi]
end

nBadOmega = sum(omega < 0 | omega > 2*pi);
if nBadOmega > 0
    warning('convert6DSMOLMtoMovieInfo:omegaRange', ...
        '%d localizations have omega outside [0, 2*pi]. Clamping.', nBadOmega);
    omega = max(0, min(2*pi, omega));
end

%% --- Build movieInfo struct array ---
numFrames = maxFrame;  % Include all frames from 1 to maxFrame

% Pre-allocate struct array with empty fields
movieInfo = repmat(struct('xCoord',[], 'yCoord',[], 'zCoord',[], ...
    'amp',[], 'theta',[], 'phi',[], 'omega',[], 'nnDist',[]), numFrames, 1);

% Fill each frame
for f = 1:numFrames
    idx = (frames == f);
    nPart = sum(idx);
    
    if nPart == 0
        % Empty frame: store empty [0 x 2] arrays
        movieInfo(f).xCoord = zeros(0, 2);
        movieInfo(f).yCoord = zeros(0, 2);
        movieInfo(f).zCoord = zeros(0, 2);
        movieInfo(f).amp    = zeros(0, 2);
        movieInfo(f).theta  = zeros(0, 2);
        movieInfo(f).phi    = zeros(0, 2);
        movieInfo(f).omega  = zeros(0, 2);
        movieInfo(f).nnDist = zeros(0, 1);
    else
        % Position: [value, uncertainty]
        movieInfo(f).xCoord = [x_nm(idx), repmat(opts.UncertaintyX, nPart, 1)];
        movieInfo(f).yCoord = [y_nm(idx), repmat(opts.UncertaintyY, nPart, 1)];
        movieInfo(f).zCoord = [z_nm(idx), repmat(opts.UncertaintyZ, nPart, 1)];
        
        % Amplitude
        movieInfo(f).amp = [intensity(idx), zeros(nPart, 1)];
        
        % Orientation: [value, uncertainty]
        movieInfo(f).theta = [theta(idx), repmat(opts.UncertaintyTheta, nPart, 1)];
        movieInfo(f).phi   = [phi(idx), repmat(opts.UncertaintyPhi, nPart, 1)];
        movieInfo(f).omega = [omega(idx), repmat(opts.UncertaintyOmega, nPart, 1)];
        
        % Nearest-neighbor distances (required by linkFeaturesKalmanSparse)
        if nPart == 1
            movieInfo(f).nnDist = Inf;
        else
            coords = [movieInfo(f).xCoord(:,1), movieInfo(f).yCoord(:,1), movieInfo(f).zCoord(:,1)];
            D = pdist2(coords, coords);
            D(D == 0) = Inf;  % Ignore self-distance
            movieInfo(f).nnDist = min(D, [], 2);  % [nPart x 1]
        end
    end
end

%% --- Summary ---
if opts.Verbose
    nFramesWithParticles = sum(arrayfun(@(s) size(s.xCoord,1), movieInfo) > 0);
    avgParticles = numLoc / numFrames;
    fprintf('--- convert6DSMOLMtoMovieInfo summary ---\n');
    fprintf('  Total localizations: %d\n', numLoc);
    fprintf('  Total frames:        %d\n', numFrames);
    fprintf('  Frames with data:    %d (%.1f%%)\n', nFramesWithParticles, ...
        100*nFramesWithParticles/numFrames);
    fprintf('  Avg particles/frame: %.1f\n', avgParticles);
    fprintf('  movieInfo size:      [%d x 1] struct\n', numFrames);
    fprintf('----------------------------------------\n');
end

end
