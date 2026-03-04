function [combined, stats] = combineTrajectories6D(trajFiles, varargin)
%COMBINETRAJECTORIES6D Combine 6D trajectory files from multiple movies.
%
%   Extended from trjR_combined_nastya.m (Landes Research Group) to handle
%   6D orientation data. Combines trajectory matrices from multiple movies
%   or crop regions, renumbering track IDs to prevent collisions.
%
%   SYNOPSIS:
%       [combined, stats] = combineTrajectories6D(trajFiles)
%       [combined, stats] = combineTrajectories6D(trajFiles, 'Name', Value, ...)
%
%   INPUT:
%       trajFiles : Cell array of file paths (.csv or .mat), OR
%                   cell array of [N x 9] trajectory matrices.
%                   Each with columns:
%                   [frame, trackID, x_nm, y_nm, z_nm, theta, phi, omega, intensity]
%
%   NAME-VALUE PARAMETERS:
%       'AddMovieID'   : Logical. Add a 10th column with movie index. Default: true.
%       'OutputFile'   : Path to save combined output. Default: '' (no save).
%       'VariableName' : Variable name to look for in MAT files. Default: auto-detect.
%       'Verbose'      : Print stats. Default: true.
%
%   OUTPUT:
%       combined : [M x 9] or [M x 10] combined trajectory matrix.
%                  Track IDs renumbered sequentially across all movies.
%                  If AddMovieID = true, column 10 = movie index.
%       stats    : Struct with combination statistics.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'trajFiles', @iscell);
addParameter(p, 'AddMovieID', true, @islogical);
addParameter(p, 'OutputFile', '', @ischar);
addParameter(p, 'VariableName', '', @ischar);
addParameter(p, 'Verbose', true, @islogical);
parse(p, trajFiles, varargin{:});

opts = p.Results;

%% --- Load and combine ---
nMovies = length(trajFiles);
allTrajs = cell(nMovies, 1);
nTracks = zeros(nMovies, 1);
nLocs = zeros(nMovies, 1);

maxTrackID = 0;  % Running maximum for renumbering

for iMov = 1:nMovies
    % Load data
    if ischar(trajFiles{iMov}) || isstring(trajFiles{iMov})
        fpath = char(trajFiles{iMov});
        [~, ~, ext] = fileparts(fpath);
        
        switch lower(ext)
            case '.csv'
                data = readmatrix(fpath);
            case '.mat'
                S = load(fpath);
                fnames = fieldnames(S);
                if ~isempty(opts.VariableName)
                    data = S.(opts.VariableName);
                elseif length(fnames) == 1
                    data = S.(fnames{1});
                else
                    % Try to find large matrix
                    for k = 1:length(fnames)
                        if isnumeric(S.(fnames{k})) && size(S.(fnames{k}),2) >= 9
                            data = S.(fnames{k});
                            break;
                        end
                    end
                end
            otherwise
                error('combineTrajectories6D:unsupported', ...
                    'Unsupported file format: %s', ext);
        end
        
        if istable(data)
            data = table2array(data);
        end
    elseif isnumeric(trajFiles{iMov})
        data = trajFiles{iMov};
    else
        error('combineTrajectories6D:invalidInput', ...
            'Entry %d must be a file path or numeric matrix.', iMov);
    end
    
    if isempty(data)
        allTrajs{iMov} = zeros(0, 9);
        continue;
    end
    
    % Renumber track IDs to be globally unique
    uniqueIDs = unique(data(:, 2));
    nTracks(iMov) = length(uniqueIDs);
    nLocs(iMov) = size(data, 1);
    
    % Create mapping: old ID -> new sequential ID
    for iID = 1:length(uniqueIDs)
        data(data(:, 2) == uniqueIDs(iID), 2) = maxTrackID + iID;
    end
    maxTrackID = maxTrackID + length(uniqueIDs);
    
    % Add movie ID column if requested
    if opts.AddMovieID
        data = [data, repmat(iMov, size(data, 1), 1)]; %#ok<AGROW>
    end
    
    allTrajs{iMov} = data;
end

%% --- Concatenate ---
combined = vertcat(allTrajs{:});

%% --- Statistics ---
stats = struct();
stats.nMovies     = nMovies;
stats.nTotalTracks = sum(nTracks);
stats.nTotalLocs  = sum(nLocs);
stats.tracksPerMovie = nTracks;
stats.locsPerMovie = nLocs;

if opts.Verbose
    fprintf('--- combineTrajectories6D summary ---\n');
    fprintf('  Movies combined:     %d\n', nMovies);
    fprintf('  Total tracks:        %d\n', sum(nTracks));
    fprintf('  Total localizations: %d\n', sum(nLocs));
    for iMov = 1:nMovies
        fprintf('  Movie %d: %d tracks, %d localizations\n', ...
            iMov, nTracks(iMov), nLocs(iMov));
    end
    fprintf('-------------------------------------\n');
end

%% --- Save if requested ---
if ~isempty(opts.OutputFile)
    [~, ~, ext] = fileparts(opts.OutputFile);
    switch lower(ext)
        case '.csv'
            if opts.AddMovieID
                header = 'frame,trackID,x_nm,y_nm,z_nm,theta,phi,omega,intensity,movieID';
            else
                header = 'frame,trackID,x_nm,y_nm,z_nm,theta,phi,omega,intensity';
            end
            fid = fopen(opts.OutputFile, 'w');
            fprintf(fid, '%s\n', header);
            fclose(fid);
            dlmwrite(opts.OutputFile, combined, '-append', 'delimiter', ',', 'precision', '%.6f');
            
        case '.mat'
            trajCombined = combined; %#ok<NASGU>
            save(opts.OutputFile, 'trajCombined');
            
        otherwise
            writematrix(combined, opts.OutputFile);
    end
    
    if opts.Verbose
        fprintf('  Saved to: %s\n', opts.OutputFile);
    end
end

end
