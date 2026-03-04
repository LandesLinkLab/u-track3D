function trajData = extractTracks6D(tracksFinal, movieInfo, varargin)
%EXTRACTTRACKS6D Convert u-track tracksFinal output to 6D trajectory matrix.
%
%   Converts the tracksFinal struct array from trackCloseGapsKalmanSparse
%   into a flat trajectory matrix with orientation data reattached from
%   the original movieInfo.
%
%   SYNOPSIS:
%       trajData = extractTracks6D(tracksFinal, movieInfo)
%       trajData = extractTracks6D(tracksFinal, movieInfo, 'Name', Value, ...)
%
%   INPUT:
%       tracksFinal : Output from trackCloseGapsKalmanSparse. Struct array
%                     with fields: .tracksCoordAmpCG, .tracksFeatIndxCG,
%                     .seqOfEvents
%       movieInfo   : Original movieInfo struct with orientation fields.
%
%   NAME-VALUE PARAMETERS:
%       'MinTrackLen' : Minimum track length to include. Default: 1.
%       'Verbose'     : Print stats. Default: true.
%
%   OUTPUT:
%       trajData : [N x 9] matrix with columns:
%                  [frame, trackID, x_nm, y_nm, z_nm, theta, phi, omega, intensity]
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addRequired(p, 'tracksFinal');
addRequired(p, 'movieInfo');
addParameter(p, 'MinTrackLen', 1, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
parse(p, tracksFinal, movieInfo, varargin{:});

minLen = p.Results.MinTrackLen;
verbose = p.Results.Verbose;

%% --- Extract tracks ---
nTracks = length(tracksFinal);
allRows = cell(nTracks, 1);
nValid = 0;

for iTrack = 1:nTracks
    % Get track coordinates: columns cycle as [x,y,z,amp,dx,dy,dz,damp] per frame
    coordAmp = tracksFinal(iTrack).tracksCoordAmpCG;
    featIndx = tracksFinal(iTrack).tracksFeatIndxCG;
    seqEvents = tracksFinal(iTrack).seqOfEvents;
    
    % For compound tracks (merges/splits), only take the first sub-track
    % (row 1 of the matrices)
    if size(coordAmp, 1) > 1
        coordAmp = coordAmp(1, :);
        featIndx = featIndx(1, :);
    end
    
    % Determine frame range from seqOfEvents
    startFrame = seqEvents(1, 1);
    endFrame = seqEvents(end, 1);
    nFrames = endFrame - startFrame + 1;
    
    % Parse coordinates (8 columns per frame)
    nCol = 8;
    frames = (startFrame:endFrame)';
    
    x = coordAmp(1:nCol:end)';
    y = coordAmp(2:nCol:end)';
    z = coordAmp(3:nCol:end)';
    amp = coordAmp(4:nCol:end)';
    
    % Truncate to actual track length
    nActual = min(nFrames, length(x));
    frames = frames(1:nActual);
    x = x(1:nActual);
    y = y(1:nActual);
    z = z(1:nActual);
    amp = amp(1:nActual);
    featIndx = featIndx(1:nActual);
    
    % Remove NaN entries (gaps that were not closed)
    valid = ~isnan(x);
    
    if sum(valid) < minLen
        continue;
    end
    
    frames = frames(valid);
    x = x(valid);
    y = y(valid);
    z = z(valid);
    amp = amp(valid);
    featIdx = featIndx(valid);
    nPts = length(frames);
    
    % Retrieve orientation from movieInfo
    theta = zeros(nPts, 1);
    phi   = zeros(nPts, 1);
    omega = zeros(nPts, 1);
    
    for iPt = 1:nPts
        f = frames(iPt);
        idx = featIdx(iPt);
        
        if f > 0 && f <= length(movieInfo) && idx > 0
            if isfield(movieInfo(f), 'theta') && idx <= size(movieInfo(f).theta, 1)
                theta(iPt) = movieInfo(f).theta(idx, 1);
                phi(iPt)   = movieInfo(f).phi(idx, 1);
                omega(iPt) = movieInfo(f).omega(idx, 1);
            end
        end
    end
    
    % Build trajectory rows
    nValid = nValid + 1;
    trackIDs = repmat(nValid, nPts, 1);
    allRows{nValid} = [frames, trackIDs, x, y, z, theta, phi, omega, amp];
end

%% --- Concatenate ---
allRows = allRows(1:nValid);

if nValid > 0
    trajData = vertcat(allRows{:});
else
    trajData = zeros(0, 9);
end

%% --- Summary ---
if verbose
    fprintf('--- extractTracks6D summary ---\n');
    fprintf('  Input tracks:  %d\n', nTracks);
    fprintf('  Output tracks: %d (>= %d frames)\n', nValid, minLen);
    if nValid > 0
        trackLens = accumarray(trajData(:,2), 1);
        fprintf('  Total points:  %d\n', size(trajData, 1));
        fprintf('  Mean length:   %.1f frames\n', mean(trackLens));
        fprintf('  Median length: %.1f frames\n', median(trackLens));
    end
    fprintf('-------------------------------\n');
end

end
