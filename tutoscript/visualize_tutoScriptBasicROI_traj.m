%% Load tracking results
trackFile = 'C:\Users\emilg2\AppData\Local\Temp\testBasicROI\endocyticPits_cropped_analysis\tracks\Channel_1_tracking_result.mat';
load(trackFile);

%% Pixel sizes in μm
pixelXY = 0.104;  % μm
pixelZ  = 0.350;  % μm

%% Visualize tracks with time-coded colors
figure('Position', [100 100 1200 500]);

nTracks = length(tracksFinal);

% Get total number of frames for colormap normalization
maxFrame = 0;
for i = 1:nTracks
    nFrames = size(tracksFinal(i).tracksFeatIndxCG, 2);
    startFrame = tracksFinal(i).seqOfEvents(1,1);
    endFrame = startFrame + nFrames - 1;
    if endFrame > maxFrame
        maxFrame = endFrame;
    end
end

% Colormap
cmap = parula(maxFrame);

% 3D view
subplot(1,2,1);
hold on;

for i = 1:min(nTracks, 100)
    coords = tracksFinal(i).tracksCoordAmpCG;
    x = coords(1:8:end) * pixelXY;  % μm
    y = coords(2:8:end) * pixelXY;  % μm
    z = coords(3:8:end) * pixelZ;   % μm
    
    startFrame = tracksFinal(i).seqOfEvents(1,1);
    nPts = length(x);
    
    for j = 1:(nPts-1)
        if ~isnan(x(j)) && ~isnan(x(j+1))
            frameIdx = startFrame + j - 1;
            plot3(x(j:j+1), y(j:j+1), z(j:j+1), '-', ...
                  'Color', cmap(frameIdx,:), 'LineWidth', 1.5);
        end
    end
end

xlabel('X (μm)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y (μm)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Z (μm)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('3D Tracks (%d total)', nTracks), 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
grid on; view(3);
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Frame';
cb.Label.FontSize = 14;
cb.Label.FontWeight = 'bold';
cb.FontSize = 12;
cb.FontWeight = 'bold';
clim([1 maxFrame]);

% 2D XY view
subplot(1,2,2);
hold on;

for i = 1:min(nTracks, 100)
    coords = tracksFinal(i).tracksCoordAmpCG;
    x = coords(1:8:end) * pixelXY;
    y = coords(2:8:end) * pixelXY;
    
    startFrame = tracksFinal(i).seqOfEvents(1,1);
    nPts = length(x);
    
    for j = 1:(nPts-1)
        if ~isnan(x(j)) && ~isnan(x(j+1))
            frameIdx = startFrame + j - 1;
            plot(x(j:j+1), y(j:j+1), '-', ...
                 'Color', cmap(frameIdx,:), 'LineWidth', 1.5);
        end
    end
end

xlabel('X (μm)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y (μm)', 'FontSize', 14, 'FontWeight', 'bold');
title('2D Projection (XY)', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
grid on; axis equal;
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Frame';
cb.Label.FontSize = 14;
cb.Label.FontWeight = 'bold';
cb.FontSize = 12;
cb.FontWeight = 'bold';
clim([1 maxFrame]);

fprintf('\nTotal tracks: %d\n', nTracks);
fprintf('Total frames: %d\n', maxFrame);
fprintf('Pixel size: XY = %.3f μm, Z = %.3f μm\n', pixelXY, pixelZ);

%% Export figure at 600 DPI
% print(gcf, 'tracks_visualization.png', '-dpng', '-r600');