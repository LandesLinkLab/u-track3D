% Debug script to understand MLE issue
% Run this from the same folder where you ran run_demo.m
% (typically the tests/ folder)

% Find the output directory
if exist('output_demo', 'dir')
    outputDir = 'output_demo';
elseif exist(fullfile('..', 'output_demo'), 'dir')
    outputDir = fullfile('..', 'output_demo');
elseif exist(fullfile(fileparts(mfilename('fullpath')), 'output_demo'), 'dir')
    outputDir = fullfile(fileparts(mfilename('fullpath')), 'output_demo');
else
    % Ask user
    outputDir = uigetdir(pwd, 'Select output_demo folder');
    if outputDir == 0
        error('Please run run_demo.m first, or select the output_demo folder');
    end
end

fprintf('Using output directory: %s\n', outputDir);

% Load data
trajFile = fullfile(outputDir, 'trajectories_6D_filtered.mat');
simFile = fullfile(outputDir, 'sim_brownian.mat');

if ~exist(trajFile, 'file')
    error('Cannot find %s\nPlease run run_demo.m first.', trajFile);
end
if ~exist(simFile, 'file')
    error('Cannot find %s\nPlease run run_demo.m first.', simFile);
end

load(trajFile, 'trajFiltered');
load(simFile, 'groundTruth', 'est');

dt = 0.05;  % Frame time
trueDtrans = groundTruth.Dtrans(1);  % True D
trueDrot = groundTruth.Drot(1);

fprintf('\n========================================\n');
fprintf('=== DEBUG MLE DIFFUSION ESTIMATION ===\n');
fprintf('========================================\n');

fprintf('\nSimulation parameters:\n');
fprintf('  True D_trans = %.0f nm^2/s\n', trueDtrans);
fprintf('  True D_rot = %.4f rad^2/s\n', trueDrot);
fprintf('  Frame time dt = %.3f s\n', dt);

% Expected MSD
fprintf('\nExpected MSD per step (no gaps):\n');
fprintf('  Spatial: 6*D_t*dt = 6*%.0f*%.3f = %.1f nm^2\n', trueDtrans, dt, 6*trueDtrans*dt);
fprintf('  Angular: 2*D_r*dt = 2*%.4f*%.3f = %.6f rad^2\n', trueDrot, dt, 2*trueDrot*dt);

% Check raw simulation output directly from groundTruth
fprintf('\n--- Checking groundTruth.trajectories{1} directly ---\n');
gt1 = groundTruth.trajectories{1};
fprintf('  Size: %d x %d\n', size(gt1,1), size(gt1,2));
fprintf('  Columns: [frame, signal, x, y, z, theta_deg, phi_deg, gamma, trackID, isInterp]\n');

% Extract positions from ground truth
frames_gt = gt1(:,1);
x_gt = gt1(:,3);
y_gt = gt1(:,4);
z_gt = gt1(:,5);
theta_gt = gt1(:,6) * pi/180;  % Convert to rad
phi_gt = gt1(:,7) * pi/180;

% Compute spatial displacements
dx_gt = diff(x_gt);
dy_gt = diff(y_gt);
dz_gt = diff(z_gt);
dr2_gt = dx_gt.^2 + dy_gt.^2 + dz_gt.^2;
dFrames_gt = diff(frames_gt);

fprintf('\nGround truth statistics (particle 1):\n');
fprintf('  Number of frames: %d\n', size(gt1,1));
fprintf('  Frame range: %d to %d\n', min(frames_gt), max(frames_gt));
fprintf('  Mean frame gap: %.2f (should be 1 with no blinking)\n', mean(dFrames_gt));
fprintf('  Num gaps > 1: %d\n', sum(dFrames_gt > 1));

fprintf('\nSpatial displacements:\n');
fprintf('  Mean dr^2 = %.2f nm^2 (expected: %.2f)\n', mean(dr2_gt), 6*trueDtrans*dt);
fprintf('  Std dr^2 = %.2f nm^2\n', std(dr2_gt));
fprintf('  Mean dr^2/dFrames = %.2f nm^2\n', mean(dr2_gt ./ dFrames_gt));

% MLE(1) estimate
D_est_gt = mean(dr2_gt ./ dFrames_gt) / (6 * dt);
fprintf('\nMLE estimate from groundTruth:\n');
fprintf('  D_trans = %.0f nm^2/s (true: %.0f)\n', D_est_gt, trueDtrans);
fprintf('  Ratio: %.2f\n', D_est_gt / trueDtrans);

% Angular displacements
dAlpha_gt = zeros(length(theta_gt)-1, 1);
for i = 1:length(dAlpha_gt)
    % Compute angular distance using dipole symmetry
    n1 = [sin(theta_gt(i))*cos(phi_gt(i)), sin(theta_gt(i))*sin(phi_gt(i)), cos(theta_gt(i))];
    n2 = [sin(theta_gt(i+1))*cos(phi_gt(i+1)), sin(theta_gt(i+1))*sin(phi_gt(i+1)), cos(theta_gt(i+1))];
    dAlpha_gt(i) = acos(min(abs(dot(n1, n2)), 1));
end
dAlpha2_gt = dAlpha_gt.^2;

fprintf('\nAngular displacements:\n');
fprintf('  Mean dAlpha^2 = %.6f rad^2 (expected: %.6f)\n', mean(dAlpha2_gt), 2*trueDrot*dt);
fprintf('  Std dAlpha^2 = %.6f rad^2\n', std(dAlpha2_gt));

D_rot_est_gt = mean(dAlpha2_gt ./ dFrames_gt) / (2 * dt);
fprintf('\nMLE estimate for D_rot:\n');
fprintf('  D_rot = %.4f rad^2/s (true: %.4f)\n', D_rot_est_gt, trueDrot);
fprintf('  Ratio: %.2f\n', D_rot_est_gt / trueDrot);

fprintf('\n--- Checking trajFiltered ---\n');
trackIDs = unique(trajFiltered(:,2));
fprintf('  Number of tracks: %d\n', length(trackIDs));

traj = trajFiltered(trajFiltered(:,2) == trackIDs(1), :);
traj = sortrows(traj, 1);
fprintf('  Track 1 length: %d frames\n', size(traj,1));

frames = traj(:,1);
x = traj(:,3);
y = traj(:,4);
z = traj(:,5);

dx = diff(x);
dy = diff(y);
dz = diff(z);
dr2 = dx.^2 + dy.^2 + dz.^2;
dFrames = diff(frames);

fprintf('\nTrajFiltered statistics (track 1):\n');
fprintf('  Mean frame gap: %.2f\n', mean(dFrames));
fprintf('  Mean dr^2 = %.2f nm^2\n', mean(dr2));
fprintf('  MLE D_trans = %.0f nm^2/s\n', mean(dr2 ./ dFrames) / (6 * dt));

fprintf('\nFirst 10 displacements:\n');
fprintf('  Step   dFrames   dx        dy        dz        dr^2\n');
for i = 1:min(10, length(dr2))
    fprintf('  %3d    %3d     %8.2f  %8.2f  %8.2f  %8.2f\n', ...
        i, dFrames(i), dx(i), dy(i), dz(i), dr2(i));
end

% === NEW: Compare actual positions ===
fprintf('\n--- COMPARING POSITIONS: groundTruth vs trajFiltered ---\n');
fprintf('Ground truth particle 1, first 5 positions:\n');
fprintf('  Frame    x_gt       y_gt       z_gt\n');
for i = 1:5
    fprintf('  %3d    %8.2f   %8.2f   %8.2f\n', ...
        gt1(i,1), gt1(i,3), gt1(i,4), gt1(i,5));
end

fprintf('\nTrajFiltered track 1, first 5 positions:\n');
fprintf('  Frame    x_filt     y_filt     z_filt\n');
for i = 1:5
    fprintf('  %3d    %8.2f   %8.2f   %8.2f\n', ...
        traj(i,1), traj(i,3), traj(i,4), traj(i,5));
end

% Check if they match the same frames
fprintf('\nFrame alignment check:\n');
fprintf('  GT frames 1-5: %s\n', mat2str(gt1(1:5,1)'));
fprintf('  Filt frames 1-5: %s\n', mat2str(traj(1:5,1)'));

% Also check the raw observed data (est) from simulation
fprintf('\n--- Checking raw simulation output (est) ---\n');
fprintf('  est size: %d x %d\n', size(est,1), size(est,2));
fprintf('  Columns: [frame, signal, x, y, z, theta_deg, phi_deg, gamma]\n');
% Get frame 1 detections
frame1_est = est(est(:,1) == 1, :);
fprintf('  Frame 1 detections: %d\n', size(frame1_est, 1));
fprintf('  Frame 1 positions:\n');
for i = 1:size(frame1_est, 1)
    fprintf('    x=%.2f, y=%.2f, z=%.2f\n', frame1_est(i,3), frame1_est(i,4), frame1_est(i,5));
end

fprintf('\n========================================\n');
