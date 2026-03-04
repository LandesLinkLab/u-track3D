function visualizeKalmanScaling(varargin)
%VISUALIZEKALMANSCALING Show how pseudo-Kalman parameters affect cost tolerance.
%
%   Visualizes how the orientation cost function changes with:
%   - Time gap (dt)
%   - Rotational diffusion coefficient (D_rot)
%   - Orientation measurement uncertainty (sigma_orient / CRLB)
%
%   This helps understand why larger gaps or higher uncertainty leads to
%   more "tolerance" for angular mismatches during tracking.
%
%   USAGE:
%       visualizeKalmanScaling()
%       visualizeKalmanScaling('Drot', 0.1, 'FrameTime', 0.05)
%
%   PARAMETERS:
%       'Drot'        : Rotational diffusion coefficient [rad²/s]. Default: 0.05
%       'FrameTime'   : Frame time [s]. Default: 0.05
%       'OrientCRLB'  : Orientation precision [rad]. Default: 0.087 (5°)
%
%   Emil Gillett, Landes Research Group, 2026

%% Parse inputs
p = inputParser;
addParameter(p, 'Drot', 0.05, @isnumeric);
addParameter(p, 'FrameTime', 0.05, @isnumeric);
addParameter(p, 'OrientCRLB', 5 * pi/180, @isnumeric);  % 5 degrees default
parse(p, varargin{:});

Drot = p.Results.Drot;
frameTime = p.Results.FrameTime;
orientCRLB = p.Results.OrientCRLB;

%% Create figure
figure('Position', [100, 100, 1400, 500], 'Color', 'w', ...
    'Name', 'Pseudo-Kalman Orientation Cost Scaling');

%% Panel 1: Cost vs Angular Distance for different time gaps
subplot(1, 3, 1);
hold on;

dAngle_deg = linspace(0, 90, 100);  % Angular distance in degrees
dAngle_rad = dAngle_deg * pi/180;

timeGaps = [1, 2, 5, 10];  % Frame gaps
colors = lines(length(timeGaps));

for i = 1:length(timeGaps)
    gap = timeGaps(i);
    dt = gap * frameTime;
    
    % Pseudo-Kalman variance scaling
    P_predicted = 2 * Drot * dt;           % Process variance (grows with time)
    R_meas = 2 * orientCRLB^2;             % Measurement variance (doubled for both endpoints)
    S_total = P_predicted + R_meas;
    
    % Orientation cost
    C_orient = dAngle_rad.^2 ./ S_total;
    
    plot(dAngle_deg, C_orient, 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('Gap = %d frames (dt=%.0f ms)', gap, dt*1000));
end

xlabel('Angular Distance (°)');
ylabel('Orientation Cost');
title(sprintf('Cost vs Angular Distance\n(D_{rot}=%.3f rad²/s, σ_{orient}=%.1f°)', ...
    Drot, orientCRLB*180/pi));
legend('Location', 'northwest');
grid on;
xlim([0 90]);

% Add annotation explaining the effect
text(50, max(ylim)*0.8, {'Larger gaps →', 'Lower cost for', 'same angle'}, ...
    'FontSize', 10, 'Color', [0.5 0.5 0.5]);

hold off;

%% Panel 2: "Effective Tolerance" vs Time Gap
subplot(1, 3, 2);
hold on;

gaps = 1:20;  % Frame gaps
dt_vals = gaps * frameTime;

% Calculate the angular distance that gives cost = 1 (arbitrary threshold)
costThreshold = 1;
tolerance_deg = zeros(size(gaps));

for i = 1:length(gaps)
    dt = dt_vals(i);
    P_predicted = 2 * Drot * dt;
    R_meas = 2 * orientCRLB^2;
    S_total = P_predicted + R_meas;
    
    % Solve: dAngle² / S_total = costThreshold
    % → dAngle = sqrt(costThreshold * S_total)
    tolerance_rad = sqrt(costThreshold * S_total);
    tolerance_deg(i) = tolerance_rad * 180/pi;
end

% Also show sqrt(MSD) expectation from pure diffusion
expected_angle_deg = sqrt(2 * Drot * dt_vals) * 180/pi;

plot(gaps, tolerance_deg, 'b-', 'LineWidth', 2, 'DisplayName', 'Pseudo-Kalman tolerance (C=1)');
plot(gaps, expected_angle_deg, 'r--', 'LineWidth', 2, 'DisplayName', 'Expected √MSD from diffusion');

xlabel('Time Gap (frames)');
ylabel('Angular Tolerance (°)');
title(sprintf('Effective Angular Tolerance vs Gap\n(for cost threshold = %.1f)', costThreshold));
legend('Location', 'northwest');
grid on;
xlim([1 20]);

% Add annotation
text(12, max(tolerance_deg)*0.5, {'Tolerance grows', 'with √time'}, ...
    'FontSize', 10, 'Color', [0.5 0.5 0.5]);

hold off;

%% Panel 3: Comparison of Direct vs Pseudo-Kalman modes
subplot(1, 3, 3);
hold on;

% Parameters for comparison
Bo = 1.0;  % Bo normalization factor (from diffusion)
wOrient = 1.0;

dAngle_deg = linspace(0, 45, 100);
dAngle_rad = dAngle_deg * pi/180;

% Gap = 1 frame
dt_1 = 1 * frameTime;
P_1 = 2 * Drot * dt_1;
R_1 = 2 * orientCRLB^2;
S_1 = P_1 + R_1;

% Gap = 5 frames
dt_5 = 5 * frameTime;
P_5 = 2 * Drot * dt_5;
R_5 = 2 * orientCRLB^2;
S_5 = P_5 + R_5;

% Direct mode (no time scaling)
C_direct = wOrient * Bo * dAngle_rad.^2;

% Pseudo-Kalman mode
C_kalman_1 = wOrient * dAngle_rad.^2 ./ S_1;
C_kalman_5 = wOrient * dAngle_rad.^2 ./ S_5;

plot(dAngle_deg, C_direct, 'k-', 'LineWidth', 2, 'DisplayName', 'Direct mode (no time scaling)');
plot(dAngle_deg, C_kalman_1, 'b-', 'LineWidth', 2, 'DisplayName', 'Pseudo-Kalman (gap=1)');
plot(dAngle_deg, C_kalman_5, 'r-', 'LineWidth', 2, 'DisplayName', 'Pseudo-Kalman (gap=5)');

xlabel('Angular Distance (°)');
ylabel('Orientation Cost');
title('Direct vs Pseudo-Kalman Modes');
legend('Location', 'northwest');
grid on;
xlim([0 45]);

% Add key insight
annotation('textbox', [0.72, 0.15, 0.25, 0.12], ...
    'String', {'Key insight:', 'Pseudo-Kalman automatically', 'relaxes penalty for larger gaps'}, ...
    'FontSize', 9, 'BackgroundColor', [1 1 0.9], 'EdgeColor', [0.8 0.8 0.8]);

hold off;

%% Add super title
sgtitle(sprintf('Pseudo-Kalman Orientation Cost Scaling\nD_{rot}=%.3f rad²/s, Frame time=%.0f ms, σ_{orient}=%.1f°', ...
    Drot, frameTime*1000, orientCRLB*180/pi), 'FontWeight', 'bold');

%% Print summary to console
fprintf('\n=== PSEUDO-KALMAN ORIENTATION SCALING SUMMARY ===\n\n');
fprintf('Parameters:\n');
fprintf('  D_rot = %.4f rad²/s\n', Drot);
fprintf('  Frame time = %.0f ms\n', frameTime*1000);
fprintf('  σ_orient (CRLB) = %.2f° (%.4f rad)\n\n', orientCRLB*180/pi, orientCRLB);

fprintf('For gap = 1 frame (dt = %.0f ms):\n', frameTime*1000);
fprintf('  P_predicted = 2 × D_rot × dt = %.6f rad²\n', 2*Drot*frameTime);
fprintf('  R_measurement = 2 × σ² = %.6f rad²\n', 2*orientCRLB^2);
fprintf('  S_total = %.6f rad²\n', 2*Drot*frameTime + 2*orientCRLB^2);
fprintf('  → Cost for 10° mismatch: %.4f\n\n', (10*pi/180)^2 / (2*Drot*frameTime + 2*orientCRLB^2));

fprintf('For gap = 5 frames (dt = %.0f ms):\n', 5*frameTime*1000);
fprintf('  P_predicted = 2 × D_rot × dt = %.6f rad²\n', 2*Drot*5*frameTime);
fprintf('  R_measurement = 2 × σ² = %.6f rad²\n', 2*orientCRLB^2);
fprintf('  S_total = %.6f rad²\n', 2*Drot*5*frameTime + 2*orientCRLB^2);
fprintf('  → Cost for 10° mismatch: %.4f\n\n', (10*pi/180)^2 / (2*Drot*5*frameTime + 2*orientCRLB^2));

fprintf('INTERPRETATION:\n');
fprintf('  The pseudo-Kalman mode makes the cost SMALLER for larger gaps,\n');
fprintf('  reflecting that we EXPECT more angular change over longer times.\n');
fprintf('  This prevents the tracker from over-penalizing normal diffusion\n');
fprintf('  during gap closing.\n\n');

end
