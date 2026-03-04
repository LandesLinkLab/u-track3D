function [traj] = DirectedMotion(t, D, kdm1, kdm2)
%DIRECTEDMOTION Simulate directed motion with fluctuating direction.
%
%   Combines Brownian diffusion with directed transport. The direction of
%   transport evolves over time via a Wiener process, creating realistic
%   trajectories that are not perfectly straight.
%
%   Speed is parameterized as: speed = k · √D
%   Direction angles (θ, φ) evolve via Wiener processes.
%
%   Reference:
%     Chatterjee, et al. "Feature Selection and Hyperparameter Optimization
%     for Machine Learned Classification of 3D Single-Particle Tracking."
%     Chemical & Biomedical Imaging (2025).
%
%   INPUT:
%     t     : [T x 1] Time vector in seconds
%     D     : Diffusion coefficient (µm²/s in original, we use nm²/s)
%     kdm1  : Minimum speed coefficient (speed = k·√D)
%     kdm2  : Maximum speed coefficient
%
%   OUTPUT:
%     traj  : [T-1 x 4] Trajectory matrix [x, y, z, time]
%
%   NOTES:
%     - k ≈ 1-5: weak directed motion (diffusion dominates)
%     - k ≈ 10: balanced directed motion
%     - k ≈ 20-50: strong directed motion
%     - Direction fluctuations use variance D/100
%
%   EXAMPLE:
%     t = linspace(0, 3, 100)';
%     traj = DirectedMotion(t, 0.1, 10, 10);  % Fixed k = 10
%
%   Emil Gillett & Jagriti Chatterjee, Landes Research Group, UIUC, 2026.

%% Initialize
dt = diff(t);

%% Random speed coefficient in range [kdm1, kdm2]
coef = kdm1 + (kdm2 - kdm1) * rand;

%% Speed module calculation: speed = k · √D
speed = coef * sqrt(D);

%% Initialize with Brownian diffusion base
traj = NormalDiffusion_DM(t, D);

%% Initial direction (random unit vector)
vel = randn([3, 1]);
v = vel ./ norm(vel);

% Convert to spherical coordinates
phi = atan2(v(2), v(1));          % Azimuthal angle
theta = acos(v(3) / sqrt(sum(v.^2)));  % Polar angle

%% Direction evolution via Wiener process
% Angular velocities with variance D/100
omega_phi = Wiener(t, 0, D/100);
omega_theta = Wiener(t, 0, D/100);

% Cumulative direction change
phi = phi + cumsum(omega_phi .* dt);
theta = theta + cumsum(omega_theta .* dt);

%% Compute velocity vector from evolving direction
% vel = speed · [cos(φ)sin(θ), sin(φ)sin(θ), cos(θ)]
vel = speed * [cos(phi) .* sin(theta), ...
               sin(phi) .* sin(theta), ...
               cos(theta)];

%% Add directed motion to Brownian base
traj(:, 1:3) = traj(:, 1:3) + cumsum(vel .* dt, 1);

end

%% Local helper: Normal diffusion
function [traj] = NormalDiffusion_DM(t, D)
%NORMALDIFFUSION_DM Simple 3D Brownian diffusion

N = length(t);
dt = diff(t);

traj = zeros([N-1, 3+1]);
traj(:, end) = t(2:end);

% Normal distribution: σ = √(2·D·dt)
mu = 0;
sigma = (2*D*dt).^0.5;

du_x = mu + sigma .* randn(N-1, 1);
du_y = mu + sigma .* randn(N-1, 1);
du_z = mu + sigma .* randn(N-1, 1);

du_3d = [du_x, du_y, du_z];
traj(:, 1:3) = cumsum(du_3d);

end

%% Local helper: Wiener process
function [res] = Wiener(time, drift, variance)
%WIENER Simulates a Wiener process (non-differentiable random motion)
%
%   INPUT:
%     time     : Time vector
%     drift    : Drift coefficient (typically 0)
%     variance : Variance coefficient
%
%   OUTPUT:
%     res      : [N-1 x 1] Cumulative Wiener process values

N = length(time);
X = randn(N-1, 1);          % Gaussian increments
dt = diff(time(:));         % Time steps

vel = drift .* dt + variance .* X .* sqrt(dt);
res = cumsum(vel);

% Ensure single-column output
if size(res, 2) ~= 1
    error('Wiener process output is not a single-column vector.');
end

end
