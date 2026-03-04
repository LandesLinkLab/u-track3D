function [traj, B] = ConfinedDiffusion(t, D, Bmin, Bmax)
%CONFINEDDIFFUSION Simulate confined diffusion using B-parameter rejection sampling.
%
%   Simulates Brownian diffusion confined within a spherical region. Uses
%   rejection sampling with sub-steps to ensure physical confinement.
%
%   The confinement is characterized by the dimensionless B parameter:
%     B = V_ellipsoid / V_sphere
%   where the sphere radius is computed from diffusion parameters.
%
%   Reference:
%     Chatterjee, et al. "Feature Selection and Hyperparameter Optimization
%     for Machine Learned Classification of 3D Single-Particle Tracking."
%     Chemical & Biomedical Imaging (2025).
%
%   INPUT:
%     t     : [T x 1] Time vector in seconds
%     D     : Diffusion coefficient (µm²/s in original, we use nm²/s)
%     Bmin  : Minimum B parameter (dimensionless, typically 1-10)
%     Bmax  : Maximum B parameter
%
%   OUTPUT:
%     traj  : [T-1 x 4] Trajectory matrix [x, y, z, time]
%     B     : Actual B parameter (ratio of ellipsoid to sphere volume)
%
%   NOTES:
%     - Larger B = stronger confinement (smaller region)
%     - B ≈ 1: weak confinement
%     - B ≈ 3-5: moderate confinement
%     - B > 10: strong confinement
%     - Uses dt/100 sub-steps for accurate rejection sampling
%
%   EXAMPLE:
%     t = linspace(0, 3, 100)';
%     [traj, B] = ConfinedDiffusion(t, 0.1, 3, 3);  % Fixed B = 3
%
%   Emil Gillett & Jagriti Chatterjee, Landes Research Group, UIUC, 2026.

%% Initialize
N = length(t);              % Number of time points
dt = mean(diff(t));         % Time step for steps
ddt = mean(diff(t)) / 100;  % Smaller time step for sub-steps (dt/100)
t_mini = (0:ddt:dt)';       % Sub-step time vector

%% Randomly choose B between Bmin and Bmax
B_param = Bmin + rand() * (Bmax - Bmin);

%% Calculate the radius of confinement from B_param
% r = √(D·T·dt) / B^(1/3)
r = sqrt(D * (N-1) * dt) / B_param^(1/3);

%% Initialize trajectory
traj = zeros(N-1, 4);           % [x, y, z, time]
traj_mini_res = zeros(N-1, 3);  % Accumulated sub-step displacements

%% Main simulation loop with rejection sampling
k = 1;
maxAttempts = 1000;  % Prevent infinite loops

while k <= N-1
    attempts = 0;
    accepted = false;
    
    while ~accepted && attempts < maxAttempts
        % Generate sub-steps using normal diffusion
        traj_mini = NormalDiffusion_substep(t_mini, D);
        
        % Test displacement: add to current position
        traj_mini_res(k, :) = traj_mini(end, 1:3);
        traj_test = cumsum(traj_mini_res(1:k, :), 1);
        
        % Check if within confinement sphere
        len = sqrt(traj_test(end, 1)^2 + traj_test(end, 2)^2 + traj_test(end, 3)^2);
        
        if len <= r
            accepted = true;
            k = k + 1;  % Accept the step
        else
            % Reject: reset this step's contribution
            traj_mini_res(k, :) = [0, 0, 0];
            attempts = attempts + 1;
        end
    end
    
    % If max attempts reached, force acceptance with boundary reflection
    if ~accepted
        % Project back onto sphere surface
        traj_test = cumsum(traj_mini_res(1:k, :), 1);
        len = sqrt(traj_test(end, 1)^2 + traj_test(end, 2)^2 + traj_test(end, 3)^2);
        if len > r
            scale = r / len * 0.99;  % Slightly inside
            traj_mini_res(k, :) = traj_mini_res(k, :) * scale;
        end
        k = k + 1;
    end
end

%% Compute cumulative trajectory
traj(:, 1:3) = cumsum(traj_mini_res, 1);
traj(:, 4) = t(2:end);

%% Calculate actual B parameter from trajectory
% Radius of gyration as approximation for ellipsoid fitting
traj_center = mean(traj(:, 1:3), 1);
rg = sqrt(mean(sum((traj(:, 1:3) - traj_center).^2, 2)));  % Radius of gyration

% Approximate volume of smallest ellipsoid
V_ell = (4/3) * pi * rg^3;

% Define the ratio B
B = V_ell / ((4/3) * pi * r^3);

end

%% Local helper function for sub-step diffusion
function [traj] = NormalDiffusion_substep(t, D)
%NORMALDIFFUSION_SUBSTEP Simple Brownian diffusion for sub-steps

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
