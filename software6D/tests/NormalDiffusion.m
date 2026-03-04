function [traj] = NormalDiffusion(t, D)
%NORMALDIFFUSION Simulate standard 3D Brownian diffusion.
%
%   Generates a trajectory following the standard diffusion equation:
%     ⟨r²⟩ = 6·D·t  (3D mean squared displacement)
%
%   Each displacement is drawn from a Gaussian distribution:
%     Δx ~ N(0, √(2·D·dt))
%
%   Reference:
%     Chatterjee, et al. "Feature Selection and Hyperparameter Optimization
%     for Machine Learned Classification of 3D Single-Particle Tracking."
%     Chemical & Biomedical Imaging (2025).
%
%   INPUT:
%     t : [T x 1] Time vector in seconds
%     D : Diffusion coefficient (µm²/s in original, we use nm²/s)
%
%   OUTPUT:
%     traj : [T-1 x 4] Trajectory matrix [x, y, z, time]
%            Positions are cumulative displacements from origin
%
%   NOTES:
%     - For D in nm²/s, typical values: 100-10000 nm²/s
%     - For D in µm²/s, typical values: 0.01-1.0 µm²/s
%     - 1 µm²/s = 1,000,000 nm²/s
%
%   EXAMPLE:
%     t = linspace(0, 3, 100)';
%     traj = NormalDiffusion(t, 0.1);  % D = 0.1 µm²/s
%
%   Emil Gillett & Jagriti Chatterjee, Landes Research Group, UIUC, 2026.

%% Initialize
N = length(t);    % N positions → N-1 displacements
dt = diff(t);     % Time differential between positions

traj = zeros([N-1, 3+1]);  % [x, y, z, time]
traj(:, end) = t(2:end);   % Fill time column

%% Generate displacements from Gaussian distribution
% Standard deviation: σ = √(2·D·dt)
mu = 0;
sigma = (2 * D * dt).^0.5;

du_x = mu + sigma .* randn(N-1, 1);
du_y = mu + sigma .* randn(N-1, 1);
du_z = mu + sigma .* randn(N-1, 1);

du_3d = [du_x, du_y, du_z];

%% Cumulative sum gives positions
traj(:, 1:3) = cumsum(du_3d);

end
