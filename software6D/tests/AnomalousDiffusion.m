function [traj] = AnomalousDiffusion(t, D, alphamin, alphamax)
%ANOMALOUSDIFFUSION Simulate anomalous (sub-diffusive) motion using Weierstrass-Mandelbrot function.
%
%   Uses the Weierstrass-Mandelbrot function as described by Saxton for
%   sub-diffusive motion simulation. The sum is taken from n = -8 to +48.
%
%   Reference:
%     Chatterjee, et al. "Feature Selection and Hyperparameter Optimization
%     for Machine Learned Classification of 3D Single-Particle Tracking."
%     Chemical & Biomedical Imaging (2025).
%
%   INPUT:
%     t         : [T x 1] Time vector in seconds
%     D         : Diffusion coefficient (µm²/s in original, we use nm²/s)
%     alphamin  : Minimum anomalous exponent (0 < α < 1 for sub-diffusion)
%     alphamax  : Maximum anomalous exponent
%
%   OUTPUT:
%     traj      : [T-1 x 4] Trajectory matrix [x, y, z, time]
%                 Positions are displacements from origin (cumulative)
%
%   NOTES:
%     - α < 1: sub-diffusion (slower than Brownian)
%     - α = 1: normal Brownian diffusion
%     - α > 1: super-diffusion (faster than Brownian)
%     - For typical sub-diffusive motion, use α ≈ 0.4-0.7
%
%   EXAMPLE:
%     t = linspace(0, 3, 100)';
%     traj = AnomalousDiffusion(t, 0.1, 0.4, 0.4);  % Fixed α = 0.4
%
%   Emil Gillett & Jagriti Chatterjee, Landes Research Group, UIUC, 2026.

%% Initialize
T = length(t);          % T - Number of positions, T-1 - number of displacements
dt = diff(t);           % Time differential between positions
alpha = alphamin + (alphamax - alphamin) * rand;  % Random α in range

traj = zeros([T-1, 3 + 1]);  % Initialize the trajectory [x, y, z, time]
traj(:, end) = t(2:end, :);  % Fill in time column

%% Weierstrass-Mandelbrot parameters
n = -8:48;              % Summation range as described by Saxton
gamma = sqrt(pi);       % Scaling factor
t_ = 2*pi / max(t) * t(2:end, :);  % Normalized time

phi = 2*pi * rand([3, length(n)]);  % Random phase for each dimension and n

%% Evaluate Weierstrass-Mandelbrot function
% W(t) = Σ_n [cos(φ) - cos(t·γⁿ + φ)] · γ^(αn/2)
for d = 1:3
    % Numerator: cos(φ) - cos(t·γⁿ + φ)
    num = cos(phi(d, :)) - cos(t_ * gamma.^n + phi(d, :));
    
    % Denominator: γ^(-αn/2)
    den = gamma.^(-alpha * n / 2);
    
    % Sum over all n
    W = sum(num .* den, 2);  % W(t) for this dimension
    
    % Append to trajectory
    traj(:, d) = W;
end

%% Rescale such that ⟨r₁²⟩ = 6 D dt
% This ensures the trajectory has the correct diffusion coefficient scaling
sqdisp = mean(sum((traj(2:end, 1:3) - traj(1:end-1, 1:3)).^2, 2));
if sqdisp > 0
    traj(:, 1:3) = traj(:, 1:3) * sqrt(6 * D * mean(dt) / sqdisp);
end

end
