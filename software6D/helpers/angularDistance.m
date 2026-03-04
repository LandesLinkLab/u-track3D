function [dAngle, dOmega, dTotal] = angularDistance(theta1, phi1, omega1, theta2, phi2, omega2, varargin)
%ANGULARDISTANCE Compute angular distance between dipole orientations.
%
%   Calculates the angular distance between dipole orientations accounting
%   for 180° head-tail symmetry. Optionally includes wobble (Omega) distance.
%
%   SYNOPSIS:
%       [dAngle, dOmega, dTotal] = angularDistance(theta1, phi1, omega1, ...
%                                                   theta2, phi2, omega2)
%       [...] = angularDistance(..., 'Name', Value, ...)
%
%   INPUT:
%       theta1, theta2 : Polar angles in radians [0, π/2]
%       phi1, phi2     : Azimuthal angles in radians [-π, π]
%       omega1, omega2 : Wobble angles in steradians [0, 2π]
%
%   NAME-VALUE PARAMETERS:
%       'wAngle'        : Weight for mean direction distance (default: 1.0)
%       'wOmega'        : Weight for wobble Ω difference (default: 0.5)
%       'useGammaWeight': Apply γ-confidence weighting (default: false)
%       'normalize'     : Normalize output to [0, 1] (default: false)
%
%   OUTPUT:
%       dAngle : Angular distance for mean dipole direction [0, π/2]
%       dOmega : Absolute difference in wobble angle
%       dTotal : Weighted combination: w_angle * dAngle + w_omega * dOmega
%
%   ALGORITHM:
%       For dipole orientations with head-tail symmetry (μ̂ ≡ -μ̂):
%           μ̂ = [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
%           dAngle = acos(|μ̂₁ · μ̂₂|)
%       Range: [0, π/2] radians (90° max due to symmetry)
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% --- Parse inputs ---
p = inputParser;
addParameter(p, 'wAngle', 1.0);
addParameter(p, 'wOmega', 0.5);
addParameter(p, 'useGammaWeight', false);
addParameter(p, 'normalize', false);
parse(p, varargin{:});

wAngle = p.Results.wAngle;
wOmega = p.Results.wOmega;
useGammaWeight = p.Results.useGammaWeight;
doNormalize = p.Results.normalize;

%% --- Compute dipole unit vectors ---
mu1_x = sin(theta1) .* cos(phi1);
mu1_y = sin(theta1) .* sin(phi1);
mu1_z = cos(theta1);

mu2_x = sin(theta2) .* cos(phi2);
mu2_y = sin(theta2) .* sin(phi2);
mu2_z = cos(theta2);

%% --- Compute dot product ---
dotProd = mu1_x .* mu2_x + mu1_y .* mu2_y + mu1_z .* mu2_z;

%% --- Angular distance with 180° symmetry ---
absDot = abs(dotProd);
absDot = min(max(absDot, 0), 1);
dAngle = acos(absDot);

%% --- Wobble distance ---
dOmega = abs(omega1 - omega2);

%% --- Optional gamma weighting ---
if useGammaWeight
    gamma1 = 1 - 3*omega1/(4*pi) + omega1.^2/(8*pi^2);
    gamma2 = 1 - 3*omega2/(4*pi) + omega2.^2/(8*pi^2);
    gammaAvg = (gamma1 + gamma2) / 2;
    dAngle = dAngle .* gammaAvg;
end

%% --- Total distance ---
dTotal = wAngle * dAngle + wOmega * dOmega;

%% --- Optional normalization ---
if doNormalize
    dAngle = dAngle / (pi/2);
    dOmega = dOmega / (2*pi);
    dTotal = dTotal / (wAngle * pi/2 + wOmega * 2*pi);
end

end
