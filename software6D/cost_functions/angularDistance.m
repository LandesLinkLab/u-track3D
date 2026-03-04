function [dAngle, dOmega, dTotal] = angularDistance(theta1, phi1, omega1, ...
    theta2, phi2, omega2, wOrient, wOmega, useGammaWeight)
%ANGULARDISTANCE Compute 6D orientation distance between dipole orientations.
%
%   Combines two terms:
%     1. Angular distance between mean dipole directions (theta, phi)
%     2. Wobble angle similarity (Omega)
%   with optional gamma-confidence weighting.
%
%   Accounts for 180-degree head-tail symmetry of fluorescent dipoles:
%   the dipole direction mu_hat and -mu_hat are physically equivalent,
%   so the distance metric satisfies d(mu_hat, -mu_hat) = 0.
%
%   SYNOPSIS:
%       [dAngle, dOmega, dTotal] = angularDistance(theta1, phi1, omega1, ...
%           theta2, phi2, omega2, wOrient, wOmega, useGammaWeight)
%
%   INPUT:
%       theta1, theta2    : Polar angle in radians, range [0, pi/2].
%                           (Measured from the optical axis.)
%       phi1, phi2        : Azimuthal angle in radians, range [-pi, pi].
%       omega1, omega2    : Wobble solid angle in steradians, range [0, 2*pi].
%                           (Set to [] or NaN to ignore wobble.)
%       wOrient           : Weight for angular distance (default: 1.0).
%       wOmega            : Weight for wobble distance (default: 0.5).
%       useGammaWeight    : If true, weight D_total by average gamma
%                           confidence factor (default: false).
%
%   OUTPUT:
%       dAngle : Angular distance between mean dipole orientations.
%                Range [0, pi/2]. Accounts for 180-degree symmetry.
%       dOmega : Absolute difference in wobble angle |Omega1 - Omega2|.
%                Range [0, 2*pi]. Returns 0 if omega values are empty/NaN.
%       dTotal : Combined orientation distance:
%                  dTotal = wOrient * dAngle + wOmega * dOmega
%                If useGammaWeight = true:
%                  dTotal *= gamma_avg = (gamma1 + gamma2) / 2
%                where gamma = 1 - 3*Omega/(4*pi) + Omega^2/(8*pi^2)
%
%   FORMULA:
%       D_angle  = acos(|mu_hat_1 . mu_hat_2|)   in [0, pi/2]
%       D_wobble = |Omega_1 - Omega_2|            in [0, 2*pi]
%       D_ori    = w_angle * D_angle + w_wobble * D_wobble
%
%       mu_hat = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
%
%       Optional gamma-confidence weighting:
%         gamma = 1 - 3*Omega/(4*pi) + Omega^2/(8*pi^2)
%         gamma ~ 1: Fixed dipole (Omega ~ 0) -> high confidence
%         gamma ~ 0: Large wobble -> low confidence -> reduce penalty
%         D_ori *= (gamma_1 + gamma_2) / 2
%
%   PHYSICAL INTERPRETATION:
%       - D_angle:  How different are the mean dipole directions?
%                   Same molecule should point similarly frame-to-frame.
%       - D_wobble: How different are the wobble cone sizes?
%                   Same molecule should have consistent Omega.
%       - gamma weighting: When Omega is large, the mean orientation is
%                   uncertain, so don't penalize angular changes as heavily.
%
%   All inputs can be arrays of the same size for vectorized computation.
%
%   REFERENCES:
%       Backer, A.S. & Moerner, W.E. (2014). J. Phys. Chem. B.
%       Jaqaman, K. et al. (2008). Nature Methods 5, 695-698.
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

% --- Default weights ---
if nargin < 7 || isempty(wOrient)
    wOrient = 1.0;
end
if nargin < 8 || isempty(wOmega)
    wOmega = 0.5;
end
if nargin < 9 || isempty(useGammaWeight)
    useGammaWeight = false;
end

% --- Convert to unit vectors on the sphere ---
% mu_hat = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
nx1 = sin(theta1) .* cos(phi1);
ny1 = sin(theta1) .* sin(phi1);
nz1 = cos(theta1);

nx2 = sin(theta2) .* cos(phi2);
ny2 = sin(theta2) .* sin(phi2);
nz2 = cos(theta2);

% --- Dot product ---
cosAngle = nx1.*nx2 + ny1.*ny2 + nz1.*nz2;

% --- Enforce 180-degree symmetry: d(mu_hat, -mu_hat) = 0 ---
% Taking |cos| maps both mu_hat and -mu_hat to the same distance.
cosAngle = abs(cosAngle);

% --- Clamp to valid range to avoid complex values from numerical noise ---
cosAngle = min(max(cosAngle, 0), 1);

% --- Angular distance in [0, pi/2] ---
dAngle = acos(cosAngle);

% --- Wobble distance ---
if nargin >= 6 && ~isempty(omega1) && ~isempty(omega2)
    % Handle NaN entries (missing wobble data)
    dOmega = abs(omega1 - omega2);
    dOmega(isnan(omega1) | isnan(omega2)) = 0;
else
    dOmega = zeros(size(dAngle));
end

% --- Weighted total ---
dTotal = wOrient .* dAngle + wOmega .* dOmega;

% --- Optional gamma-confidence weighting ---
% gamma = 1 - 3*Omega/(4*pi) + Omega^2/(8*pi^2)
% gamma ~ 1 means fixed dipole (high confidence in orientation)
% gamma ~ 0 means large wobble (low confidence -> reduce penalty)
if useGammaWeight && nargin >= 6 && ~isempty(omega1) && ~isempty(omega2)
    gamma1 = 1 - 3*omega1./(4*pi) + omega1.^2./(8*pi^2);
    gamma2 = 1 - 3*omega2./(4*pi) + omega2.^2./(8*pi^2);
    
    % Clamp gamma to [0, 1]
    gamma1 = max(0, min(1, gamma1));
    gamma2 = max(0, min(1, gamma2));
    
    % Average gamma as confidence weight
    gammaAvg = (gamma1 + gamma2) / 2;
    
    % Handle NaN
    gammaAvg(isnan(gammaAvg)) = 1;
    
    dTotal = gammaAvg .* dTotal;
end

end
