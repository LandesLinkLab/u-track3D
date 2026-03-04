function [data, groundTruth] = generate6DTrajectories(varargin)
%GENERATE6DTRAJECTORIES Generate synthetic 6D-SMOLM trajectories in Fusion output format.
%
%   Produces 8-column data matching the mDeepSTORM3D + mDeep-SMOLM fusion
%   output: [frame, signal, x_nm, y_nm, z_nm, theta_deg, phi_deg, gamma].
%
%   Supports multiple spatial and rotational motion models, realistic
%   photoblinking, and localization noise. Use this to test the full
%   tracking pipeline before running on real experimental data.
%
%   SYNOPSIS:
%       [data, groundTruth] = generate6DTrajectories()
%       [data, groundTruth] = generate6DTrajectories('Name', Value, ...)
%
%   NAME-VALUE PARAMETERS:
%
%     --- Simulation geometry ---
%     'NumParticles'    : Number of molecules. Default: 5.
%     'NumFrames'       : Total frames in movie. Default: 500.
%     'FieldOfView'     : [xMin xMax yMin yMax] in nm. Default: [0 6000 0 6000].
%     'ZRange'          : [zMin zMax] in nm. Default: [-400 400].
%
%     --- Frame timing ---
%     'FrameTime'       : Exposure time per frame in seconds. Default: 0.05.
%
%     --- Spatial motion model ---
%     'SpatialModel'    : 'brownian' | 'confined' | 'directed' | 'hop'
%                         Default: 'brownian'.
%       Brownian:  free 3D diffusion.
%       Confined:  Brownian within a spherical confinement radius.
%       Directed:  Brownian + constant drift velocity.
%       Hop:       Brownian within a confined region, with occasional
%                  hops to a new center (models surface diffusion with
%                  trapping sites).
%     'Dtrans'          : Translational diffusion coefficient (nm^2/s).
%                         Default: 500. Can be a vector [nParticles x 1]
%                         for heterogeneous populations.
%     'DriftVelocity'   : [vx, vy, vz] in nm/s (for 'directed'). Default: [0 0 0].
%     'ConfinementRadius': Radius in nm (for 'confined'/'hop'). Default: 200.
%     'HopRate'         : Probability of hopping per frame (for 'hop'). Default: 0.01.
%     'HopDistance'     : Mean hop distance in nm (for 'hop'). Default: 300.
%
%     --- Rotational motion model ---
%     'RotationalModel' : 'brownian' | 'confined' | 'fixed'
%                         Default: 'brownian'.
%       Brownian:  Free rotational diffusion of the dipole orientation
%                  on the unit sphere. Uses the geometric integration
%                  scheme of Höfling & Straube, Phys. Rev. Research 7,
%                  043034 (2025). Random rotations are applied in the
%                  tangent plane via Rodrigues' formula. This avoids
%                  coordinate singularities at the poles, preserves
%                  |u|=1 exactly, and satisfies detailed balance.
%       Confined:  Rotational diffusion with a restoring torque toward
%                  a preferred axis (center of ThetaRange, PhiRange).
%                  The torque T_ext = kappa * (u x n_pref) derives from
%                  potential V(u) = -kappa * (u . n_pref), yielding
%                  the von Mises-Fisher equilibrium distribution.
%                  Reference: Höfling & Straube (2025), Eq. (28).
%       Fixed:     Orientation doesn't change (static dipole).
%     'Drot'            : Rotational diffusion coefficient (rad^2/s).
%                         Default: 0.05. Can be a vector.
%                         Rotational correlation time: tau_rot = 1/(6*Drot).
%                         For slow biomolecules in crowded environments,
%                         typical Drot ~ 0.01-0.1 rad^2/s gives tau_rot
%                         ~ 1-17 s, producing correlated orientations
%                         across ms-timescale frames.
%     'ThetaRange'      : Allowed theta range in degrees. Default: [0 90].
%                         Dipole head-tail symmetry is enforced when
%                         ThetaRange is within [0, 90].
%     'PhiRange'        : Allowed phi range in degrees. Default: [-180 180].
%
%     --- Wobble (gamma) model ---
%     'GammaModel'      : 'constant' | 'fluctuating'
%                         Default: 'constant'.
%     'GammaMean'       : Mean gamma value [0,1]. Default: 0.7.
%     'GammaStd'        : Std of gamma fluctuations (for 'fluctuating'). Default: 0.05.
%
%     --- Photophysics ---
%     'Signal'          : Mean photon count. Default: 1500.
%     'SignalStd'       : Std of photon count (shot noise). Default: 300.
%     'BlinkOnRate'     : Probability of turning on per frame (when dark). Default: 0.7.
%     'BlinkOffRate'    : Probability of turning off per frame (when on). Default: 0.1.
%     'BleachRate'      : Probability of permanent bleaching per frame. Default: 0.0005.
%
%     --- Localization noise ---
%     'AddNoise'        : Add localization noise to output. Default: false.
%     'NoiseXY'         : Localization noise std in nm (x,y). Default: 0.
%     'NoiseZ'          : Localization noise std in nm (z). Default: 0.
%     'NoiseTheta'      : Angular noise std in degrees (theta). Default: 0.
%     'NoisePhi'        : Angular noise std in degrees (phi). Default: 0.
%     'NoiseGamma'      : Gamma noise std. Default: 0.
%
%     --- Output ---
%     'OutputFile'      : Path to save .mat file. Default: auto-saves to
%                         the tests/ folder with a descriptive filename like
%                         sim_brownian_brownian_5part_500fr_20260209_143022.mat
%     'Seed'            : Random seed for reproducibility. Default: [] (no set).
%     'Verbose'         : Print summary. Default: true.
%
%   OUTPUT:
%     data        : [N x 8] matrix in Fusion output format:
%                   [frame, signal, x_nm, y_nm, z_nm, theta_deg, phi_deg, gamma]
%                   This includes localization noise — what the detector sees.
%
%     groundTruth : Struct with fields:
%       .trajectories : Cell array {nParticles x 1}. Each cell is [T x 10]:
%                       [frame, signal, x, y, z, theta_deg, phi_deg, gamma, trackID, isInterp]
%                       These are the TRUE positions (no noise).
%       .params       : Copy of all simulation parameters.
%       .Dtrans       : [nParticles x 1] true D_trans per particle.
%       .Drot         : [nParticles x 1] true D_rot per particle.
%
%   EXAMPLES:
%     % Simple Brownian
%     data = generate6DTrajectories('NumParticles', 3, 'Dtrans', 300);
%
%     % Confined diffusion with low wobble
%     data = generate6DTrajectories('SpatialModel', 'confined', ...
%         'ConfinementRadius', 150, 'GammaMean', 0.9);
%
%     % Directed motion with fast rotation
%     data = generate6DTrajectories('SpatialModel', 'directed', ...
%         'DriftVelocity', [50 0 0], 'Drot', 0.5);
%
%     % Heterogeneous population
%     data = generate6DTrajectories('NumParticles', 4, ...
%         'Dtrans', [100; 500; 1000; 2000], 'Drot', [0.01; 0.05; 0.1; 0.5]);
%
%   REFERENCES:
%     Rotational diffusion uses the geometric integration scheme of:
%       Höfling, F. & Straube, A. V. (2025). Rotational Brownian motion:
%       trajectory, orientation dynamics, and geometric integration.
%       Phys. Rev. Research, 7, 043034.
%       https://doi.org/10.1103/PhysRevResearch.7.043034
%
%   Emil Gillett, Landes Research Group, UIUC, 2026.

%% ========================================================================
%  PARSE INPUTS
%  ========================================================================
p = inputParser;

% Geometry
addParameter(p, 'NumParticles', 5, @isnumeric);
addParameter(p, 'NumFrames', 500, @isnumeric);
addParameter(p, 'FieldOfView', [0 6000 0 6000], @isnumeric);
addParameter(p, 'ZRange', [-400 400], @isnumeric);

% Timing
addParameter(p, 'FrameTime', 0.05, @isnumeric);

% Spatial motion
addParameter(p, 'SpatialModel', 'brownian', @ischar);
addParameter(p, 'Dtrans', 500, @isnumeric);
addParameter(p, 'DriftVelocity', [0 0 0], @isnumeric);
addParameter(p, 'ConfinementRadius', 200, @isnumeric);
addParameter(p, 'HopRate', 0.01, @isnumeric);
addParameter(p, 'HopDistance', 300, @isnumeric);

% Rotational motion
addParameter(p, 'RotationalModel', 'brownian', @ischar);
addParameter(p, 'Drot', 0.05, @isnumeric);
addParameter(p, 'ThetaRange', [0 90], @isnumeric);
addParameter(p, 'PhiRange', [-180 180], @isnumeric);

% Wobble
addParameter(p, 'GammaModel', 'constant', @ischar);
addParameter(p, 'GammaMean', 0.7, @isnumeric);
addParameter(p, 'GammaStd', 0.05, @isnumeric);

% Photophysics
addParameter(p, 'Signal', 1500, @isnumeric);
addParameter(p, 'SignalStd', 300, @isnumeric);
addParameter(p, 'BlinkOnRate', 0.7, @isnumeric);
addParameter(p, 'BlinkOffRate', 0.1, @isnumeric);
addParameter(p, 'BleachRate', 0.0005, @isnumeric);

% Localization noise
addParameter(p, 'AddNoise', false, @islogical);
addParameter(p, 'NoiseXY', 0, @isnumeric);
addParameter(p, 'NoiseZ', 0, @isnumeric);
addParameter(p, 'NoiseTheta', 0, @isnumeric);
addParameter(p, 'NoisePhi', 0, @isnumeric);
addParameter(p, 'NoiseGamma', 0, @isnumeric);

% Starting positions
addParameter(p, 'StartSpread', [], @isnumeric);  % If set, particles start within this radius of FOV center

% Output
addParameter(p, 'OutputFile', '', @ischar);
addParameter(p, 'Seed', [], @isnumeric);
addParameter(p, 'Verbose', true, @islogical);

parse(p, varargin{:});
o = p.Results;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
if ~isempty(o.Seed)
    rng(o.Seed);
end

nPart = o.NumParticles;
nFrames = o.NumFrames;
dt = o.FrameTime;
fov = o.FieldOfView;
zRange = o.ZRange;

% Expand scalar D values to per-particle vectors
if isscalar(o.Dtrans), Dtrans = repmat(o.Dtrans, nPart, 1);
else, Dtrans = o.Dtrans(:); end

if isscalar(o.Drot), Drot = repmat(o.Drot, nPart, 1);
else, Drot = o.Drot(:); end

assert(length(Dtrans) == nPart, 'Dtrans length must match NumParticles');
assert(length(Drot) == nPart, 'Drot length must match NumParticles');

%% ========================================================================
%  GENERATE TRAJECTORIES
%  ========================================================================

allData = [];           % Observed data (with noise)
trajCells = cell(nPart, 1);  % Ground truth per particle

for iPart = 1:nPart
    
    % --- Initial state ---
    if ~isempty(o.StartSpread)
        % Start particles within StartSpread radius of FOV center
        fovCenterX = (fov(1) + fov(2)) / 2;
        fovCenterY = (fov(3) + fov(4)) / 2;
        zCenterVal = (zRange(1) + zRange(2)) / 2;
        
        % Random offset within StartSpread for this particle
        angle = 2 * pi * rand;
        radius = o.StartSpread * sqrt(rand);  % sqrt for uniform distribution in disk
        x = fovCenterX + radius * cos(angle);
        y = fovCenterY + radius * sin(angle);
        z = zCenterVal + (rand - 0.5) * min(o.StartSpread, zRange(2) - zRange(1));
    else
        % Original behavior: random across entire FOV
        x = fov(1) + rand * (fov(2) - fov(1));
        y = fov(3) + rand * (fov(4) - fov(3));
        z = zRange(1) + rand * (zRange(2) - zRange(1));
    end
    
    theta_deg = o.ThetaRange(1) + rand * (o.ThetaRange(2) - o.ThetaRange(1));
    phi_deg   = o.PhiRange(1) + rand * (o.PhiRange(2) - o.PhiRange(1));
    gamma     = max(0, min(1, o.GammaMean + o.GammaStd * randn));
    
    % Confined diffusion center
    xCenter = x; yCenter = y; zCenter = z;
    
    % Photophysics state: 1 = on, 0 = dark, -1 = bleached
    photoState = 1;
    
    D_t = Dtrans(iPart);
    D_r = Drot(iPart);
    
    truthRows = [];
    
    for iFr = 1:nFrames
        
        %% --- Photophysics ---
        if photoState == -1
            continue;  % Permanently bleached
        end
        
        % Bleaching check
        if rand < o.BleachRate
            photoState = -1;
            continue;
        end
        
        % Blinking transitions
        if photoState == 1
            if rand < o.BlinkOffRate
                photoState = 0;
            end
        else  % photoState == 0
            if rand < o.BlinkOnRate
                photoState = 1;
            end
        end
        
        %% --- Spatial propagation (always, even when dark) ---
        sigma_step = sqrt(2 * D_t * dt);
        
        switch lower(o.SpatialModel)
            case 'brownian'
                x = x + sigma_step * randn;
                y = y + sigma_step * randn;
                z = z + sigma_step * randn;
                
            case 'confined'
                x = x + sigma_step * randn;
                y = y + sigma_step * randn;
                z = z + sigma_step * randn;
                
                % Reflect off confinement boundary
                R = o.ConfinementRadius;
                dr = sqrt((x-xCenter)^2 + (y-yCenter)^2 + (z-zCenter)^2);
                if dr > R
                    % Project back onto sphere surface
                    scale = R / dr;
                    x = xCenter + (x - xCenter) * scale;
                    y = yCenter + (y - yCenter) * scale;
                    z = zCenter + (z - zCenter) * scale;
                end
                
            case 'directed'
                vDrift = o.DriftVelocity(:)';
                x = x + sigma_step * randn + vDrift(1) * dt;
                y = y + sigma_step * randn + vDrift(2) * dt;
                z = z + sigma_step * randn + vDrift(3) * dt;
                
            case 'hop'
                x = x + sigma_step * randn;
                y = y + sigma_step * randn;
                z = z + sigma_step * randn;
                
                % Confine within current trap
                R = o.ConfinementRadius;
                dr = sqrt((x-xCenter)^2 + (y-yCenter)^2 + (z-zCenter)^2);
                if dr > R
                    scale = R / dr;
                    x = xCenter + (x - xCenter) * scale;
                    y = yCenter + (y - yCenter) * scale;
                    z = zCenter + (z - zCenter) * scale;
                end
                
                % Hop to new trap center
                if rand < o.HopRate
                    hopDir = randn(1, 3);
                    hopDir = hopDir / norm(hopDir);
                    hopDist = o.HopDistance * (0.5 + rand);
                    xCenter = xCenter + hopDir(1) * hopDist;
                    yCenter = yCenter + hopDir(2) * hopDist;
                    zCenter = zCenter + hopDir(3) * hopDist;
                    x = xCenter; y = yCenter; z = zCenter;
                end
                
            otherwise
                error('Unknown SpatialModel: %s', o.SpatialModel);
        end
        
        % Clamp to FOV and z-range
        x = max(fov(1), min(fov(2), x));
        y = max(fov(3), min(fov(4), y));
        z = max(zRange(1), min(zRange(2), z));
        
        %% --- Rotational propagation ---
        % Uses geometric rotational Brownian dynamics on S^2 following:
        %   Höfling & Straube, Phys. Rev. Research 7, 043034 (2025).
        %   "Rotational Brownian motion: trajectory, orientation dynamics,
        %    and geometric integration"
        %
        % The algorithm applies infinitesimal random rotations in the
        % tangent plane via Rodrigues' formula, avoiding coordinate
        % singularities at the poles and preserving |u|=1 exactly.
        
        switch lower(o.RotationalModel)
            case 'brownian'
                % Convert current orientation to unit vector
                u = orientation_to_uvec(theta_deg, phi_deg);
                
                % Apply geometric Brownian step on sphere
                u = rotate_brownian_step(u, D_r, dt);
                
                % Enforce dipole symmetry and angular bounds
                [u, theta_deg, phi_deg] = enforce_orientation_bounds(u, o.ThetaRange, o.PhiRange);
                
            case 'confined'
                % Orientational confinement via restoring potential.
                % Uses restoring torque T_ext = kappa * (u x n_pref)
                % derived from potential V(u) = -kappa * (u . n_pref).
                % Equilibrium: von Mises-Fisher distribution.
                %
                % Reference: Höfling & Straube (2025), Eq. (28) and Table I.
                
                % Convert current orientation to unit vector
                u = orientation_to_uvec(theta_deg, phi_deg);
                
                % Preferred axis: center of angular range
                n_pref_theta = mean(o.ThetaRange);
                n_pref_phi = mean(o.PhiRange);
                n_pref = orientation_to_uvec(n_pref_theta, n_pref_phi);
                
                % Dimensionless restoring strength kappa/(k_B T)
                % Higher values = tighter confinement to preferred axis
                kappa_over_kBT = 3.0;  % moderate confinement (~30-45 deg cone)
                
                % Effective torque: T_eff = (kappa/kBT) * D_r * (u x n_pref)
                T_eff = kappa_over_kBT * D_r * cross(u, n_pref);
                
                % Apply biased geometric Brownian step
                u = rotate_brownian_step_torque(u, D_r, dt, T_eff);
                
                % Enforce dipole symmetry and angular bounds
                [u, theta_deg, phi_deg] = enforce_orientation_bounds(u, o.ThetaRange, o.PhiRange);
                
            case 'fixed'
                % No change - orientation remains constant
                
            otherwise
                error('Unknown RotationalModel: %s', o.RotationalModel);
        end
        
        %% --- Wobble (gamma) propagation ---
        switch lower(o.GammaModel)
            case 'constant'
                % gamma stays at initial value (set above)
                
            case 'fluctuating'
                gamma = gamma + o.GammaStd * randn * 0.3;  % Slow fluctuation
                gamma = max(0, min(1, gamma));
                
            otherwise
                error('Unknown GammaModel: %s', o.GammaModel);
        end
        
        %% --- Record if emitting ---
        if photoState == 1
            % Signal with shot noise
            sig = max(100, o.Signal + o.SignalStd * randn);
            
            % Ground truth (no noise)
            truthRows = [truthRows; ...
                iFr, sig, x, y, z, theta_deg, phi_deg, gamma, iPart, 0]; %#ok<AGROW>
            
            % Observed data (with or without localization noise)
            if o.AddNoise
                x_obs     = x + o.NoiseXY * randn;
                y_obs     = y + o.NoiseXY * randn;
                z_obs     = z + o.NoiseZ * randn;
                theta_obs = max(0, min(90, theta_deg + o.NoiseTheta * randn));
                phi_obs   = mod(phi_deg + o.NoisePhi * randn + 180, 360) - 180;
                gamma_obs = max(0, min(1, gamma + o.NoiseGamma * randn));
            else
                % No noise - output matches ground truth
                x_obs = x;
                y_obs = y;
                z_obs = z;
                theta_obs = theta_deg;
                phi_obs = phi_deg;
                gamma_obs = gamma;
            end
            
            allData = [allData; ...
                iFr, sig, x_obs, y_obs, z_obs, theta_obs, phi_obs, gamma_obs]; %#ok<AGROW>
        end
        
    end  % frames
    
    trajCells{iPart} = truthRows;
    
end  % particles

%% ========================================================================
%  SORT BY FRAME (mimics real Fusion output)
%  ========================================================================
data = sortrows(allData, 1);

%% ========================================================================
%  GROUND TRUTH
%  ========================================================================
groundTruth.trajectories = trajCells;
groundTruth.params = o;
groundTruth.Dtrans = Dtrans;
groundTruth.Drot = Drot;
groundTruth.gammaToOmega = @(g) pi * (3 - sqrt(1 + 8*g));  % conversion function

%% ========================================================================
%  SAVE
%  ========================================================================

% Auto-generate output path if not specified
if isempty(o.OutputFile)
    thisDir = fileparts(mfilename('fullpath'));
    if isempty(thisDir), thisDir = pwd; end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    autoName = sprintf('sim_%s_%s_%dpart_%dfr_%s.mat', ...
        o.SpatialModel, o.RotationalModel, nPart, nFrames, timestamp);
    o.OutputFile = fullfile(thisDir, autoName);
end

est = data; %#ok<NASGU>  % "est" matches Fusion naming convention
save(o.OutputFile, 'est', 'groundTruth');

if o.Verbose
    fprintf('  Saved: %s\n', o.OutputFile);
end

%% ========================================================================
%  SUMMARY
%  ========================================================================
if o.Verbose
    fprintf('\n--- generate6DTrajectories summary ---\n');
    fprintf('  Particles:       %d\n', nPart);
    fprintf('  Frames:          %d\n', nFrames);
    fprintf('  Spatial model:   %s\n', o.SpatialModel);
    fprintf('  Rotational model: %s\n', o.RotationalModel);
    fprintf('  Total detections: %d (%.1f per frame avg)\n', ...
        size(data,1), size(data,1)/nFrames);
    
    % Per-particle stats
    for iPart = 1:nPart
        tr = trajCells{iPart};
        if ~isempty(tr)
            nDet = size(tr, 1);
            fracOn = nDet / nFrames;
            fprintf('  Particle %d: D_t=%.0f nm^2/s, D_r=%.3f rad^2/s, %d dets (%.0f%% on)\n', ...
                iPart, Dtrans(iPart), Drot(iPart), nDet, fracOn*100);
        else
            fprintf('  Particle %d: bleached/no detections\n', iPart);
        end
    end
    
    fprintf('  Output columns: [frame, signal, x_nm, y_nm, z_nm, theta_deg, phi_deg, gamma]\n');
    fprintf('--------------------------------------\n');
end

end

%% ========================================================================
%  GEOMETRIC ROTATIONAL BROWNIAN DYNAMICS HELPER FUNCTIONS
%  ========================================================================
%  Reference: Höfling & Straube, Phys. Rev. Research 7, 043034 (2025).
%  "Rotational Brownian motion: trajectory, orientation dynamics, and
%   geometric integration"
%
%  These functions implement singularity-free integration of rotational
%  Brownian motion on the unit sphere via infinitesimal random rotations
%  applied through Rodrigues' formula. This approach:
%    - Avoids coordinate singularities at the poles
%    - Preserves |u|=1 exactly
%    - Satisfies detailed balance
%  ========================================================================

function u = orientation_to_uvec(theta_deg, phi_deg)
%ORIENTATION_TO_UVEC  Convert (theta, phi) in degrees to a 3x1 unit vector.
%   theta: polar angle from +z axis [0, 90] deg (upper hemisphere)
%   phi:   azimuthal angle [-180, 180] deg
    th = theta_deg * pi/180;
    ph = phi_deg   * pi/180;
    u  = [sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];
end

function [u, theta_out, phi_out] = enforce_orientation_bounds(u, theta_range, phi_range)
%ENFORCE_ORIENTATION_BOUNDS  Reflect orientation into permissible ranges.
%   Uses reflecting boundaries in (theta, phi) to confine the orientation
%   to user-specified angular ranges. If theta_range is entirely within
%   [0, 90], dipole head-tail symmetry is applied first (u_z < 0 -> flip).
%
%   Reflecting boundaries: when theta or phi exceeds a boundary, the
%   excess is "bounced" back. Multiple reflections are handled iteratively.
%
%   u           : 3x1 unit vector
%   theta_range : [theta_min theta_max] in degrees
%   phi_range   : [phi_min phi_max] in degrees

    % Step 1: Dipole head-tail symmetry only if range is within upper hemisphere
    if theta_range(2) <= 90 && u(3) < 0
        u = -u;
    end

    % Step 2: Convert to angles (full [0, 180] range for theta)
    theta = acosd(max(-1, min(1, u(3))));   % [0, 180]
    phi   = atan2d(u(2), u(1));             % [-180, 180]

    % Step 3: Reflect theta into theta_range
    th_lo = theta_range(1);
    th_hi = theta_range(2);
    for iter = 1:50
        if theta >= th_lo && theta <= th_hi, break; end
        if theta < th_lo
            theta = 2*th_lo - theta;   % reflect off lower bound
        elseif theta > th_hi
            theta = 2*th_hi - theta;   % reflect off upper bound
        end
    end
    theta = max(th_lo, min(th_hi, theta));

    % Step 4: Reflect phi into phi_range
    ph_lo = phi_range(1);
    ph_hi = phi_range(2);
    for iter = 1:50
        if phi >= ph_lo && phi <= ph_hi, break; end
        if phi < ph_lo
            phi = 2*ph_lo - phi;       % reflect off lower bound
        elseif phi > ph_hi
            phi = 2*ph_hi - phi;       % reflect off upper bound
        end
    end
    phi = max(ph_lo, min(ph_hi, phi));

    % Step 5: Reconstruct unit vector from reflected angles
    th_r = theta * pi/180;
    ph_r = phi   * pi/180;
    u = [sin(th_r)*cos(ph_r); sin(th_r)*sin(ph_r); cos(th_r)];
    theta_out = theta;
    phi_out   = phi;
end

function u_new = rotate_brownian_step(u, Drot, dt)
%ROTATE_BROWNIAN_STEP  One step of geometric rotational BD on S^2.
%   Applies a random rotation in the tangent plane of u with variance
%   2*Drot*dt per component, using Rodrigues' formula.
%
%   Reference: Höfling & Straube, Phys. Rev. Research 7, 043034 (2025),
%   Eq. (24)-(27) and Table I.
%
%   u    : current orientation (3x1 unit vector)
%   Drot : rotational diffusion coefficient [rad^2/s]
%   dt   : time step [s]

    sigma = sqrt(2 * Drot * dt);  % std dev of each tangent component [rad]

    % Build an orthonormal basis for the tangent plane at u.
    % Choose e1 perpendicular to u:
    if abs(u(3)) < 0.9
        e1 = cross(u, [0; 0; 1]);
    else
        e1 = cross(u, [1; 0; 0]);
    end
    e1 = e1 / norm(e1);
    e2 = cross(u, e1);  % e2 = u x e1, already unit length

    % Random rotation vector in tangent plane
    dOmega = sigma * randn * e1 + sigma * randn * e2;

    % Apply rotation via Rodrigues' formula
    angle = norm(dOmega);
    if angle < 1e-15
        u_new = u;
    else
        u_new = rodrigues_rotate(u, dOmega / angle, angle);
    end
end

function u_new = rotate_brownian_step_torque(u, Drot, dt, T_eff)
%ROTATE_BROWNIAN_STEP_TORQUE  Geometric rotational BD with external torque.
%   Same as rotate_brownian_step, but the tangent-plane Gaussian
%   coefficients have nonzero means determined by the effective angular
%   velocity T_eff = (1/zeta_R) * T_ext projected onto the tangent plane.
%
%   For orientational confinement with restoring potential
%   V(u) = -kappa * (u . n_pref):
%     T_eff = (kappa/kBT) * D_R * (u x n_pref)
%
%   Reference: Höfling & Straube, Phys. Rev. Research 7, 043034 (2025),
%   Eq. (28) and Table I.
%
%   u     : current orientation (3x1 unit vector)
%   Drot  : rotational diffusion coefficient [rad^2/s]
%   dt    : time step [s]
%   T_eff : effective angular velocity = zeta_R^{-1} * T_ext (3x1 vector)

    sigma = sqrt(2 * Drot * dt);  % std dev per tangent component [rad]

    % Build an orthonormal basis for the tangent plane at u.
    if abs(u(3)) < 0.9
        e1 = cross(u, [0; 0; 1]);
    else
        e1 = cross(u, [1; 0; 0]);
    end
    e1 = e1 / norm(e1);
    e2 = cross(u, e1);

    % Deterministic mean from torque: project T_eff onto tangent basis
    % <xi_i> = e_i . (zeta_R^{-1} T_ext) * dt  [Eq. (28)]
    mean1 = dot(e1, T_eff) * dt;
    mean2 = dot(e2, T_eff) * dt;

    % Draw biased Gaussian coefficients
    xi1 = mean1 + sigma * randn;
    xi2 = mean2 + sigma * randn;

    % Random rotation vector in tangent plane (with drift)
    dOmega = xi1 * e1 + xi2 * e2;

    % Apply rotation via Rodrigues' formula
    angle = norm(dOmega);
    if angle < 1e-15
        u_new = u;
    else
        u_new = rodrigues_rotate(u, dOmega / angle, angle);
    end
end

function v_rot = rodrigues_rotate(v, k, angle)
%RODRIGUES_ROTATE  Rotate vector v about unit axis k by angle [rad].
%   Rodrigues' rotation formula following the sign convention of
%   Höfling & Straube, Phys. Rev. Research 7, 043034 (2025), Eq. (10)-(12):
%     v_rot = v*cos(a) - (k x v)*sin(a) + k*(k.v)*(1-cos(a))
%   The minus sign on the cross-product term arises from their generator
%   matrices (J_i)_{jk} = -eps_{ijk}.
    ca = cos(angle);
    sa = sin(angle);
    v_rot = v * ca - cross(k, v) * sa + k * (dot(k, v)) * (1 - ca);
    % Enforce unit norm (guard against floating-point drift)
    v_rot = v_rot / norm(v_rot);
end
