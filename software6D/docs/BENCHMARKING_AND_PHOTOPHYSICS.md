# Relating Simulation Parameters to Experimental Photophysics & Benchmarking Guide

## Part 1: Relating BlinkOnRate/BlinkOffRate to Experimental Photoblinking

### The Physics of Photoblinking

Fluorescent molecules like Rhodamine B undergo **photoblinking** - stochastic transitions between emissive ("ON") and non-emissive ("OFF") states. The OFF state is typically:
- A **triplet state** (T₁) with lifetime ~1-10 ms
- A **radical dark state** formed through the triplet, with lifetime ~10 ms to seconds
- A **charge-separated state** from electron transfer

### Experimental Timescales for Rhodamine Dyes

| Parameter | Rhodamine 6G/B | Units | Source |
|-----------|----------------|-------|--------|
| Fluorescence lifetime | 1.5-4 ns | ns | Varies with solvent |
| Triplet lifetime | 1-10 ms | ms | In polymer matrices |
| Dark state lifetime | 10 ms - 10 s | ms-s | Oxygen-dependent |
| ON time (τ_on) | 100 ms - 10 s | ms-s | Environment-dependent |
| OFF time (τ_off) | 1 ms - 1 s | ms-s | Environment-dependent |

### Converting Experimental Values to Simulation Parameters

The simulation uses **per-frame probabilities**:
- `BlinkOffRate` = Probability of turning OFF per frame (when currently ON)
- `BlinkOnRate` = Probability of turning ON per frame (when currently OFF)

**Key relationship:**
```
BlinkOffRate ≈ FrameTime / τ_on    (probability of turning off per frame)
BlinkOnRate  ≈ FrameTime / τ_off   (probability of recovering per frame)
```

Where:
- `τ_on` = mean ON time before blinking OFF
- `τ_off` = mean OFF time (dark state duration)
- `FrameTime` = camera exposure/frame interval (default: 50 ms)

### Example Calculations

**For Rhodamine B in PVA (typical SMOLM conditions):**
- τ_on ≈ 500 ms (mean time in ON state)
- τ_off ≈ 70 ms (mean dark state duration)
- FrameTime = 50 ms

```matlab
BlinkOffRate = 50 ms / 500 ms = 0.10   % 10% chance to turn OFF per frame
BlinkOnRate  = 50 ms / 70 ms  = 0.71   % 71% chance to turn ON per frame
```

**These are the defaults in run_demo.m!**

### Parameter Regimes

| Regime | BlinkOffRate | BlinkOnRate | Physical Meaning |
|--------|--------------|-------------|------------------|
| Stable emitter | 0.01-0.05 | 0.8-0.95 | Long ON times, quick recovery |
| Moderate blinking | 0.10-0.20 | 0.5-0.7 | Typical organic dyes |
| Heavy blinking | 0.3-0.5 | 0.3-0.5 | STORM/PALM probes |
| Mostly dark | 0.5-0.8 | 0.1-0.2 | Photoactivatable probes |

### Factors Affecting Experimental Blinking

1. **Oxygen concentration**: Oxygen quenches triplet states → shorter τ_off → higher BlinkOnRate
2. **Excitation intensity**: Higher intensity → more triplet formation → lower BlinkOffRate (more blinking)
3. **Polymer matrix**: Restricts oxygen diffusion → longer dark states → lower BlinkOnRate
4. **Temperature**: Higher T → faster recovery → higher BlinkOnRate
5. **Reducing agents**: (e.g., BME, MEA in STORM) → control blinking rates

---

## Part 2: Benchmarking Parameters from ISBI Challenge & u-track3D

### ISBI Particle Tracking Challenge Parameters

The ISBI 2012 Challenge established standard benchmarking conditions:

#### Signal-to-Noise Ratio (SNR)
| Level | SNR Value | Description |
|-------|-----------|-------------|
| Very Low | 1 | Methods fail significantly |
| Low | 2 | Challenging detection |
| **Critical** | **4** | Threshold where methods break down |
| High | 7 | Good performance expected |

**SNR Definition:** Peak signal above background divided by noise standard deviation

#### Particle Density
| Level | Particles | Description |
|-------|-----------|-------------|
| Low | 50-100 | Well-separated, easy linking |
| Medium | 100-500 | Moderate ambiguity |
| High | 1000+ | Dense, significant overlap |

#### Motion Scenarios
1. **Vesicles**: Brownian motion, D = 0.1-1 µm²/s
2. **Receptors**: Confined diffusion
3. **Viruses**: Directed + Brownian (constant velocity + diffusion)
4. **Microtubule tips**: Linear motion with variable velocity

### u-track3D Benchmark Parameters

From the u-track3D paper (Cell Reports Methods, 2023):

| Parameter | Test Range | Notes |
|-----------|------------|-------|
| SNR | 1-7 | Critical threshold at 4 |
| Particle density | Low/Med/High | Scenario-dependent |
| Gap closing | 1-10 frames | Maximum gap to close |
| Minimum track length | 3-10 frames | Filter short tracks |

### Performance Metrics

| Metric | Description | Formula |
|--------|-------------|---------|
| **α (alpha)** | Detection precision | TP / (TP + FP) |
| **β (beta)** | Detection recall | TP / (TP + FN) |
| **JSC** | Jaccard similarity (tracks) | Overlap measure |
| **JSCθ** | Jaccard with temporal weighting | Time-weighted overlap |
| **RMSE** | Localization error | √(mean squared position error) |

---

## Part 3: Parameters to Modify for Benchmarking in run_demo.m

### Currently Available in generate6DTrajectories.m

```matlab
%% === PHOTOPHYSICS (Blinking/Bleaching) ===
'Signal', 1500, ...             % [photons] Higher = better SNR
'SignalStd', 300, ...           % [photons] Signal variation
'BlinkOffRate', 0.10, ...       % [prob/frame] 0.01-0.5
'BlinkOnRate', 0.70, ...        % [prob/frame] 0.1-0.95
'BleachRate', 0.0005, ...       % [prob/frame] Permanent photobleaching

%% === LOCALIZATION NOISE ===
'AddNoise', true, ...           % Enable localization noise
'NoiseXY', 20, ...              % [nm] Lateral localization precision
'NoiseZ', 50, ...               % [nm] Axial localization precision  
'NoiseTheta', 5, ...            % [degrees] Polar angle precision
'NoisePhi', 5, ...              % [degrees] Azimuthal angle precision
'NoiseGamma', 0.05, ...         % Wobble parameter noise

%% === MOTION PARAMETERS ===
'Dtrans', 500, ...              % [nm²/s] Translational diffusion
'Drot', 0.05, ...               % [rad²/s] Rotational diffusion

%% === DENSITY ===
'NumParticles', 5, ...          % Increase for higher density
'FieldOfView', [0 6000 0 6000], % [nm] Decrease for higher density
```

### Recommended Benchmark Test Matrix

#### Test 1: SNR Sweep (Localization Noise)
```matlab
% Low SNR (challenging)
'Signal', 500, 'NoiseXY', 50, 'NoiseZ', 100, 'NoiseTheta', 15, 'NoisePhi', 15

% Medium SNR  
'Signal', 1500, 'NoiseXY', 20, 'NoiseZ', 50, 'NoiseTheta', 5, 'NoisePhi', 5

% High SNR (easy)
'Signal', 5000, 'NoiseXY', 10, 'NoiseZ', 20, 'NoiseTheta', 2, 'NoisePhi', 2
```

#### Test 2: Particle Density Sweep
```matlab
% Low density
'NumParticles', 3, 'FieldOfView', [0 6000 0 6000]

% Medium density
'NumParticles', 10, 'FieldOfView', [0 6000 0 6000]

% High density
'NumParticles', 30, 'FieldOfView', [0 6000 0 6000]
```

#### Test 3: Blinking Severity Sweep
```matlab
% Minimal blinking (stable emitter)
'BlinkOffRate', 0.02, 'BlinkOnRate', 0.90

% Moderate blinking (default)
'BlinkOffRate', 0.10, 'BlinkOnRate', 0.70

% Heavy blinking (STORM-like)
'BlinkOffRate', 0.30, 'BlinkOnRate', 0.40
```

#### Test 4: Motion Speed Sweep
```matlab
% Slow diffusion
'Dtrans', 100, 'Drot', 0.01

% Medium diffusion (default)
'Dtrans', 500, 'Drot', 0.05

% Fast diffusion
'Dtrans', 2000, 'Drot', 0.20
```

### SNR Calculation for Simulated Data

The effective SNR in our simulation can be approximated as:

```
SNR ≈ Signal / sqrt(Signal + NoiseXY²)
```

For `Signal=1500, NoiseXY=20`: SNR ≈ 1500/√1520 ≈ 38 (very high)

To get ISBI-like SNR=4:
```matlab
'Signal', 200, 'NoiseXY', 50  % SNR ≈ 200/√2700 ≈ 4
```

---

## Part 4: Suggested Benchmark Scenarios for 6D-SMOLM

### Scenario Matrix for Publication

| Test ID | SNR | Density | Blinking | D_trans | D_rot | Purpose |
|---------|-----|---------|----------|---------|-------|---------|
| A1 | High | Low | Light | 500 | 0.05 | Baseline/validation |
| A2 | Med | Low | Light | 500 | 0.05 | SNR effect |
| A3 | Low | Low | Light | 500 | 0.05 | SNR limit |
| B1 | High | Med | Light | 500 | 0.05 | Density effect |
| B2 | High | High | Light | 500 | 0.05 | Density limit |
| C1 | High | Low | Heavy | 500 | 0.05 | Blinking effect |
| C2 | High | Low | Light | 2000 | 0.20 | Fast motion |
| D1 | Med | Med | Moderate | 500 | 0.05 | Realistic conditions |

### Implementation in run_demo.m

```matlab
% === BENCHMARK SCENARIO D1: Realistic conditions ===
scenarioName = 'benchmark_realistic';
simParams = {
    'NumParticles', 10, ...         % Medium density
    'NumFrames', 500, ...
    'Dtrans', 500, ...              % [nm²/s]
    'Drot', 0.05, ...               % [rad²/s]
    'SpatialModel', 'brownian', ...
    'RotationalModel', 'brownian', ...
    'GammaMean', 0.7, ...
    'Signal', 1000, ...             % Medium SNR
    'SignalStd', 200, ...
    'BlinkOffRate', 0.15, ...       % Moderate blinking
    'BlinkOnRate', 0.60, ...
    'AddNoise', true, ...
    'NoiseXY', 30, ...              % [nm] Realistic localization error
    'NoiseZ', 60, ...               % [nm]
    'NoiseTheta', 8, ...            % [deg] Realistic orientation error
    'NoisePhi', 8, ...              % [deg]
    'StartSpread', 500, ...
    'Seed', 42 ...
};
```

---

## References

1. Zondervan et al. "Photoblinking of Rhodamine 6G in Poly(vinyl alcohol)" J. Phys. Chem. A (2003)
2. Chenouard et al. "Objective comparison of particle tracking methods" Nature Methods (2014)
3. Roudot et al. "u-track3D: Measuring, navigating, and validating dense particle trajectories" Cell Reports Methods (2023)
4. ISBI Particle Tracking Challenge: http://bioimageanalysis.org/track/
