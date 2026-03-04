# Orientation Cost Implementation in 6D-SMOLM Tracking

## Overview

This document explains how orientation costs are incorporated into the 
6D-SMOLM tracking pipeline. Two modes are available:

1. **DIRECT MODE** (default): Simple frame-to-frame angular distance
2. **PSEUDO-KALMAN MODE** (new): Uncertainty-weighted Mahalanobis-like costs

---

## Mode 1: Direct Comparison (Default)

### Frame-to-Frame Linking (`costMat6DSMOLMLink.m`)

For linking detections between frame `t` and frame `t+1`:

```
Total Cost = wSpatial × C_spatial + wOrient × C_orient + wOmega × C_omega
```

Where:

| Component | Source | Description |
|-----------|--------|-------------|
| `C_spatial` | u-track3D Kalman filter | Mahalanobis distance from **predicted** position |
| `C_orient` | Direct comparison | `α² / maxAngDist²` (α = angular displacement) |
| `C_omega` | Direct comparison | `Δω² / (2π)²` (Δω = wobble difference, range [0, 2π] sr) |

### Limitations of Direct Mode

| | Spatial Cost | Orientation Cost |
|---|---|---|
| **Uses Kalman prediction?** | ✅ YES | ❌ NO |
| **Accounts for motion model?** | ✅ YES (Brownian/directed) | ❌ NO |
| **Uncertainty-weighted?** | ✅ YES (Mahalanobis) | ❌ NO (raw distance) |

---

## Mode 2: Pseudo-Kalman (NEW!)

### Enabling Pseudo-Kalman Mode

In `run_demo.m`:
```matlab
kalmanParams.useOrientKalman = true;       % Enable pseudo-Kalman
kalmanParams.D_rot_prior = 0.05;           % rad²/s (rotational diffusion)
kalmanParams.orientCRLB = 5 * pi/180;      % ~5 degrees measurement uncertainty
kalmanParams.frameTime = 0.05;             % 50 ms frame time
kalmanParams.useSignalWeighting = false;   % Optional signal-dependent CRLB
```

### The Pseudo-Kalman Formula

```
C_orient = α² / (P_predicted + R_meas)
```

Where:
- `α` = angular displacement between orientations (radians)
- `P_predicted = 2 × D_rot × Δt` = expected angular variance from diffusion
- `R_meas = σ_CRLB²` = measurement variance

This is analogous to the Kalman innovation cost:
```
(z - H*x_pred)ᵀ × S⁻¹ × (z - H*x_pred)
```
Where `S = P_predicted + R_meas` is the innovation covariance.

### Physical Interpretation

| Parameter | Meaning | Typical Value |
|-----------|---------|---------------|
| `D_rot` | Rotational diffusion coefficient | 0.01-0.1 rad²/s |
| `Δt` | Frame time | 0.05 s (50 ms) |
| `P_predicted` | Expected angular spread | 2×0.05×0.05 = 0.005 rad² |
| `σ_CRLB` | Measurement uncertainty | 5° ≈ 0.087 rad |
| `R_meas` | Measurement variance | 0.087² ≈ 0.0076 rad² |

### Benefits of Pseudo-Kalman Mode

1. **Uncertainty-weighted**: 
   - Large angular distances are penalized relative to expected variance
   - If `D_rot` is high, large angular changes are more acceptable
   - If `orientCRLB` is high, measurements are trusted less

2. **Time-aware**:
   - `P_predicted` grows with `dt`, naturally handling different frame rates
   - No need for ad-hoc Bo scaling

3. **Signal-weighted** (optional):
   - When `useSignalWeighting = true`, bright molecules have lower `R_meas`
   - `R_meas(signal) = orientCRLB² × (refSignal / signal)`
   - Bright molecules contribute more reliably to linking

### Comparison: Direct vs Pseudo-Kalman

| Scenario | Direct Mode | Pseudo-Kalman |
|----------|-------------|---------------|
| Small angular change (0.1 rad) | 0.1²/(π/2)² = 0.004 | 0.1²/(0.005+0.008) = 0.77 |
| Large angular change (0.5 rad) | 0.5²/(π/2)² = 0.10 | 0.5²/(0.005+0.008) = 19.2 |
| **Ratio** | 25× | 25× |

The ratio is similar, but pseudo-Kalman provides:
- Physically meaningful normalization
- Automatic scaling with D_rot and measurement quality

---

## Angular Distance Calculation

Both modes use the same angular distance calculation:

### Variable Naming Convention

To avoid confusion with the polar angle θ (theta), we use **ψ (psi)** for angular displacement:
- **θ, φ** = polar and azimuthal angles defining dipole orientation
- **ψ** = angular displacement between two dipole orientations

### Convert (θ, φ) to Unit Vectors
```
μ̂ = [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
```

### Compute Angular Displacement with Head-Tail Symmetry
```
ψ = arccos(|μ̂₁ · μ̂₂|)    ∈ [0, π/2]
```

The absolute value `|·|` enforces that μ̂ and -μ̂ are equivalent (dipole symmetry).

### Rotational MSD

The expected mean squared angular displacement follows:
```
⟨Δψ²⟩ = 2 D_rot Δt + 2 σ_ψ²
```
Where:
- `D_rot` = rotational diffusion coefficient (rad²/s)
- `Δt` = time lag
- `σ_ψ` = orientation measurement precision (from CRLB)

---

## Gap Closing

For gap closing across multiple frames, the same principles apply with:

### Bo Shuang Scaling (Direct Mode)
```
d_angle_scaled = d_angle / √(timeGap)
```

### Time-Gap Scaling (Pseudo-Kalman Mode)
```
P_predicted = 2 × D_rot × (timeGap × dt)
```

The variance grows linearly with time gap, naturally accounting for expected diffusion.

---

## Default Parameter Values

```matlab
% Weights (both modes)
wSpatial = 1.0;        % Full spatial cost
wOrient  = 0.3;        % 30% weight on orientation
wOmega   = 0.5;        % 50% weight on wobble

% Direct mode
maxAngDist = π/2;      % 90° normalization

% Pseudo-Kalman mode
D_rot_prior = 0.05;    % rad²/s
orientCRLB = 5°;       % ~0.087 rad
frameTime = 0.05;      % 50 ms
```

---

## Future Enhancements

The pseudo-Kalman approach is a stepping stone to full orientation Kalman filtering:

### Full Kalman Filter (not yet implemented)
- State: `[θ, φ, ω_θ, ω_φ]` (angles + angular velocities)
- Quaternion representation for singularity-free tracking
- Track-specific D_rot estimation
- Innovation-based outlier rejection

### Current Limitations
- No angular velocity tracking (assumes D_rot dominates)
- Single D_rot_prior for all particles
- No track history used (stateless)

---

*Emil Gillett, Landes Research Group, UIUC, 2026*
