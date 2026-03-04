# software6D — 6D-SMOLM Tracking Pipeline

Extension of u-track3D for 6D single-molecule orientation and localization microscopy (6D-SMOLM).

## Overview

Adds orientation-aware cost functions and Landes group analysis insights to the u-track3D particle tracking framework, enabling simultaneous tracking of 3D positions (x, y, z) and 3D molecular orientations (θ, φ, Ω) of fluorescent molecules.

## Requirements

- MATLAB R2020a or later
- u-track3D (https://github.com/DanuserLab/u-track3D) on the MATLAB path
- Image Processing Toolbox, Statistics Toolbox, Optimization Toolbox

## Quick Start

```matlab
% 1. Add to path
addpath(genpath('path/to/u-track3D/software'))
addpath(genpath('path/to/software6D'))

% 2. Edit run6DTracking.m — set inputFile and outputDir
edit run6DTracking.m

% 3. Run
run6DTracking
```

## Testing (without u-track3D)

```matlab
cd software6D/tests
test_pipeline_synthetic
```

## Architecture

```
software6D/
├── 6D_tracking_pipeline_flowchart.html   ← Interactive pipeline diagram
├── run6DTracking.m                        ← Main wrapper script
├── README.md
├── preprocessing/
│   └── convert6DSMOLMtoMovieInfo.m       ← Fusion output → movieInfo
├── cost_functions/
│   ├── angularDistance.m                  ← Dipole distance with 180° symmetry
│   ├── costMat6DSMOLMLink.m              ← Cost Matrix 1: frame-to-frame
│   └── costMat6DSMOLMCloseGaps.m         ← Cost Matrix 2: gap closing (Bo's MLE)
├── postprocessing/
│   ├── extractTracks6D.m                 ← tracksFinal → trajectory matrix
│   ├── traj_filt_6D.m                    ← Filter + angular interpolation
│   └── combineTrajectories6D.m           ← Merge multi-movie trajectories
├── analysis/
│   └── estimateDiffusion6D_MLE.m         ← D_trans + D_rot with variable time lags
└── tests/
    └── test_pipeline_synthetic.m         ← Component tests with synthetic data
```

## Cost Matrix Structure

The 6D linking cost combines spatial and orientational distances:

```
cost_ij = w_pos * D_pos(i,j) + w_ori * D_ori(i,j)

where:
  D_pos    = sqrt(dx^2 + dy^2 + dz^2)             — 3D Euclidean distance
  D_ori    = w_angle * D_angle + w_wobble * D_wobble
  D_angle  = acos(|mu_i . mu_j|)                   — Dipole angle (0 to pi/2)
  D_wobble = |Omega_i - Omega_j|                   — Wobble cone similarity
  mu       = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]

Optional gamma-confidence weighting:
  gamma = 1 - 3*Omega/(4*pi) + Omega^2/(8*pi^2)
  D_ori *= (gamma_i + gamma_j) / 2
  (gamma ~ 1: fixed dipole, high confidence; gamma ~ 0: large wobble, reduce penalty)
```

## Key Design Decisions

| Decision | Value |
|----------|-------|
| Notation | θ ∈ [0, π/2], φ ∈ [-π, π], Ω ∈ [0, 2π] sr |
| Core tracker | `trackCloseGapsKalmanSparse` — unmodified |
| Kalman filter | Spatial only (x,y,z); orientation is a constraint |
| Dipole symmetry | 180° via `|cos(angle)|` in `angularDistance` |
| Bo's insight | `√(timeGap)` normalization in Cost Matrix 2 only |
| Motion models | Brownian spatial + Brownian rotational (independently selectable) |

## References

- Roudot, P. et al. Cell Reports Methods (2023) — u-track3D
- Jaqaman, K. et al. Nature Methods 5, 695–698 (2008) — u-track
- Shuang, B. et al. Langmuir 29, 228–234 (2013) — Photoblinking MLE

Emil Gillett, Landes Research Group, UIUC, 2026
