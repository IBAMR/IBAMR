# Zhang et al. (2018) Complete Test Matrix

This document defines the complete parameter space from Zhang et al. (2018).

---

## Reynolds Number Cases

Zhang tested 8 Reynolds numbers:

| Re | Regime | Characteristics | Input File Pattern |
|----|--------|-----------------|-------------------|
| 50 | Viscous-dominated | Weak vortices, diffusion | `input2d_Zhang_Re50_hXXX` |
| 500 | Low Re transitional | 2S wake pattern | `input2d_Zhang_Re500_hXXX` |
| 1000 | Transitional | 2S or 2P wake | `input2d_Zhang_Re1000_hXXX` |
| 5000 | Moderate Re | 2P wake pattern | `input2d_Zhang_Re5000_hXXX` |
| 10000 | Inertial regime | 2P pattern, stronger thrust | `input2d_Zhang_Re10000_hXXX` |
| 50000 | High Re | 2P + KH eddies | `input2d_Zhang_Re50000_hXXX` |
| 100000 | Very high Re | Multiple vortices, turbulent | `input2d_Zhang_Re100000_hXXX` |
| 200000 | Maximum Re | Complex wake, transition | `input2d_Zhang_Re200000_hXXX` |

---

## Thickness Ratio Cases

Zhang tested 6 thickness ratios (NACA 00XX foils):

| h/L | NACA Foil | Volume V_f* | Input File Pattern |
|-----|-----------|-------------|-------------------|
| 0.04 | NACA0004 | 0.03146 | `input2d_Zhang_ReXXXX_h004` |
| 0.08 | NACA0008 | 0.05887 | `input2d_Zhang_ReXXXX_h008` |
| 0.12 | NACA0012 | 0.08621 | `input2d_Zhang_ReXXXX_h012` |
| 0.16 | NACA0016 | 0.1135 | `input2d_Zhang_ReXXXX_h016` |
| 0.20 | NACA0020 | 0.1408 | `input2d_Zhang_ReXXXX_h020` |
| 0.24 | NACA0024 | 0.1681 | `input2d_Zhang_ReXXXX_h024` |

---

## Complete Parameter Matrix

**Total test cases**: 8 Re × 6 h/L = **48 simulations**

### Baseline Cases (h/L = 0.04, NACA0004)

- [x] `input2d_Zhang_Re50_h004` (created)
- [x] `input2d_Zhang_Re500_h004` (to be created)
- [x] `input2d_Zhang_Re1000_h004` (created)
- [ ] `input2d_Zhang_Re5000_h004`
- [ ] `input2d_Zhang_Re10000_h004`
- [ ] `input2d_Zhang_Re50000_h004`
- [ ] `input2d_Zhang_Re100000_h004`
- [ ] `input2d_Zhang_Re200000_h004`

### All Thickness Cases (Full Matrix)

For each Re above, create variants for:
- h004 (NACA0004)
- h008 (NACA0008)
- h012 (NACA0012)
- h016 (NACA0016)
- h020 (NACA0020)
- h024 (NACA0024)

---

## Fixed Parameters (All Cases)

These parameters are **constant** across all 48 test cases:

```
base_amplitude = 0.125
base_frequency = 0.785
envelope_c0 = 0.03125
envelope_c1 = 1.03125
envelope_power = 1.0
swimming_mode = 0.0
enable_shape_adaptation = FALSE
```

---

## Variable Parameters

### Per-Case Reynolds Number
```
Re = {50, 500, 1000, 5000, 10000, 50000, 100000, 200000}
MU = 0.785 / Re
```

### Per-Case Thickness
```
thickness_ratio = {0.04, 0.08, 0.12, 0.16, 0.20, 0.24}
```

### Per-Case Simulation Settings

**Low Re (50-1000):**
```
N = 128 (finer grid)
DT_MAX = 0.00005
CFL_MAX = 0.2
END_TIME = 20.0 (longer simulation)
```

**Moderate Re (5000-10000):**
```
N = 64
DT_MAX = 0.0001
CFL_MAX = 0.3
END_TIME = 10.0
```

**High Re (50000-200000):**
```
N = 64 (or 128 if needed)
DT_MAX = 0.00005
CFL_MAX = 0.25
END_TIME = 8.0
```

---

## Expected Outputs (from Zhang Fig. 4-6)

### Swimming Speed U₀

| Re | h/L=0.04 | h/L=0.08 | h/L=0.12 | h/L=0.16 | h/L=0.20 | h/L=0.24 |
|----|----------|----------|----------|----------|----------|----------|
| 50 | ~0.05 | ~0.03 | ~0.02 | ~0.01 | ~0.01 | ~0.01 |
| 1000 | ~0.18 | ~0.15 | ~0.12 | ~0.10 | ~0.08 | ~0.07 |
| 10000 | ~0.35 | ~0.32 | ~0.28 | ~0.25 | ~0.22 | ~0.20 |
| 100000 | ~0.50 | ~0.48 | ~0.45 | ~0.42 | ~0.38 | ~0.35 |

*(Values are approximate, read from Zhang Fig. 4)*

### Strouhal Number St

Peaks around Re ~ 1000-5000, range [0.3, 2.1] overall.

### Efficiency η

Increases with Re, optimal thickness depends on Re.

---

## Generation Script

To generate all 48 input files, use:

```bash
cd Zhang_2018/input_files
./generate_all_zhang_cases.sh
```

Or manually create from templates following the pattern in:
- `input2d_Zhang_Re50_h004`
- `input2d_Zhang_Re1000_h004`

---

## Validation Criteria

For each simulation, verify:
- [ ] U₀ matches Zhang Fig. 4 trend
- [ ] St matches Zhang Fig. 5 range
- [ ] η matches Zhang Fig. 6 trend
- [ ] Wake pattern matches regime:
  - Re 50: weak, diffuse
  - Re 500-1000: 2S
  - Re 5000-10000: 2P
  - Re ≥ 50000: 2P + KH eddies

---

## Priority Cases for Initial Validation

Recommend running these first:

1. **Re = 1000, h/L = 0.04** - Baseline transitional regime
2. **Re = 5000, h/L = 0.04** - Moderate Re, well-documented
3. **Re = 10000, h/L = 0.08** - Higher Re, thicker foil
4. **Re = 100000, h/L = 0.04** - High Re validation

Then expand to full matrix.

---

## Computational Requirements

**Per simulation estimate:**
- Low Re (50-1000): 12-24 hours on 6-12 cores
- Moderate Re (5000-10000): 6-12 hours on 6-12 cores
- High Re (≥50000): 12-48 hours on 12-24 cores (may need more resolution)

**Total for all 48 cases:**
- Approximately 600-1200 core-hours
- Can be parallelized (run multiple cases simultaneously)

---

## File Naming Convention

```
input2d_Zhang_Re[VALUE]_h[XXX]

Examples:
input2d_Zhang_Re50_h004       → Re=50, h/L=0.04
input2d_Zhang_Re1000_h012     → Re=1000, h/L=0.12
input2d_Zhang_Re100000_h024   → Re=100000, h/L=0.24
```

---

**Document Version**: 1.0
**Last Updated**: 2025-01-17
**Status**: Test Matrix Definition
