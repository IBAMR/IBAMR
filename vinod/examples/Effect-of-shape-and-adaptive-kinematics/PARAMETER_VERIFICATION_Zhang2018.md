# Parameter Verification: Zhang et al. (2018)
## "Effects of Reynolds number and thickness on an undulatory self-propelled foil"

**Paper**: Physics of Fluids 30, 071902 (2018)
**DOI**: https://doi.org/10.1063/1.5034439
**Authors**: Dong Zhang, Guang Pan, Liming Chao, Ya Zhang

---

## ✅ COMPLETE PARAMETER EXTRACTION FROM PAPER

### 1. Kinematics Model (Equation 1, Page 3)

**From Paper**:
```
y(x,t) = A_max × [(x + 0.03125)/(1.03125)] × sin[2π(x - t/T)]    (0 ≤ x ≤ 1)
```

**Where**:
- x = projected value of midline on x-axis (0 to 1)
- y = lateral displacement of midline
- A_max = **0.125 m** (maximum amplitude) ✓
- T = time period = 1/f
- f = frequency

**Amplitude Envelope**:
```
A(x) = A_max × (x + 0.03125)/(1.03125)
```

This is **LINEAR** from leading edge to trailing edge:
- At x=0 (leading edge): A(0) = 0.125 × 0.03125/1.03125 = 0.00379 m
- At x=1 (trailing edge): A(1) = 0.125 × 1.03125/1.03125 = 0.125 m

---

### 2. Reference Dimensions

| Parameter | Value | Our Implementation |
|-----------|-------|-------------------|
| Chord Length (L) | **1.0 m** | ✓ Assumed |
| Maximum Amplitude (A_max) | **0.125 m** | ✓ Correct (adaptive code) |
| Reference Velocity (U_Ref) | **1.0 L/s** | ✓ For nondimensionalization |

---

### 3. Reynolds Number Definition (Equation 2, Page 3)

```
Re = (V_max × L) / ν
```

**Where**:
- V_max = **2π × A_max × f** (maximum undulatory velocity)
- L = 1.0 m
- ν = kinematic viscosity (VARIED to keep Re constant as f varies)

**Example Calculation**:
- For f = 1.0 Hz, A_max = 0.125 m:
  - V_max = 2π × 0.125 × 1.0 = 0.785 m/s
  - For Re = 5609: ν = (0.785 × 1.0)/5609 = 1.4 × 10⁻⁴ m²/s ✓

---

### 4. Strouhal Number (Equation 3, Page 3)

```
St = (2 × A_max) / (T × U_0)
```

**Where**:
- U_0 = mean forward velocity of self-propelled foil (DEPENDENT variable)
- T = 1/f

**This is NOT a fixed value** - it varies with Re and frequency!

From Figure 5 (Page 5):
- St ranges from approximately 0.3 to 2.1 depending on Re
- St decreases with increasing Re
- At Re ≥ 10⁴, St approaches the "optimal" range of 0.2-0.4

---

### 5. Reynolds Numbers Tested (Page 3, Section III.A)

**From Paper**:
```
Re = 50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵
```

**Regimes**:
- **Viscous regime**: Re = 50
- **Intermediate regime**: Re = 500
- **Inertial regime**: Re ≥ 1000

**Our Adaptive Code**:
- Re = 1000 ✓ (matches paper)
- Re = 5609 ⚠️ (validation case, not in main study)
- Re = 10000 ✓ (matches paper)

**Missing from our implementation**:
- Re = 50, 500, 5000, 5×10⁴, 10⁵, 2×10⁵

---

### 6. Frequency Range (Page 3)

**From Paper**:
```
f varied from 0.3/s to 1.3/s with an interval of 0.1/s
```

**Tested frequencies**: 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3 Hz

**Our Implementation**:
- Base frequency = 0.785 Hz (which is π/4)
- ⚠️ We don't sweep frequency like the paper does

---

### 7. NACA Profiles and Thickness Study

#### For Reynolds Number Study (Section III.A):
**Profile**: NACA 0012 (fixed)

#### For Thickness Study (Section III.B):
**Profiles**: NACA 0004, 0008, 0012, 0016, 0020, 0024

**Thickness Ratios**: 0.04, 0.08, 0.12, 0.16, 0.20, 0.24

**Tested at Two Re Values**:
1. Re = **5×10⁴** (inertial regime)
2. Re = **500** (intermediate regime)

**Our Implementation**:
- Thickness ratios: 0.04, 0.06, 0.08 ⚠️ (incomplete range)
- Missing: 0.12, 0.16, 0.20, 0.24

---

### 8. Efficiency Definition (Equations 7-10, Page 4)

**From Paper**:
```
η = E / W

E = (1/2) × m × U_0²     (kinetic energy in forward direction)
m = ρ × V_f              (mass of foil)
W = ∫∫ -σ·n·U dS dt      (work performed by foil)
```

**Where**:
- ρ = fluid density = foil density (assumed equal)
- V_f = volume of foil
- σ = viscous stress tensor
- n = outward surface normal
- U = velocity of body surface

**Our Implementation**:
- ⚠️ Need to verify if we use this exact definition

---

### 9. Validation Case (Page 4, Section II.D)

**Parameters**:
- L = 1 m
- A_max = 0.125 m
- T = 1 s (f = 1 Hz)
- ν = **1.4 × 10⁻⁴ m²/s**
- **Re = 5609** ✓

**Body Width** (Equation 11):
```
w(x) = { sqrt(2*W_H*x - x²)           for 0 ≤ x < x_H
       { W_H * (L-x)/(L-x_H)          for x_H ≤ x < L
```
Where: W_H = x_H = **0.04 L**

**Results**:
- Mean forward velocity: U_0 = **0.61 L/s** ✓
- Mean lateral velocity: V_0 ≈ **0.02 L/s**
- Forward velocity oscillation amplitude: ∆U_0 = **10%** of U_0
- Lateral velocity oscillation amplitude: ∆V_0 = **9%** of U_0

---

### 10. Grid and Numerical Parameters (Page 4)

**Domain**: 10L × 4L (periodic)

**Grid Levels**: 4 levels with AMR
- Grid spacing ratio: 4:1 between levels
- Final grid: 40×16 cells at coarsest level = 2560×1024 equivalent

**Time Step**:
- Dynamic time step
- CFL number ≤ 0.1

---

## 🔍 COMPARISON WITH OUR IMPLEMENTATION

### Adaptive Research Code (`/home/user/Effect-of-shape-and-adaptive-kinematics/`)

| Parameter | Zhang et al. (2018) | Our Implementation | Status |
|-----------|---------------------|-------------------|--------|
| **Kinematics** | y(x,t) = A_max×[(x+0.03125)/1.03125]×sin[2π(x-t/T)] | Custom adaptive | ⚠️ DIFFERENT |
| **A_max** | 0.125 m | 0.125 m | ✅ CORRECT |
| **Length L** | 1.0 m | 1.0 m | ✅ CORRECT |
| **Frequency range** | 0.3-1.3 Hz (11 values) | Base: 0.785 Hz | ⚠️ INCOMPLETE |
| **Reynolds numbers** | 50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵ | 1000, 5609, 10000 | ⚠️ INCOMPLETE |
| **NACA profiles** | 0004, 0008, 0012, 0016, 0020, 0024 | Variable | ⚠️ INCOMPLETE |
| **Thickness ratios** | 0.04-0.24 | 0.04, 0.06, 0.08 | ⚠️ INCOMPLETE |
| **Strouhal number** | Computed (0.3-2.1 range) | Computed | ✅ CORRECT approach |
| **Efficiency** | η = E/W | Similar | ✅ LIKELY CORRECT |

---

### Validation Suite (`/home/user/Validation-Gupta2022/`)

| Parameter | Zhang et al. (2018) | Validation Suite | Status |
|-----------|---------------------|------------------|--------|
| **Paper basis** | Zhang et al. 2018 | **Gupta et al. 2022** | ℹ️ DIFFERENT PAPER |
| **Re** | Variable (50-2×10⁵) | Fixed: 5000 | ℹ️ Different study |
| **St** | Computed (variable) | Fixed: 0.6 | ℹ️ Different study |
| **Kinematics** | Linear envelope | Exponential/Quadratic | ℹ️ Different study |

**Conclusion**: Validation suite is correctly based on **different paper** (Gupta et al. 2022)

---

## ⚠️ KEY FINDINGS AND DISCREPANCIES

### 1. **Kinematics Equation** - CRITICAL DIFFERENCE ⚠️

**Zhang et al. (2018)**:
```
y(x,t) = A_max × [(x + 0.03125)/(1.03125)] × sin[2π(x - t/T)]
```
This is a **simple linear amplitude envelope**.

**Our Adaptive Code**:
Uses Reynolds-dependent adaptive kinematics with:
- Amplitude adaptation: A_adapted = A_base × f_Re(Re) × f_h(h/L)
- Frequency adaptation: f_adapted = f_base × g_Re(Re) × g_h(h/L)
- Envelope power: A(s) = A_max × (s/L)^p

**Status**: ⚠️ **OUR CODE GOES BEYOND THE PAPER** - we add adaptive features not in Zhang et al. (2018)

---

### 2. **Reynolds Number Range** - INCOMPLETE ⚠️

**Paper tests**: 50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵

**Our tests**: 1000, 5609, 10000

**Missing**: 50, 500, 5000, 5×10⁴, 10⁵, 2×10⁵

**Recommendation**: Add missing Re values to match paper exactly

---

### 3. **Frequency Sweep** - MISSING ⚠️

**Paper**: Sweeps f from 0.3 to 1.3 Hz (11 values)

**Our code**: Uses single base frequency 0.785 Hz

**Impact**: Cannot reproduce Figure 6(a,b,c,d) from paper without frequency sweep

---

### 4. **Thickness Range** - INCOMPLETE ⚠️

**Paper**: 0.04, 0.08, 0.12, 0.16, 0.20, 0.24

**Our code**: 0.04, 0.06, 0.08

**Missing**: 0.12, 0.16, 0.20, 0.24

**Impact**: Cannot reproduce thickness study results (Section III.B)

---

### 5. **Two Separate Studies in Zhang et al. (2018)**

**Study 1: Reynolds Number Effect** (Section III.A)
- Fixed: NACA 0012, A_max = 0.125 m
- Variable: Re (8 values), f (11 values)
- Results: Figures 6, 7

**Study 2: Thickness Effect** (Section III.B)
- Fixed: Re = 5×10⁴ or 500, A_max = 0.125 m, f sweep
- Variable: NACA profiles (6 types)
- Results: Figures 8-13

**Our code**: Combines both into one "adaptive" framework ⚠️

---

## 📊 KEY RESULTS FROM PAPER (For Verification)

### From Page 2 (Abstract):

> "thinner foils outperformed thicker foils in terms of the forward velocity and input power at both Re values"

> "However, the efficiency related to the conversion of input power into kinetic energy for NACA 0004 was the lowest."

> "An optimum thickness exists that depends on Re."

### From Page 6 (Section III.A.1):

**Re Effect on Velocity** (Figure 6a):
- At f=1.0 Hz:
  - Re=50: U_0 ≈ 0.15 L/s
  - Re=500: U_0 ≈ 0.30 L/s
  - Re=1000: U_0 ≈ 0.45 L/s
  - Re=5000: U_0 ≈ 0.65 L/s
  - Re=10⁴: U_0 ≈ 0.72 L/s
  - Re=5×10⁴: U_0 ≈ 0.80 L/s
  - Re=10⁵: U_0 ≈ 0.82 L/s
  - Re=2×10⁵: U_0 ≈ 0.82 L/s (asymptotic)

**Re Effect Asymptotic Beyond Re = 5×10⁴** ✓

### From Page 7 (Section III.A.2):

**Four Vortex Types**:
1. **Re=50** (viscous): Leading-edge vortex fails to shed
2. **Re=500** (intermediate): 2S reverse von Kármán
3. **Re=1000-10⁴**: 2S → 2P transition
4. **Re ≥ 5×10⁴**: 2P + Kelvin-Helmholtz instability

### From Page 8-9 (Section III.B):

**Optimal Thickness** (Figure 9e, 12e):
- At Re = 5×10⁴: **NACA 0016** (η ≈ 1.4)
- At Re = 500: **NACA 0016** (η ≈ 0.12)

**NACA 0004**:
- Highest velocity ✓
- Lowest power ✓
- **Lowest efficiency** ✓ (due to low mass → low kinetic energy)

---

## 📝 WHAT NEEDS TO BE UPDATED

### 1. Documentation Updates

#### README.md (Adaptive Code)
**Current**:
```
This code implements the research described in:
- "Effects of Reynolds number and thickness on an undulatory self-propelled foil"
- "Anguilliform and carangiform fish-inspired hydrodynamic study..."
```

**Should Be**:
```
This code is INSPIRED BY and EXTENDS:

Primary Reference:
[1] Zhang, D., Pan, G., Chao, L., & Zhang, Y. (2018).
    "Effects of Reynolds number and thickness on an undulatory self-propelled foil"
    Physics of Fluids, 30, 071902.
    DOI: https://doi.org/10.1063/1.5034439

Extensions in This Code:
- Reynolds-dependent adaptive kinematics (not in original paper)
- Adaptive amplitude and frequency modulation
- Extended parameter ranges
```

#### Add Citation
**BibTeX**:
```bibtex
@article{Zhang2018,
  title={Effects of Reynolds number and thickness on an undulatory self-propelled foil},
  author={Zhang, Dong and Pan, Guang and Chao, Liming and Zhang, Ya},
  journal={Physics of Fluids},
  volume={30},
  number={7},
  pages={071902},
  year={2018},
  publisher={AIP Publishing},
  doi={10.1063/1.5034439}
}
```

---

### 2. Parameter Alignment

**To exactly reproduce Zhang et al. (2018)**, create:

#### `input2d_Zhang2018_Re50_NACA0012`
- Re = 50
- NACA 0012
- f = 1.0 Hz
- A_max = 0.125 m
- Kinematics: y(x,t) = 0.125 × [(x+0.03125)/1.03125] × sin[2π(x-t)]

#### `input2d_Zhang2018_Re500_NACA0012`
- Re = 500
- (similar...)

#### Continue for all Re values...

---

### 3. Add Frequency Sweep Script

```bash
#!/bin/bash
# reproduce_zhang2018_fig6.sh

for RE in 50 500 1000 5000 10000 50000 100000 200000; do
    for FREQ in 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3; do
        echo "Running Re=$RE, f=$FREQ"
        # Generate input file with these parameters
        # Run simulation
    done
done
```

---

### 4. Add Thickness Study Cases

For Re = 5×10⁴:
- input2d_Zhang2018_Re50000_NACA0004
- input2d_Zhang2018_Re50000_NACA0008
- input2d_Zhang2018_Re50000_NACA0012
- input2d_Zhang2018_Re50000_NACA0016
- input2d_Zhang2018_Re50000_NACA0020
- input2d_Zhang2018_Re50000_NACA0024

For Re = 500:
- (same set with Re=500)

---

## ✅ CORRECTED PARAMETER TABLE

### Zhang et al. (2018) - Exact Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| **Chord Length** | L = 1.0 m | Fixed |
| **Max Amplitude** | A_max = 0.125 m | Fixed |
| **Amplitude Envelope** | A(x) = 0.125×(x+0.03125)/1.03125 | Linear |
| **Kinematics** | y(x,t) = A(x)×sin[2π(x-t/T)] | Equation 1 |
| **Frequency Range** | f = 0.3 to 1.3 Hz (∆f=0.1) | 11 values |
| **Period** | T = 1/f | Variable |
| **Reynolds Numbers** | 50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵ | 8 values |
| **V_max** | 2π × A_max × f | Max undulatory velocity |
| **Re Definition** | Re = V_max × L / ν | Equation 2 |
| **Kinematic Viscosity** | ν = V_max × L / Re | Adjusted per Re |
| **Strouhal Number** | St = 2×A_max/(T×U_0) | Computed output |
| **NACA Profiles** | 0012 (Re study), 0004-0024 (thickness study) | 6 for thickness |
| **Thickness Ratios** | 0.04, 0.08, 0.12, 0.16, 0.20, 0.24 | 6 values |
| **Efficiency** | η = E/W = (½mU_0²)/W | Takes mass into account |
| **Reference Velocity** | U_Ref = 1 L/s | For nondimensionalization |
| **Fluid/Foil Density** | ρ_fluid = ρ_foil | Assumed equal |

---

## 🎯 VALIDATION CHECKLIST

To verify our implementation matches Zhang et al. (2018):

### Re Study (NACA 0012, A_max=0.125m)
- [ ] Re=50, f=1.0Hz → U_0 ≈ 0.15 L/s
- [ ] Re=500, f=1.0Hz → U_0 ≈ 0.30 L/s
- [ ] Re=1000, f=1.0Hz → U_0 ≈ 0.45 L/s
- [ ] Re=5000, f=1.0Hz → U_0 ≈ 0.65 L/s
- [ ] Re=10⁴, f=1.0Hz → U_0 ≈ 0.72 L/s
- [ ] Re=5×10⁴, f=1.0Hz → U_0 ≈ 0.80 L/s
- [ ] Re=10⁵, f=1.0Hz → U_0 ≈ 0.82 L/s
- [ ] Re=2×10⁵, f=1.0Hz → U_0 ≈ 0.82 L/s (asymptotic)

### Thickness Study (Re=5×10⁴, f=1.0Hz)
- [ ] NACA 0004 → Highest U_0, lowest η
- [ ] NACA 0016 → Optimal η ≈ 1.4
- [ ] NACA 0024 → Lowest U_0, moderate η

### Wake Structures
- [ ] Re=50 → Leading-edge vortex
- [ ] Re=500 → 2S wake
- [ ] Re=10⁴ → 2P wake
- [ ] Re≥5×10⁴ → 2P + K-H instability

---

## 📚 REFERENCE INFORMATION

**Full Citation**:
```
Zhang, D., Pan, G., Chao, L., & Zhang, Y. (2018). Effects of Reynolds number
and thickness on an undulatory self-propelled foil. Physics of Fluids, 30(7), 071902.
https://doi.org/10.1063/1.5034439
```

**Repository Location**:
```
https://github.com/vinodthale/Effect-of-shape-and-adaptive-kinematics/blob/main/
Effects%20of%20Reynolds%20number%20and%20thickness%20on%20an%20undulatory%20self%20propelled%20foil.pdf
```

**Paper Location (Local)**:
```
/home/user/Effect-of-shape-and-adaptive-kinematics/
Effects of Reynolds number and thickness on an undulatory self propelled foil.pdf
```

---

**Verification Date**: 2025-11-17
**Status**: Parameters extracted and compared
**Next Steps**: Update documentation, add missing cases, verify numerical results

---
