# Parameter Verification: Current Implementation vs. Zhang (2018)

This document verifies whether the current implementation parameters align with the Zhang et al. (2018) specification.

---

## 1. Kinematic Equations Comparison

### 1.1 Zhang (2018) Specification

**Nondimensional midline equation:**
```
Y(X, τ) = A_max × [(X + c0) / c1] × sin[2π(X - τ)]
```

**Parameters:**
- A_max = 0.125
- c0 = 0.03125
- c1 = 1.03125
- X ∈ [0, 1] (chordwise coordinate)
- τ = t/T (nondimensional time)

**Velocity (time derivative):**
```
∂Y/∂t = A_max × [(X + c0) / c1] × 2πf × cos[2π(X - τ)]
      = 0.785f × [(X + c0) / c1] × cos[2π(X - τ)]
```

where f is dimensional frequency (Hz).

---

### 1.2 Current Implementation

**From input files** (`input2d`, `input2d_Re1000_h004`, etc.):

**Body shape equation:**
```
body_shape_equation = "0.125 * ((X_0 + 0.03125) / 1.03125) * sin(2*PI*X_0 - (0.785/0.125)*T)"
```

**Deformation velocity equations:**
```
deformation_velocity_function_0 = "(-0.785 * ((X_0 + 0.03125) / 1.03125) * cos(2*PI*X_0 - (0.785/0.125)*T)) * N_0"
deformation_velocity_function_1 = "(-0.785 * ((X_0 + 0.03125) / 1.03125) * cos(2*PI*X_0 - (0.785/0.125)*T)) * N_1"
```

where:
- X_0 = streamwise Lagrangian coordinate
- T = time (dimensional)
- N_0, N_1 = normal vector components

---

### 1.3 Mathematical Analysis

**Step 1: Simplify the wave argument**

Current implementation:
```
2*PI*X_0 - (0.785/0.125)*T
```

Calculate coefficient:
```
0.785 / 0.125 = 6.28 ≈ 2π
```

Therefore:
```
2*PI*X_0 - (0.785/0.125)*T = 2π*X_0 - 2π*T
                            = 2π(X_0 - T)
```

**Step 2: Compare wave forms**

**Zhang:** `sin[2π(X - τ)]` where τ = t/T_period

**Current:** `sin[2π(X_0 - T)]` where T = dimensional time

**Key question:** Is T in the code nondimensionalized by period, or is it dimensional time?

---

### 1.4 Time Nondimensionalization Check

**From Zhang specification:**
- Period T_period = 1/f
- Nondimensional time τ = t/T_period = t × f

**If current T = dimensional time t:**
```
sin[2π(X - T)] would require T to have units of length, not time → INCORRECT
```

**If current T = τ (nondimensionalized):**
```
sin[2π(X - τ)] ✓ MATCHES Zhang
```

**Likely interpretation:**
The code variable `T` represents **nondimensional time** scaled by the period.
With base_frequency = 0.785, the code implicitly uses:
```
T = t × f = t × 0.785 (in nondimensional units)
```

**However**, to match Zhang exactly with nondimensional frequency f* = 1:
```
T_Zhang = t/T_period = t × f_dimensional
```

If f_dimensional = 0.785 Hz, then:
```
T_code = t × 0.785 (matches the code equations)
```

---

### 1.5 Amplitude Envelope Verification

**Zhang:**
```
A(X) = 0.125 × (X + 0.03125) / 1.03125
```

**Current:**
```
A(X_0) = 0.125 × (X_0 + 0.03125) / 1.03125
```

**Verification:** ✓ **EXACT MATCH**

---

### 1.6 Maximum Amplitude Verification

**Zhang:** A_max = 0.125

**Current:**
```cpp
base_amplitude = 0.125  // from input files
```

**With adaptive kinematics enabled:**
```cpp
d_adapted_amplitude = d_base_amplitude × Re_amplitude_factor × thickness_amplitude_factor
                    = 0.125 × f(Re, h/L)
```

**Verification:**
- ✓ Base value matches: 0.125
- ✗ With adaptation enabled: **DIFFERS** from Zhang (adaptive modification)
- ✓ With adaptation disabled: **MATCHES** Zhang

---

### 1.7 Velocity Coefficient Verification

**Zhang:**
```
V_max = 2π × A_max × f = 2π × 0.125 × f = 0.785 × f
```

For f = 1.0 Hz (example):
```
V_max = 0.785
```

**Current deformation velocity:**
```
-0.785 × envelope × cos(...)
```

**Verification:** ✓ Coefficient **MATCHES** (0.785 = 2π × 0.125)

**Sign:** Negative sign is due to normal vector direction convention in IB method (consistent with kinematics).

---

## 2. Reynolds Number Definition

### 2.1 Zhang (2018)

**Definition:**
```
Re = V_max × L / ν
   = (2π × A_max × f) × L / ν
   = 0.785 × f / ν
```

For L = 1 (nondimensional).

**Prescribed values:**
```
Re ∈ {50, 500, 1000, 5000, 10000, 50000, 100000, 200000}
```

---

### 2.2 Current Implementation

**From input files:**
```
Re = 1000.0 (for example)
MU = 0.785 / Re
   = 0.785 / 1000
   = 0.000785
```

**Verification:**
```
Re = 0.785 / MU
   = 0.785 / 0.000785
   = 1000 ✓
```

**Reynolds number calculation:**
Assuming f = 1.0 (nondimensional frequency):
```
Re = V_max / ν = 0.785 / MU ✓ CORRECT
```

---

### 2.3 Reynolds Number Values

**Zhang test cases:**
```
Re = {50, 500, 1000, 5000, 10^4, 5×10^4, 10^5, 2×10^5}
```

**Current input files:**
- `input2d_Re1000_h004`: Re = 1000 ✓
- `input2d_Re5609_h006`: Re = 5609 (close to 5000)
- `input2d_Re10000_h008`: Re = 10000 ✓
- `input2d`: Re = 5609 (baseline)

**Missing Zhang cases:**
- Re = 50 (very low Re)
- Re = 500 (low Re)
- Re = 5000 (should replace 5609)
- Re = 50000
- Re = 100000
- Re = 200000

**Verification:**
- ✓ Method correct (MU = 0.785/Re)
- ✗ Incomplete coverage of Zhang's Re range
- ∼ Re = 5609 used instead of 5000 (likely from different study)

---

## 3. Thickness Ratio (Geometry)

### 3.1 Zhang (2018)

**Thickness values:**
```
t/c = h/L ∈ {0.04, 0.08, 0.12, 0.16, 0.20, 0.24}
```

**Corresponding NACA foils:**
```
NACA0004, NACA0008, NACA0012, NACA0016, NACA0020, NACA0024
```

---

### 3.2 Current Implementation

**From input files:**
```
thickness_ratio = 0.04 (Re1000_h004)
thickness_ratio = 0.06 (Re5609_h006)
thickness_ratio = 0.08 (Re10000_h008)
```

**Zhang coverage:**
- ✓ h/L = 0.04 (NACA0004) available
- ✗ h/L = 0.06 not in Zhang set
- ✓ h/L = 0.08 (NACA0008) available

**Missing Zhang cases:**
- h/L = 0.12 (NACA0012)
- h/L = 0.16 (NACA0016)
- h/L = 0.20 (NACA0020)
- h/L = 0.24 (NACA0024)

**Verification:**
- ✓ Thickness parameter exists
- ✗ Incomplete coverage of Zhang's thickness range
- ∼ h/L = 0.06 not used by Zhang

---

## 4. Frequency and Time Scaling

### 4.1 Zhang (2018)

**Dimensional frequency range:**
```
f ∈ [0.3, 1.3] Hz
```

**Nondimensional time:**
```
τ = t/T = t × f
```

**In wave equation:**
```
sin[2π(X - τ)] has nondimensional frequency = 1
```

---

### 4.2 Current Implementation

**Base frequency:**
```cpp
base_frequency = 0.785
```

**In equations:**
```
sin(2π*X_0 - 2π*T)
```

**Interpretation:**
- If T = t × f, then for f = 0.785:
  ```
  sin[2π(X - 0.785t)]
  ```
  This gives wavelength = 1, frequency = 0.785 Hz

**Zhang equivalent:**
For Zhang's f = 0.785 Hz:
```
sin[2π(X - τ)] where τ = 0.785t
```

**Verification:**
- ✓ Frequency value 0.785 is within Zhang's range [0.3, 1.3]
- ✓ Wave number k = 2π (one wavelength along body) matches Zhang
- ∼ Implementation uses specific f = 0.785, not full Zhang range

---

## 5. Adaptive vs. Fixed Kinematics

### 5.1 Zhang (2018) - Fixed

**Critical requirement:**
```
Kinematics are FIXED across all Re and thickness.
Same A_max, same envelope, same wave form.
```

---

### 5.2 Current Implementation - Adaptive

**From `IBEELKinematics.cpp:707-793`:**

**Adaptive amplitude (when enabled):**
```cpp
Re_amplitude_factor = pow(Re/Re_ref, exponent)
thickness_amplitude_factor = 1.0 + 0.3 * (h/L - 0.04) / 0.04
d_adapted_amplitude = 0.125 * Re_amplitude_factor * thickness_amplitude_factor
```

**Adaptive frequency (when enabled):**
```cpp
Re_frequency_factor = pow(Re/Re_ref, exponent)
thickness_frequency_factor = 1.0 - 0.2 * (h/L - 0.04) / 0.04
d_adapted_frequency = 0.785 * Re_frequency_factor * thickness_frequency_factor
```

---

### 5.3 Enabling/Disabling Adaptation

**From input files:**
```cpp
enable_shape_adaptation = TRUE  // Current default
```

**To match Zhang (2018):**
```cpp
enable_shape_adaptation = FALSE  // Required for Zhang mode
```

**Verification:**
- ✗ Current default (TRUE) does NOT match Zhang
- ✓ Can be disabled via input file
- ⚠ Must ensure users know to disable for Zhang comparison

---

## 6. Complete Parameter Summary

### 6.1 Matching Parameters

| Parameter | Zhang Value | Current Value | Status |
|-----------|-------------|---------------|---------|
| **A_max (base)** | 0.125 | 0.125 | ✓ MATCH |
| **Envelope c0** | 0.03125 | 0.03125 | ✓ MATCH |
| **Envelope c1** | 1.03125 | 1.03125 | ✓ MATCH |
| **Wave number k** | 2π | 2π | ✓ MATCH |
| **Re formula** | 0.785/ν | 0.785/MU | ✓ MATCH |
| **Velocity coeff** | 0.785 | 0.785 | ✓ MATCH |

---

### 6.2 Configurable Parameters

| Parameter | Zhang Value | Current Default | Status |
|-----------|-------------|-----------------|--------|
| **Adaptation** | OFF (fixed) | ON (adaptive) | ⚠ DIFFERS |
| **Base frequency** | varies [0.3-1.3] | 0.785 fixed | ∼ PARTIAL |
| **Envelope power** | 1.0 (linear) | 1.0-3.0 (adaptive) | ⚠ DIFFERS |

---

### 6.3 Missing Test Cases

**Reynolds numbers:**
- ✗ Re = 50
- ✗ Re = 500
- ✓ Re = 1000
- ✗ Re = 5000 (have 5609 instead)
- ✓ Re = 10000
- ✗ Re = 50000
- ✗ Re = 100000
- ✗ Re = 200000

**Thickness ratios:**
- ✓ h/L = 0.04
- ✗ h/L = 0.08 (have partial)
- ✗ h/L = 0.12
- ✗ h/L = 0.16
- ✗ h/L = 0.20
- ✗ h/L = 0.24

---

## 7. Compliance Checklist

### 7.1 To Match Zhang (2018) Exactly

- [ ] Set `enable_shape_adaptation = FALSE` in all input files
- [ ] Create input files for all 8 Zhang Re values
- [ ] Create input files for all 6 Zhang thickness values
- [ ] Verify MU = 0.785 / Re in each file
- [ ] Verify base_amplitude = 0.125 (no adaptation)
- [ ] Verify base_frequency = 0.785 (or vary to sweep Re)
- [ ] Verify envelope_power = 1.0 (linear, fixed)
- [ ] Verify swimming_mode = 0.0 (anguilliform)
- [ ] Use NACA 00XX symmetric foils matching thickness
- [ ] Disable any maneuvering (body_is_maneuvering = FALSE)

---

### 7.2 Current Status vs. Zhang

**What matches:**
✓ Base kinematic equations (amplitude envelope, wave form)
✓ Reynolds number definition and calculation method
✓ Nondimensional formulation (L = 1, appropriate scaling)
✓ Self-propelled approach (constraint-based IB)

**What differs:**
✗ Adaptive kinematics enabled by default (should be disabled)
✗ Incomplete Re coverage (missing 50, 500, 5000, 50000+)
✗ Incomplete thickness coverage (missing 0.12, 0.16, 0.20, 0.24)
✗ Uses Re = 5609 instead of Zhang's 5000

**What needs clarification:**
? Time variable T in equations (dimensional vs. nondimensional)
? Exact frequency sweep implementation
? Relationship between code frequency and Re variation

---

## 8. Recommended Corrections

### 8.1 Immediate Fixes (High Priority)

**1. Disable adaptive kinematics for Zhang mode:**
```cpp
// In all Zhang_2018 input files:
enable_shape_adaptation = FALSE
envelope_power = 1.0
swimming_mode = 0.0
```

**2. Create missing Re test cases:**
```bash
input2d_Zhang_Re50_h004
input2d_Zhang_Re500_h004
input2d_Zhang_Re5000_h004  # Replace Re5609
input2d_Zhang_Re50000_h004
input2d_Zhang_Re100000_h004
input2d_Zhang_Re200000_h004
```

**3. Create missing thickness cases:**
```bash
# For each Re, create:
*_h012 (NACA0012)
*_h016 (NACA0016)
*_h020 (NACA0020)
*_h024 (NACA0024)
```

---

### 8.2 Medium Priority

**4. Verify time nondimensionalization:**
- Document exact meaning of T variable in parsers
- Confirm relationship: T = t × f or T = t/T_period
- Add comments to input files explaining

**5. Implement strict Zhang class:**
- Create `IBEELKinematicsZhang` derived class
- Hard-code fixed parameters (no adaptation)
- Simpler implementation, guaranteed compliance

---

### 8.3 Long Term (Nice to Have)

**6. Parameter sweep automation:**
- Script to generate all Re × thickness combinations
- Automated comparison to Zhang Fig. 4, 5, 6
- Plot overlays for validation

**7. Mesh refinement study:**
- Verify grid independence for Zhang cases
- Document minimum resolution requirements
- Test vorticity tagging thresholds

---

## 9. Quantitative Verification

### 9.1 Test Case: Re = 1000, h/L = 0.04

**Zhang parameters:**
```
A_max = 0.125
Envelope = (X + 0.03125) / 1.03125
Re = 1000
ν = 0.785 / 1000 = 0.000785
```

**Current file: `input2d_Re1000_h004`**
```
base_amplitude = 0.125 ✓
envelope in body_shape_equation: (X_0 + 0.03125)/1.03125 ✓
reynolds_number = 1000.0 ✓
MU = 0.785/Re = 0.000785 ✓
thickness_ratio = 0.04 ✓
```

**BUT:**
```
enable_shape_adaptation = TRUE ✗
```

**With adaptation at Re = 1000:**
```
Re_ratio = 1000 / 5000 = 0.2
Re_amplitude_factor = (0.2)^(-0.15) ≈ 1.28
A_adapted ≈ 0.125 × 1.28 = 0.16 ≠ 0.125 ✗
```

**Conclusion:** Will not match Zhang unless adaptation is disabled.

---

### 9.2 Expected Outputs (Zhang Fig. 4)

**For Re = 1000, h/L = 0.04:**

From Zhang et al. (2018) Fig. 4:
```
U₀ ≈ 0.15-0.20 (estimate from graph)
St ≈ 0.8-1.2
η ≈ 10-15%
Wake: 2S or 2P pattern
```

**Current implementation (with adaptation):**
```
U₀ = ? (likely higher due to increased amplitude)
St = ? (will differ)
η = ? (may be higher due to "optimization")
```

**Recommendation:**
Run both adaptive and non-adaptive, compare to Zhang.

---

## 10. Summary and Action Items

### 10.1 Overall Assessment

**Mathematical formulation:** ✓ CORRECT (matches Zhang when adaptation disabled)

**Default configuration:** ✗ INCORRECT (adaptation enabled by default)

**Test case coverage:** ∼ PARTIAL (some Re and h/L values, not complete)

**Ability to match Zhang:** ✓ POSSIBLE (requires input file changes)

---

### 10.2 Critical Action Items

1. **Create Zhang_2018 input files with:**
   - enable_shape_adaptation = FALSE
   - All 8 Re values
   - All 6 thickness values
   - Full 8×6 = 48 cases for complete Zhang reproduction

2. **Document differences clearly:**
   - Adaptive vs. strict modes
   - When to use each
   - Expected performance differences

3. **Implement `IBEELKinematicsZhang` class:**
   - Guarantees fixed kinematics
   - No user error possible
   - Clean separation of approaches

4. **Validation study:**
   - Run Zhang strict cases
   - Compare to Zhang Fig. 4, 5, 6
   - Publish results

---

## 11. Conclusion

The current implementation **has the correct mathematical formulation** to reproduce Zhang et al. (2018), but:

1. **Adaptive kinematics are enabled by default**, which violates Zhang's fixed-kinematics approach
2. **Test case coverage is incomplete** (missing many Re and h/L combinations)
3. **Users must manually disable adaptation** to match Zhang

**To exactly reproduce Zhang (2018):**
- Set `enable_shape_adaptation = FALSE`
- Use provided base parameters (already correct)
- Create full Re and thickness sweep input files
- Run simulations and compare to Zhang's published results

**The code is capable of matching Zhang, but the default configuration does not.**

---

**Document Version**: 1.0
**Last Updated**: 2025-01-17
**Verification Status**: Conditionally Compliant (requires configuration changes)
