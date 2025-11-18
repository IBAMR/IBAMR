# Comparison: Adaptive Kinematics vs. Zhang et al. (2018) Strict Implementation

This document analyzes the key differences between the current **adaptive kinematics** implementation and the **strict Zhang et al. (2018)** specification.

---

## Executive Summary

| Aspect | Current Implementation (Adaptive) | Zhang (2018) Strict |
|--------|-----------------------------------|---------------------|
| **Amplitude** | Re-dependent, thickness-dependent | **Fixed**: A_max = 0.125 |
| **Frequency** | Re-dependent, thickness-dependent | **Fixed in kinematics** (varied to change Re) |
| **Envelope** | Swimming mode-dependent power law | **Fixed**: linear (X + 0.03125)/1.03125 |
| **Re variation** | Modify kinematics + vary ν or f | **Only vary ν or f** |
| **Strouhal** | Can be prescribed or emergent | **Emergent output** |
| **Philosophy** | Bio-inspired adaptation | Controlled parameter study |

---

## 1. Kinematic Law Comparison

### 1.1 Zhang (2018) - Fixed Kinematics

**Midline equation:**
```
Y(X, τ) = 0.125 × [(X + 0.03125) / 1.03125] × sin[2π(X - τ)]
```

**Properties:**
- A_max = 0.125 (constant across all Re and thickness)
- Linear envelope from head to tail
- One wavelength along body (k = 2π/L = 2π)
- **No adaptation**: Same kinematics for Re = 50 and Re = 2×10⁵

---

### 1.2 Current Implementation - Adaptive Kinematics

**Location in code:** `IBEELKinematics.cpp:707-793`

**Adaptive amplitude:**
```cpp
// Lines 720-737
Re_amplitude_factor = pow(Re/Re_ref, exponent)
thickness_amplitude_factor = 1.0 + 0.3 * (h/L - 0.04) / 0.04
d_adapted_amplitude = 0.125 * Re_amplitude_factor * thickness_amplitude_factor
```

**Scaling behavior:**
- **Low Re (< 5000)**: Amplitude ∝ Re^(-0.15) → amplitude **increases** as Re decreases
- **High Re (≥ 5000)**: Amplitude ∝ Re^(-0.08) → amplitude **decreases** as Re increases
- **Thickness effect**: +30% amplitude for h/L = 0.08 vs. 0.04

**Adaptive frequency:**
```cpp
// Lines 741-756
Re_frequency_factor = pow(Re/Re_ref, exponent)
thickness_frequency_factor = 1.0 - 0.2 * (h/L - 0.04) / 0.04
d_adapted_frequency = 0.785 * Re_frequency_factor * thickness_frequency_factor
```

**Scaling behavior:**
- **Low Re (< 5000)**: Frequency ∝ Re^(+0.12) → frequency **decreases** at low Re
- **High Re (≥ 5000)**: Frequency ∝ Re^(+0.05) → frequency **increases** at high Re
- **Thickness effect**: -20% frequency for h/L = 0.08 vs. 0.04

**Adaptive envelope:**
```cpp
// Lines 758-764
d_envelope_power = 1.0 + 2.0 * swimming_mode
// Anguilliform (mode=0): power = 1.0 (linear)
// Carangiform (mode=1): power = 3.0 (cubic)
```

---

## 2. Physical Interpretation Differences

### 2.1 Zhang Philosophy: Controlled Observation

**Approach:**
1. Fix the kinematic pattern (shape of undulation)
2. Systematically vary flow regime (Re) and geometry (thickness)
3. Observe how the *same motion* performs differently across regimes
4. **Question**: "How does Re affect performance of a given swimming pattern?"

**Analogy:**
Testing a fixed wing design (NACA 0012) in wind tunnels at different speeds.

**Scientific value:**
- Isolates effect of flow regime from behavioral adaptation
- Clean parameter study
- Direct comparison across Re (same kinematics)
- Matches classical fluid mechanics approach

---

### 2.2 Adaptive Philosophy: Bio-Inspired Optimization

**Approach:**
1. Recognize that real fish change their swimming pattern with speed/size
2. Model kinematic adaptation as function of Re and thickness
3. Each Re uses "optimized" kinematics for that regime
4. **Question**: "What kinematics are best for each Re and thickness?"

**Analogy:**
Birds changing wing beat frequency and amplitude based on flight speed.

**Scientific value:**
- Models realistic biological behavior
- May achieve better performance at each Re
- Explores adaptive strategies
- Closer to natural swimming

---

## 3. Quantitative Comparison

### 3.1 Example: Re = 1000, h/L = 0.04

**Zhang (2018) strict:**
```
A_max = 0.125
f_nondim = 1.0 (in wave equation)
Envelope: (X + 0.03125) / 1.03125
```

**Adaptive implementation:**
```
Re_ratio = 1000 / 5000 = 0.2
Re_amplitude_factor = (0.2)^(-0.15) = 1.28
thickness_amplitude_factor = 1.0 (at reference)
A_adapted = 0.125 × 1.28 × 1.0 = 0.160  (+28%)

Re_frequency_factor = (0.2)^(0.12) = 0.78
thickness_frequency_factor = 1.0 (at reference)
f_adapted = 0.785 × 0.78 × 1.0 = 0.612  (-22%)

Envelope power = 1.0 + 2.0 × 0.0 = 1.0 (same)
```

**Result:** Adaptive uses **larger amplitude**, **lower frequency** at low Re.

---

### 3.2 Example: Re = 100,000, h/L = 0.08

**Zhang (2018) strict:**
```
A_max = 0.125
f_nondim = 1.0
Envelope: (X + 0.03125) / 1.03125
```

**Adaptive implementation:**
```
Re_ratio = 100000 / 5000 = 20.0
Re_amplitude_factor = (20.0)^(-0.08) = 0.79
thickness_amplitude_factor = 1.0 + 0.3 × (0.08 - 0.04) / 0.04 = 1.30
A_adapted = 0.125 × 0.79 × 1.30 = 0.129  (+3% net)

Re_frequency_factor = (20.0)^(0.05) = 1.17
thickness_frequency_factor = 1.0 - 0.2 × 0.04 / 0.04 = 0.80
f_adapted = 0.785 × 1.17 × 0.80 = 0.735  (-6%)

Envelope power = 1.0 + 2.0 × swimming_mode  (depends on mode)
```

**Result:** At high Re, thickness effect dominates; still differs from Zhang.

---

### 3.3 Summary Table

| Re | h/L | Zhang A_max | Adaptive A | Δ% | Zhang f | Adaptive f | Δ% |
|----|-----|-------------|------------|-----|---------|------------|-----|
| 500 | 0.04 | 0.125 | ~0.175 | +40% | 1.0 | ~0.56 | -44% |
| 1000 | 0.04 | 0.125 | ~0.160 | +28% | 1.0 | ~0.61 | -39% |
| 5000 | 0.04 | 0.125 | 0.125 | 0% | 1.0 | 0.785 | -22% |
| 10000 | 0.08 | 0.125 | ~0.136 | +9% | 1.0 | ~0.74 | -26% |
| 100000 | 0.08 | 0.125 | ~0.129 | +3% | 1.0 | ~0.74 | -26% |

**Key observations:**
- Largest differences at **low Re** (viscous regime)
- Adaptive always uses **lower frequency** (nondimensionally)
- Adaptive compensates low Re with **higher amplitude**

---

## 4. Reynolds Number Variation

### 4.1 Zhang (2018) Approach

**Fixed kinematics, vary Re by:**

**Method A: Vary dimensional frequency**
```
Re = V_max L / ν = 2π A_max f L / ν
For fixed A_max, L, ν:
  Re ∝ f
  f = Re × ν / (2π × 0.125 × L) = Re × ν / 0.785
```

**Method B: Vary kinematic viscosity**
```
For fixed A_max, f, L:
  ν = 2π A_max f L / Re = 0.785 f / Re
```

**Zhang likely uses Method B in simulations** (easier to keep f constant).

---

### 4.2 Adaptive Approach

**Modified kinematics + vary Re:**

The adaptive code changes BOTH kinematics AND Re:
1. Set desired Re via input file
2. Calculate adapted A and f based on Re
3. Use different kinematic pattern at each Re
4. Vary ν or dimensional f to achieve target Re

**Problem:** This conflates two effects:
- Effect of Re on flow physics (what Zhang studies)
- Effect of different swimming patterns (behavioral adaptation)

---

## 5. Expected Performance Differences

### 5.1 Thrust and Propulsion

**Zhang (fixed):**
- At low Re: Weaker thrust (small amplitude, viscous losses)
- At high Re: Better thrust (inertial regime)
- Performance varies **only due to flow regime change**

**Adaptive:**
- At low Re: Compensates with larger amplitude → may maintain thrust
- At high Re: Smaller amplitude → potentially better efficiency
- Performance varies due to **flow regime + kinematic adaptation**

**Expected result:**
- Adaptive likely shows **flatter performance curves** across Re
- Zhang shows **stronger Re dependence** (by design)

---

### 5.2 Strouhal Number

**Zhang (fixed):**
```
St = 2 A_max f / U_0
Since A_max and f are fixed, St varies only if U_0 changes.
St ∈ [0.3, 2.1] observed range
```

**Adaptive:**
```
St = 2 A_adapted(Re) f_adapted(Re) / U_0
Both numerator and denominator vary with Re.
Could target specific St range (e.g., optimal 0.2-0.4).
```

**Key difference:**
- Zhang: St is pure output (emergent)
- Adaptive: Can implicitly target St via adaptation rules

---

### 5.3 Efficiency

**Zhang (fixed):**
- Efficiency η(Re) reflects how well **this fixed pattern** works at each Re
- Peak efficiency at "natural" Re for this pattern
- Mismatch at low/high Re

**Adaptive:**
- Efficiency η(Re) reflects **optimized pattern** at each Re
- Potentially high efficiency across broader Re range
- But: Is this what we want to study?

---

## 6. Scientific Questions Answered

### 6.1 Zhang (2018) Answers

✓ How does **flow regime** (Re) affect a given swimming pattern?
✓ What is the **optimal thickness** for a fixed kinematic strategy?
✓ How do **wake structures** change with Re for BCF swimmers?
✓ What Re range transitions from viscous to inertial propulsion?
✓ Fundamental fluid-structure interaction at different Re

---

### 6.2 Adaptive Implementation Answers

✓ How should **kinematics adapt** to different Re and thickness?
✓ What are **bio-inspired scaling laws** for swimming?
✓ Can we achieve **robust performance** across Re ranges?
✓ How do real fish **optimize their swimming** in different regimes?
✓ Adaptive control strategies for bio-inspired robots

---

## 7. When to Use Each Approach

### 7.1 Use Zhang (2018) Strict Mode When:

✓ Studying **fundamental Re effects** on fluid dynamics
✓ Comparing to **Zhang's published results**
✓ Need **controlled parameter study** (isolate Re effect)
✓ Benchmarking and validation against literature
✓ Understanding transition regimes (viscous → inertial)
✓ Clean separation of geometry vs. flow regime effects

**Example applications:**
- Academic validation studies
- Fundamental research on Re scaling
- Comparing different numerical methods
- Reproducing Zhang's Fig. 4, Fig. 6, etc.

---

### 7.2 Use Adaptive Kinematics Mode When:

✓ Designing **bio-inspired robots** that adapt to conditions
✓ Modeling **realistic fish behavior** across speeds
✓ Optimizing swimming for **specific performance metrics**
✓ Exploring **control strategies** for AUVs
✓ Engineering applications (maximize efficiency)
✓ Studying behavioral adaptation in swimming animals

**Example applications:**
- Autonomous underwater vehicle design
- Robotic fish control algorithms
- Biomimetic propulsion systems
- Energy-efficient swimming strategies

---

## 8. Code Implementation Differences

### 8.1 Zhang (2018) Strict Implementation

**Required changes to current code:**

```cpp
// In IBEELKinematics.cpp

void IBEELKinematics::calculateAdaptiveKinematics(const double time)
{
    // ZHANG MODE: No adaptation, use base values
    d_adapted_amplitude = d_base_amplitude;   // Always 0.125
    d_adapted_frequency = d_base_frequency;    // Always 0.785 (or 1.0 nondim)
    d_adapted_wavelength = 1.0;                // Always 1.0
    d_envelope_power = 1.0;                    // Always linear

    // No Re-dependent or thickness-dependent modification
}
```

**Input file changes:**
```
# Zhang (2018) strict mode
enable_shape_adaptation = FALSE
base_amplitude = 0.125
base_frequency = 0.785  # or 1.0 in nondimensional time
envelope_power = 1.0
swimming_mode = 0.0  # Anguilliform (linear envelope)

# Vary Re by changing viscosity
reynolds_number = 1000  # Target Re
# Solver calculates: viscosity = 0.785 / Re
```

---

### 8.2 Current Adaptive Implementation

**Key functions:**

**`calculateAdaptiveKinematics()` (lines 707-793):**
- Calculates Re_amplitude_factor
- Calculates thickness_amplitude_factor
- Computes adapted amplitude
- Calculates Re_frequency_factor
- Calculates thickness_frequency_factor
- Computes adapted frequency
- Adjusts envelope power based on swimming_mode

**Input parameters:**
```cpp
d_reynolds_number         // Target Re
d_thickness_ratio         // h/L
d_base_amplitude          // Reference amplitude (0.125)
d_base_frequency          // Reference frequency (0.785)
d_swimming_mode           // 0=anguilliform, 1=carangiform
d_enable_shape_adaptation // TRUE = adaptive, FALSE = fixed
```

---

## 9. Validation and Benchmarking

### 9.1 Zhang (2018) Strict Mode Validation

**Target results from Zhang et al. (2018):**

**Figure 4 - Swimming speed vs. Re:**
- U₀ increases with Re
- Different curves for each thickness
- Transition around Re ~ 10³-10⁴

**Figure 5 - Strouhal number vs. Re:**
- St peaks around Re ~ 10³
- St ∈ [0.3, 2.1]
- Thickness affects St

**Figure 6 - Efficiency vs. Re:**
- η increases with Re (to a point)
- Optimal thickness depends on Re
- Transition from viscous to inertial regime

**Wake structures (Figure 7):**
- Re = 50: Weak, diffuse
- Re = 500-1000: 2S pattern
- Re ≥ 10⁴: 2P pattern
- Re ≥ 5×10⁴: Multiple vortices, KH instability

**Validation criteria:**
- [ ] U₀(Re) matches Fig. 4 trends
- [ ] St(Re) matches Fig. 5 peaks
- [ ] η(Re) matches Fig. 6 efficiency curves
- [ ] Wake patterns match Re regimes
- [ ] Thickness effects match experimental observations

---

### 9.2 Adaptive Mode Validation

**No direct Zhang comparison possible** (different approach).

**Validation against:**
- General fish swimming data (Videler, Lauder, etc.)
- Strouhal number optimization (St ~ 0.2-0.4 for efficiency)
- Biological scaling laws (size-speed relationships)
- Engineering performance metrics (COT, efficiency)

---

## 10. Recommended Workflow

### Phase 1: Validate with Zhang Strict Mode
1. Implement strict Zhang mode (disable adaptation)
2. Run Re sweep: {50, 500, 1000, 5000, 10⁴, 10⁵, 2×10⁵}
3. Run thickness sweep: {0.04, 0.08, 0.12, 0.16, 0.20, 0.24}
4. Compare U₀, St, η to Zhang Fig. 4-6
5. Verify wake patterns
6. **Purpose:** Establish baseline, validate solver

### Phase 2: Explore Adaptive Kinematics
1. Enable adaptive mode
2. Same Re and thickness sweeps
3. Compare performance to strict mode
4. Quantify improvement from adaptation
5. Identify optimal adaptation strategies
6. **Purpose:** Investigate bio-inspired adaptation

### Phase 3: Combined Analysis
1. Analyze when adaptation helps most (low Re? high thickness?)
2. Determine regime-dependent strategies
3. Develop control algorithms for robots
4. Publish comparative study
5. **Purpose:** Advance knowledge of both approaches

---

## 11. Recommendations

### For Fundamental Research:
**Use Zhang (2018) strict mode**
- Isolates Re effects
- Clean comparison to literature
- Better for publishing fundamental studies

### For Engineering Applications:
**Use adaptive kinematics**
- Better performance across operating range
- Realistic for autonomous vehicles
- Matches biological systems

### For Comprehensive Study:
**Implement both, compare**
- Create `IBEELKinematicsZhang` class (strict)
- Keep `IBEELKinematics` class (adaptive)
- Run both on same cases
- Quantify benefit of adaptation
- Publish comparative paper

---

## 12. Summary Table

| Criterion | Zhang Strict | Adaptive | Winner |
|-----------|--------------|----------|--------|
| **Literature validation** | ✓✓✓ | ✗ | Zhang |
| **Bio-realism** | ✗ | ✓✓✓ | Adaptive |
| **Parameter isolation** | ✓✓✓ | ✗ | Zhang |
| **Engineering performance** | ✓ | ✓✓✓ | Adaptive |
| **Fundamental science** | ✓✓✓ | ✓ | Zhang |
| **Robotic control** | ✗ | ✓✓✓ | Adaptive |
| **Ease of comparison** | ✓✓✓ | ✗ | Zhang |
| **Broader Re range** | ✓ | ✓✓✓ | Adaptive |

---

## References

1. Zhang, C., Huang, H., & Lu, X.-Y. (2018). Effects of Reynolds number and thickness on an undulatory self-propelled foil. *Physics of Fluids*, **30**, 071902.

2. Videler, J. J. (1993). *Fish Swimming*. Chapman & Hall.

3. Lauder, G. V., & Tytell, E. D. (2006). Hydrodynamics of undulatory propulsion. *Fish Biomechanics*, 425-468.

4. Taylor, G. K., Nudds, R. L., & Thomas, A. L. R. (2003). Flying and swimming animals cruise at a Strouhal number tuned for high power efficiency. *Nature*, **425**, 707-711.

---

**Document Version**: 1.0
**Last Updated**: 2025-01-17
**Status**: Technical Comparison
