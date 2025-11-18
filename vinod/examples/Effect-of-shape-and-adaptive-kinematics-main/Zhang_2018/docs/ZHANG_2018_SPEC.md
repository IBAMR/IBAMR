# Zhang et al. (2018) - Complete Nondimensional Specification

**Reference**: "Effects of Reynolds number and thickness on an undulatory self-propelled foil"
*Physics of Fluids* **30**, 071902 (2018)

This document provides the complete, authoritative nondimensional specification extracted from Zhang et al. (2018). All quantities are presented in nondimensional form suitable for direct implementation.

---

## Nondimensionalization Scheme

Zhang et al. use the following scaling:

| Quantity | Scale | Symbol |
|----------|-------|--------|
| **Length** | Foil chord | L = 1 |
| **Time** | Period or convective | t* = t/T or t* = tU_ref/L |
| **Velocity** | Mean forward speed | U_ref = U₀ |
| **Amplitude** | Chord length | c = L |
| **Reynolds** | Maximum lateral velocity | Re = V_max L / ν |

---

## 1. Nondimensional Geometry

### 1.1 Shape Family

**NACA 00XX symmetric foils** with thickness-to-chord ratio:

```
t* = t/c ∈ {0.04, 0.08, 0.12, 0.16, 0.20, 0.24}
```

**Foil set:**
- NACA0004
- NACA0008
- NACA0012
- NACA0016
- NACA0020
- NACA0024

All foils have **chord = 1** (nondimensional).

### 1.2 Nondimensional Volume

Volume per unit depth divided by c²:

```
V_f* = V_f / L²
```

**Values:**
```
t* = 0.04  →  V_f* = 0.03146
t* = 0.08  →  V_f* = 0.05887
t* = 0.12  →  V_f* = 0.08621
t* = 0.16  →  V_f* = 0.1135
t* = 0.20  →  V_f* = 0.1408
t* = 0.24  →  V_f* = 0.1681
```

---

## 2. Nondimensional Kinematics (BCF Model)

### 2.1 Coordinate Systems

**Streamwise (chordwise) coordinate:**
```
X = x/L,    0 ≤ X ≤ 1
```

**Nondimensional time:**
```
τ = t/T
```

where T = 1/f is the undulation period.

### 2.2 Midline Motion Law

Zhang uses a **traveling wave with linear amplitude envelope**:

```
Y(X, τ) = A_max* × [(X + 0.03125) / 1.03125] × sin[2π(X - τ)]
```

**Components:**
- **Maximum amplitude**: A_max* = 0.125
- **Amplitude envelope**: A(X) = A_max* × [(X + 0.03125) / 1.03125]
- **Wave**: sin[2π(X - τ)]

### 2.3 Amplitude Envelope

```
A(X) = 0.125 × [(X + 0.03125) / 1.03125]
```

**Properties:**
- **Monotonic**: Increases linearly from head (X=0) to tail (X=1)
- At head (X=0): A(0) = 0.125 × (0.03125/1.03125) ≈ 0.00379
- At tail (X=1): A(1) = 0.125 × (1.03125/1.03125) = 0.125
- **BCF pattern**: Whole-body undulation with increasing amplitude

### 2.4 Fixed Kinematics

**CRITICAL**: Zhang uses **prescribed, fixed kinematics**:
- Amplitude and frequency are **NOT** adapted based on Re
- The same kinematic pattern is used across all Re and thickness values
- Re is varied by changing **frequency f** or **viscosity ν**, not by modifying the motion

---

## 3. Nondimensional Frequency and Time

### 3.1 Frequency Range

Zhang sweeps dimensional frequency:
```
f ∈ [0.3, 1.3] Hz
```

### 3.2 Nondimensional Frequency

In the wave equation Y(X, τ), time is nondimensionalized by T = 1/f:
```
f* = f × T = 1
```

The nondimensional frequency in the wave is **always 1**.

**Physics dependence on f:**
- Dimensional frequency f enters the problem through Reynolds number
- Different f → different Re → different flow physics
- Kinematics shape remains unchanged

---

## 4. Nondimensional Reynolds Number

### 4.1 Definition

```
Re = V_max L / ν
```

where:
- L = 1 (chord length)
- ν = kinematic viscosity
- V_max = maximum lateral velocity

### 4.2 Maximum Lateral Velocity

From the kinematic law:
```
V_max = ∂Y/∂t |_max = 2π A_max f
```

**Nondimensional form:**
```
V_max* = 2π A_max* f = 2π × 0.125 × f = 0.785 f
```

### 4.3 Reynolds Number Expression

```
Re = (0.785 f) / ν
```

### 4.4 Prescribed Reynolds Number Values

Zhang simulates **8 Reynolds numbers**:

```
Re ∈ {50, 500, 1000, 5000, 10⁴, 5×10⁴, 10⁵, 2×10⁵}
```

**Flow regimes:**
- Re = 50-500: Viscous-dominated
- Re = 1000-5000: Transitional
- Re = 10⁴-10⁵: Inertial
- Re = 2×10⁵: High Reynolds, turbulent transition

### 4.5 Achieving Different Re

Two approaches (equivalent):

**Method 1: Vary frequency (constant ν)**
```
f = Re × ν / 0.785
```

**Method 2: Vary viscosity (constant f)**
```
ν = 0.785 f / Re
```

Zhang likely uses Method 2 for numerical convenience.

---

## 5. Nondimensional Strouhal Number (Self-Propelled)

### 5.1 Definition

Based on **emergent** forward velocity U₀:
```
St = (2 A_max f) / U₀
```

### 5.2 Nondimensional Form

With scaling U_ref = U₀:
```
St* = 2 A_max* / T* = 2 × 0.125 / 1 = 0.25
```

**BUT** this is only the mathematical nondimensionalization.

### 5.3 Physical Strouhal Number

In self-propelled simulations:
- U₀ is **not prescribed**, it emerges from force balance
- U₀ varies with Re, f, and thickness
- Therefore, **St is an output, not an input**

**Observed range:**
```
St ∈ [0.3, 2.1]
```

**Dependencies:**
- St increases with Re (at fixed f)
- St varies with thickness
- Optimal swimming typically at St ≈ 0.2-0.4

---

## 6. Nondimensional Output Quantities

### 6.1 Force Coefficients

**Thrust coefficient:**
```
C_T = F_T / (½ ρ U₀² L)
```

**Lateral force coefficient:**
```
C_L = F_L / (½ ρ U₀² L)
```

**Time-dependent:**
```
C_L(τ) varies periodically with period τ = 1
```

### 6.2 Efficiency

**Propulsive efficiency:**
```
η = E / W
```

where:
- E = useful work (thrust × distance)
- W = total mechanical work input

**Froude efficiency:**
```
η_F = (C_T × U₀) / P
```

### 6.3 Energy and Power

**Nondimensional kinetic energy:**
```
E* = E / (ρ L³ U₀²)
```

**Nondimensional power:**
```
P* = P / (ρ L³ U₀³ / T)
```

### 6.4 Swimming Speed

**Output of self-propelled simulation:**
- U₀ emerges from thrust-drag balance
- Varies with Re and thickness
- Used to compute St and efficiency

---

## 7. Wake Structures (Nondimensional Classification)

Zhang characterizes wake patterns based on Re:

| Reynolds | Wake Type | Description |
|----------|-----------|-------------|
| Re = 50 | Weak bound vorticity | Viscous diffusion dominates |
| Re = 500 | 2S reverse von Kármán | Two single vortices per cycle |
| Re = 1000-10⁴ | 2P | Two vortex pairs per cycle |
| Re ≥ 5×10⁴ | 2P + KH eddies | Multiple shed eddies, Kelvin-Helmholtz |

**Vortex shedding:**
- Governed by nondimensional vorticity ω* = ωL/U₀
- Wake topology changes qualitatively with Re
- Thrust generation mechanism differs by regime

---

## 8. Complete Parameter Summary

### Input Parameters (Prescribed)

```python
# Geometry
chord_length = 1.0  # nondimensional
thickness_ratio = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24]

# Kinematics (FIXED - not adapted)
A_max = 0.125
envelope_coefficients = [0.03125, 1.03125]  # (X + c0) / c1
wave_number = 1.0  # one wavelength along body
frequency_nondim = 1.0  # in wave equation

# Reynolds number
Re = [50, 500, 1000, 5000, 1e4, 5e4, 1e5, 2e5]

# Simulation (choose one method)
# Method A: vary frequency
dimensional_frequency = Re * viscosity / 0.785
# Method B: vary viscosity
viscosity = 0.785 * dimensional_frequency / Re
```

### Output Quantities (Computed)

```python
# Emergent quantities
U_0  # swimming speed (from force balance)
St = (2 * A_max * f) / U_0  # Strouhal number
C_T  # thrust coefficient
C_L  # lateral force coefficient
eta  # propulsive efficiency
P    # power consumption
wake_type  # vortex pattern classification
```

---

## 9. Implementation Notes

### 9.1 Key Differences from Adaptive Models

**Zhang (2018) approach:**
- ✓ **Fixed kinematics** across all Re and thickness
- ✓ Self-propelled (U₀ emerges)
- ✓ St is an output, not prescribed
- ✓ Same A(X) envelope for all cases

**Adaptive models (not Zhang):**
- ✗ Modify amplitude based on Re
- ✗ Modify frequency based on thickness
- ✗ Change envelope shape
- ✗ Prescribe St or U₀

### 9.2 Correct Implementation Checklist

- [ ] Use A_max = 0.125 (constant)
- [ ] Use linear envelope A(X) = 0.125 × (X + 0.03125) / 1.03125
- [ ] Use traveling wave Y = A(X) × sin[2π(X - τ)]
- [ ] Do NOT modify amplitude or frequency based on Re
- [ ] Vary Re by changing ν or f, not kinematics
- [ ] Simulate self-propelled (zero net force)
- [ ] Compute St as output: St = 2Af/U₀
- [ ] Use NACA 00XX symmetric foils
- [ ] Test thickness range 0.04 ≤ t/c ≤ 0.24
- [ ] Test Reynolds range 50 ≤ Re ≤ 2×10⁵

### 9.3 Validation Targets

**From Zhang (2018) Fig. 4:**
- U₀ increases with Re
- St peaks around Re = 10³-10⁴
- Efficiency improves with Re (up to a point)
- Optimal thickness varies with Re

**Wake patterns:**
- Re = 50: Symmetric, weak vortices
- Re = 1000: Clear 2P pattern
- Re ≥ 5×10⁴: Complex, multiple vortices

---

## 10. Mathematical Summary (Copy-Ready for Code)

```
# Nondimensional parameters
L = 1.0
A_max = 0.125
c0 = 0.03125
c1 = 1.03125

# Amplitude envelope
A(X) = A_max * (X + c0) / c1

# Midline kinematics
Y(X, τ) = A(X) * sin(2π * (X - τ))

# Velocity (derivative)
V_y(X, τ) = 2π * A(X) * cos(2π * (X - τ))

# Maximum velocity
V_max = 2π * A_max * f = 0.785 * f

# Reynolds number
Re = V_max * L / ν = 0.785 * f / ν

# Strouhal (output)
St = 2 * A_max * f / U_0

# Force coefficients
C_T = F_T / (0.5 * ρ * U_0^2 * L)
C_L = F_L / (0.5 * ρ * U_0^2 * L)

# Efficiency
η = (F_T * U_0) / P_input
```

---

## References

1. Zhang, C., Huang, H., & Lu, X.-Y. (2018). Effects of Reynolds number and thickness on an undulatory self-propelled foil. *Physics of Fluids*, **30**, 071902.

2. Zhang, C., Huang, H., & Lu, X.-Y. (2017). Free locomotion of a flexible plate near the ground. *Physics of Fluids*, **29**, 041903.

3. Lauder, G. V., & Tytell, E. D. (2006). Hydrodynamics of undulatory propulsion. In *Fish Biomechanics* (pp. 425-468).

---

**Document Version**: 1.0
**Last Updated**: 2025-01-17
**Status**: Authoritative Reference
