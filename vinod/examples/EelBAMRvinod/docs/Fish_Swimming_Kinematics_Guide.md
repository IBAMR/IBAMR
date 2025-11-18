# Fish Swimming Kinematics: Anguilliform vs Carangiform

## Understanding the Swimming Wave Equation

### **Basic Equation: y(x,t) = A(x) · sin(kx - ωt)**

This equation describes the **lateral displacement** of the fish body during swimming.

#### **Components Explained:**

| Component | Symbol | Meaning | Units |
|-----------|--------|---------|-------|
| **Lateral displacement** | y(x,t) | How far the body moves sideways at position x and time t | m |
| **Amplitude envelope** | A(x) | Maximum lateral displacement at each position along body | m |
| **Wave number** | k = 2π/λ | Spatial frequency of the wave | rad/m |
| **Angular frequency** | ω = 2πf | Temporal frequency | rad/s |
| **Wavelength** | λ | Length of one complete wave | m |
| **Position along body** | x | Distance from head (x=0) to tail (x=L) | m |
| **Time** | t | Current time | s |

---

## The Key Difference: Amplitude Envelope A(x)

The **amplitude envelope A(x)** determines HOW MUCH each part of the body bends during swimming.

- **Head (x=0)**: Usually small amplitude (rigid head)
- **Mid-body**: Depends on swimming mode
- **Tail (x=L)**: Usually largest amplitude (propulsion)

---

## 1. ANGUILLIFORM Swimming (Eel-like)

### **Amplitude Envelope:**
```
A(X) = A₀ · exp[α(X - 1)]
```

Where:
- `X = x/L`: Normalized position (0 = head, 1 = tail)
- `A₀`: Base amplitude (typically 0.1L = 10% of body length)
- `α`: Growth rate parameter (typically 1.0 to 2.0)

### **Characteristics:**

| Property | Value | Description |
|----------|-------|-------------|
| **Amplitude growth** | Exponential | Smooth increase from head to tail |
| **Body participation** | Whole body | > 80% of body length undulates |
| **Wavelength** | λ < L | Multiple waves on body (0.6-0.8L) |
| **Wave number** | k > 2π/L | More than one complete wave |
| **Frequency** | 1-3 Hz | Moderate oscillation rate |
| **Swimming speed** | 0.5-1.5 BL/s | Body lengths per second |

### **Example Parameters:**
```
L = 1.0 m          # Body length
A₀ = 0.1 m         # Base amplitude (10% of length)
α = 1.0            # Exponential growth rate
λ = 0.8 m          # Wavelength (0.8L)
f = 2.0 Hz         # Frequency
```

### **Amplitude Distribution:**
```
Position (X)  |  0.0  |  0.2  |  0.4  |  0.6  |  0.8  |  1.0
Amplitude (m) | 0.037 | 0.045 | 0.055 | 0.067 | 0.082 | 0.100
```

The body gradually increases its lateral motion from head to tail.

### **Animals Using Anguilliform:**
- Eels (Anguilla)
- Lampreys
- Elongated fish
- Sea snakes

### **Advantages:**
✅ Excellent maneuverability
✅ Can swim in confined spaces
✅ Good for slow, precise movements
✅ Effective in vegetation

### **Disadvantages:**
❌ Lower efficiency (~60-65%)
❌ Slower maximum speed
❌ Higher energy cost per distance

---

## 2. CARANGIFORM Swimming (Tuna-like)

### **Amplitude Envelope:**
```
A(X) = a₀ + a₁X + a₂X²
```

Where:
- `X = x/L`: Normalized position (0 = head, 1 = tail)
- `a₀`: Constant term (typically ~0.01L)
- `a₁`: Linear term (typically ~-0.05L)
- `a₂`: Quadratic term (typically ~0.14L)

### **Characteristics:**

| Property | Value | Description |
|----------|-------|-------------|
| **Amplitude growth** | Quadratic | Rapid increase at posterior region |
| **Body participation** | Posterior body | < 50% of body length undulates significantly |
| **Wavelength** | λ > L | Less than one wave on body (1.2-2.0L) |
| **Wave number** | k < 2π/L | Fraction of one complete wave |
| **Frequency** | 2-5 Hz | Higher oscillation rate |
| **Swimming speed** | 1.0-4.0 BL/s | Body lengths per second |

### **Example Parameters:**
```
L = 1.0 m          # Body length
a₀ = 0.01 m        # Constant term
a₁ = -0.05 m       # Linear term (negative!)
a₂ = 0.14 m        # Quadratic term
λ = 1.5 m          # Wavelength (1.5L - longer than body!)
f = 3.0 Hz         # Frequency
```

### **Amplitude Distribution:**
```
Position (X)  |  0.0  |  0.2  |  0.4  |  0.6  |  0.8  |  1.0
Amplitude (m) | 0.010 | 0.012 | 0.022 | 0.040 | 0.066 | 0.100
```

The anterior body is nearly rigid, with most bending concentrated at the tail.

### **Animals Using Carangiform:**
- Tuna (Thunnus)
- Mackerel
- Jacks (Caranx)
- Many sharks
- Most fast-swimming fish

### **Advantages:**
✅ High efficiency (~75-85%)
✅ High maximum speed
✅ Lower energy cost per distance
✅ Effective for cruising

### **Disadvantages:**
❌ Reduced maneuverability
❌ Requires open water
❌ Less effective at low speeds
❌ Cannot navigate tight spaces

---

## Mathematical Comparison

### **Complete Swimming Equation:**

#### Anguilliform:
```
y(x,t) = [0.1 · exp(α(x/L - 1))] · sin(2π·x/0.8 - 2π·2·t)
```

#### Carangiform:
```
y(x,t) = [0.01 - 0.05(x/L) + 0.14(x/L)²] · sin(2π·x/1.5 - 2π·3·t)
```

---

## Velocity Calculation (for IBAMR)

The velocity of each point on the body is:

```
V(x,t) = ∂y/∂t = -A(x) · ω · cos(kx - ωt)
```

Where:
- `ω = 2πf`: Angular frequency
- The negative sign indicates the direction of motion

### **Implementation in IBAMR:**

```cpp
// Velocity calculation
double A = calculateAmplitude(X);  // A(X) based on swimming mode
double omega = 2.0 * M_PI * frequency;
double k = 2.0 * M_PI / wavelength;
double phase = k * x - omega * time;

// Lateral velocity (perpendicular to swimming direction)
double velocity_y = -A * omega * cos(phase);
```

---

## Performance Comparison

| Metric | Anguilliform | Carangiform |
|--------|--------------|-------------|
| **Efficiency** | 60-65% | 75-85% |
| **Max speed (BL/s)** | 1.5 | 4.0 |
| **Turn radius** | 0.2 BL | 1.5 BL |
| **Energy (J/m)** | Higher | Lower |
| **Acceleration** | Moderate | High |
| **Stability** | Moderate | High |

*BL = Body Length*

---

## Hydrodynamic Considerations

### **Anguilliform:**
- **Recoil motion**: Significant lateral head motion wastes energy
- **Added mass**: Entire body accelerates fluid
- **Wake structure**: Complex, unsteady vortices along entire body
- **Thrust production**: Distributed along body length

### **Carangiform:**
- **Minimal recoil**: Head stays nearly straight (efficient!)
- **Added mass**: Only posterior body accelerates fluid
- **Wake structure**: Concentrated reverse Kármán vortex street at tail
- **Thrust production**: Concentrated at tail (high efficiency)

---

## Choosing Swimming Mode for Simulation

### **Use Anguilliform when studying:**
- Low-speed maneuvering
- Confined spaces navigation
- Prey capture kinematics
- Elongated body mechanics
- Undulatory propulsion fundamentals

### **Use Carangiform when studying:**
- High-speed swimming performance
- Energy efficiency
- Open-water cruising
- Thrust optimization
- Tail propulsion mechanics

---

## Implementation in IBAMR ConstraintIB

### **Key Steps:**

1. **Define amplitude envelope function `A(X)`**
   - Anguilliform: `A(X) = A₀ · exp[α(X-1)]`
   - Carangiform: `A(X) = a₀ + a₁X + a₂X²`

2. **Set wave parameters**
   - Wavelength λ
   - Frequency f
   - Derived: k = 2π/λ, ω = 2πf

3. **Calculate position: `y(x,t) = A(x) · sin(kx - ωt)`**

4. **Calculate velocity: `V(x,t) = -A(x) · ω · cos(kx - ωt)`**

5. **Apply velocities to IB points in `setKinematicsVelocity()`**

---

## Parameter Tuning Guide

### **Frequency (f):**
- Higher f → faster swimming but more energy
- Typical range: 1-5 Hz
- Start with 2 Hz (anguilliform) or 3 Hz (carangiform)

### **Wavelength (λ):**
- Shorter λ → more waves on body
- Anguilliform: 0.6-0.8 L
- Carangiform: 1.2-2.0 L

### **Amplitude (A):**
- Higher A → more thrust but more drag
- Typical tail amplitude: 0.1-0.15 L
- Start with 0.1 L

### **Alpha (α) for Anguilliform:**
- Higher α → more concentration at tail
- Typical range: 0.5-2.0
- Start with 1.0

### **Quadratic coefficients for Carangiform:**
- Adjust a₀, a₁, a₂ to match experimental data
- Ensure A(0) is small (rigid head)
- Ensure A(1) = desired tail amplitude

---

## References

1. **Lighthill, M. J. (1971).** "Large-amplitude elongated-body theory of fish locomotion." *Proceedings of the Royal Society B*, 179(1055), 125-138.

2. **Videler, J. J., & Hess, F. (1984).** "Fast continuous swimming of two pelagic predators, saithe (Pollachius virens) and mackerel (Scomber scombrus): a kinematic analysis." *Journal of Experimental Biology*, 109(1), 209-228.

3. **Sfakiotakis, M., Lane, D. M., & Davies, J. B. C. (1999).** "Review of fish swimming modes for aquatic locomotion." *IEEE Journal of Oceanic Engineering*, 24(2), 237-252.

4. **Lauder, G. V., & Tytell, E. D. (2005).** "Hydrodynamics of undulatory propulsion." *Fish Physiology*, 23, 425-468.

5. **Borazjani, I., & Sotiropoulos, F. (2008).** "Numerical investigation of the hydrodynamics of carangiform swimming in the transitional and inertial flow regimes." *Journal of Experimental Biology*, 211(10), 1541-1558.

---

## Summary Table

| Aspect | Anguilliform | Carangiform |
|--------|--------------|-------------|
| **Equation** | A(X) = A₀·exp[α(X-1)] | A(X) = a₀+a₁X+a₂X² |
| **Body participation** | Whole body (>80%) | Posterior (< 50%) |
| **Wavelength** | Short (0.6-0.8L) | Long (1.2-2.0L) |
| **Frequency** | Moderate (1-3 Hz) | High (2-5 Hz) |
| **Speed** | Moderate (0.5-1.5 BL/s) | High (1.0-4.0 BL/s) |
| **Efficiency** | Lower (60-65%) | Higher (75-85%) |
| **Maneuverability** | Excellent | Moderate |
| **Best for** | Confined spaces | Open water cruising |
| **Examples** | Eels, lampreys | Tuna, mackerel, sharks |

---

## Next Steps

1. Review the **anguilliform example** in `examples/anguilliform/`
2. Review the **carangiform example** in `examples/carangiform/`
3. Choose the swimming mode appropriate for your research
4. Modify the amplitude envelope parameters
5. Generate geometry using MATLAB scripts
6. Compile and run IBAMR simulation
7. Visualize results in VisIt
8. Analyze thrust, efficiency, and flow patterns

---

**Happy swimming simulation! 🐟🐍**
