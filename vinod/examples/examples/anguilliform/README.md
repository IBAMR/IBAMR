# Anguilliform (Eel-like) Swimming Simulation

## Overview

This example demonstrates **anguilliform swimming** using IBAMR's ConstraintIB method. Anguilliform swimming is characterized by:

- **Whole body undulation**: The entire body participates in the swimming motion
- **Exponential amplitude envelope**: Amplitude increases exponentially from head to tail
- **Multiple waves**: More than one wavelength fits on the body (λ < L)
- **Examples in nature**: Eels, lampreys, sea snakes

## Kinematics

### Amplitude Envelope
```
A(X) = 0.1 × exp[α(X - 1)]
```

Where:
- `X = x/L`: Normalized position along body (0 = head, 1 = tail)
- `0.1`: Base amplitude (10% of body length)
- `α = 1.0`: Exponential growth rate

### Complete Swimming Equation
```
y(x,t) = A(x) · sin(kx - ωt)
      = [0.1 · exp(α(x/L - 1))] · sin(2π·x/λ - 2π·f·t)
```

### Default Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Frequency | f | 2.0 Hz | Tail beat frequency |
| Wavelength | λ | 0.8 L | 80% of body length |
| Base amplitude | A₀ | 0.1 L | 10% of body length |
| Growth rate | α | 1.0 | Exponential parameter |
| Body length | L | 1.0 m | Fish length |
| Body width | w_max | 0.04 m | 4% of length |

### Amplitude Distribution

```
Position (X):    0.0   0.2   0.4   0.6   0.8   1.0
Amplitude (m): 0.037 0.045 0.055 0.067 0.082 0.100
```

The body gradually increases its lateral motion from head to tail.

## File Structure

```
anguilliform/
├── IBEELKinematics_Anguilliform.cpp   # C++ kinematics implementation
├── IBEELKinematics.h                   # Header file
├── generate_eel_geometry.m             # MATLAB geometry generator
├── input2d_anguilliform                # IBAMR input configuration
└── README.md                           # This file
```

## Usage Instructions

### Step 1: Generate Geometry

Run MATLAB to create the eel body geometry:

```bash
matlab -batch "generate_eel_geometry"
```

This creates:
- `eel2d.vertex`: Lagrangian mesh file (256 points)
- `eel_geometry.png`: Visualization of body shape

### Step 2: Set Up IBAMR Project

Copy files to your IBAMR project:

```bash
# Create project directory
mkdir -p ~/ibamr-projects/anguilliform
cd ~/ibamr-projects/anguilliform

# Copy source files
cp IBEELKinematics_Anguilliform.cpp IBEELKinematics.cpp
cp IBEELKinematics.h .
cp input2d_anguilliform input2d
cp eel2d.vertex .

# Copy example main program (you need example.cpp from IBAMR examples)
cp $IBAMR_ROOT/examples/ConstraintIB/eel2d/example.cpp .
```

### Step 3: Create CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.10)
project(anguilliform)

find_package(IBAMR REQUIRED)

add_executable(main2d
    example.cpp
    IBEELKinematics.cpp
)

target_link_libraries(main2d IBAMR::IBAMR2d)
```

### Step 4: Build

```bash
mkdir build
cd build
cmake ..
make -j4
```

### Step 5: Run Simulation

```bash
cd ..  # Back to project root
mpirun -np 4 ./build/main2d input2d
```

**Simulation details:**
- Duration: 10 seconds (20 swimming cycles at 2 Hz)
- Time step: 1 ms
- Output interval: Every 100 steps (0.1 seconds)
- Computational time: ~2-4 hours on 4 cores

### Step 6: Visualize Results

Use VisIt to view results:

```bash
visit -o viz_anguilliform/dumps.visit
```

**What to visualize:**
- **Velocity magnitude**: Shows flow field around swimming eel
- **Vorticity**: Reveals vortex shedding patterns
- **Pressure**: Shows pressure distribution
- **Mesh**: Shows deforming Lagrangian structure

## Parameter Tuning

### Frequency (f)

Controls how fast the eel beats its tail:

```
f = 1.5 Hz  →  Slow swimming  →  More efficient
f = 2.0 Hz  →  Moderate (default)
f = 3.0 Hz  →  Fast swimming  →  Higher thrust
```

**In input2d:** Modify `kinematics_velocity_function { frequency = 2.0 }`

### Wavelength (λ)

Controls number of waves on body:

```
λ = 0.6 L  →  1.67 waves on body  →  More undulation
λ = 0.8 L  →  1.25 waves (default)
λ = 1.0 L  →  1.0 wave on body
```

**In input2d:** Modify `kinematics_velocity_function { wavelength = 0.8 }`

### Amplitude (A₀)

Controls maximum tail deflection:

```
A₀ = 0.08 L  →  Smaller amplitude  →  Less thrust
A₀ = 0.10 L  →  Default (10%)
A₀ = 0.12 L  →  Larger amplitude  →  More thrust
```

**In input2d:** Modify `kinematics_velocity_function { amplitude = 0.1 }`

### Exponential Growth Rate (α)

Controls how amplitude increases along body:

```
α = 0.5  →  Gentle increase  →  More whole-body motion
α = 1.0  →  Default
α = 2.0  →  Rapid increase   →  More tail-concentrated motion
```

**In input2d:** Modify `kinematics_velocity_function { alpha = 1.0 }`

### Amplitude Profiles for Different α

```
α = 0.5:
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.061 0.068 0.075 0.082 0.091 0.100

α = 1.0 (default):
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.037 0.045 0.055 0.067 0.082 0.100

α = 2.0:
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.014 0.020 0.030 0.045 0.067 0.100
```

Higher α concentrates more motion at the tail.

## Output Files

### Force and Power Data

`anguilliform_output/eel_*.dat`: Contains time-series data

Columns:
1. Time (s)
2. Drag force (N)
3. Lateral force (N)
4. Power (W)
5. Center of mass X (m)
6. Center of mass Y (m)
7. Forward velocity (m/s)
8. Lateral velocity (m/s)

### Visualization Files

`viz_anguilliform/`: VisIt database files
- Velocity field
- Pressure field
- Vorticity
- Lagrangian structure

## Expected Results

### Performance Metrics

| Metric | Expected Value | Description |
|--------|----------------|-------------|
| Forward velocity | 0.8-1.2 m/s | 0.8-1.2 BL/s |
| Lateral recoil | 0.03-0.05 m | 3-5% of length |
| Thrust | 0.5-1.0 N | Forward force |
| Power | 1-3 W | Mechanical power |
| Efficiency | 60-65% | Froude efficiency |
| Strouhal number | 0.25-0.35 | St = fA/U |

### Flow Features

1. **Vortex wake**: Alternating vortices shed along entire body
2. **Reverse Kármán street**: Vortices arranged for thrust
3. **Pressure waves**: Traveling along body surface
4. **Boundary layer**: Thin attached flow

### Motion Characteristics

1. **Whole-body undulation**: Visible wave from head to tail
2. **Lateral recoil**: Head oscillates laterally (unavoidable)
3. **Forward progression**: Steady swimming velocity
4. **Smooth motion**: No abrupt changes

## Troubleshooting

### Problem: Simulation crashes early

**Possible causes:**
- Time step too large
- Insufficient grid resolution near body

**Solutions:**
```
// In input2d:
dt_max = 0.0005  // Reduce from 0.001
cfl = 0.2        // Reduce from 0.3
max_levels = 4   // Increase from 3
```

### Problem: Eel doesn't move forward

**Possible causes:**
- Wrong sign in velocity calculation
- Missing reaction force computation
- Boundary conditions too restrictive

**Solutions:**
- Check `velocity[1] = -amplitude * omega * cos(phase);` has negative sign
- Ensure `calculate_translational_momentum = TRUE`
- Make domain larger: `x_lo = -3.0, x_up = 3.0`

### Problem: Unrealistic oscillations

**Possible causes:**
- Amplitude too large
- Frequency too high
- Wrong α parameter

**Solutions:**
```
// Reduce amplitude or frequency:
amplitude = 0.08   // Was 0.1
frequency = 1.5    // Was 2.0
```

### Problem: Slow performance

**Solutions:**
- Reduce simulation time: `end_time = 5.0`
- Reduce output frequency: `viz_dump_interval = 200`
- Use fewer processors efficiently: Try `-np 8` or `-np 16`
- Reduce grid refinement: `max_levels = 2`

## Analysis Scripts

### Plot Swimming Trajectory

```matlab
data = load('anguilliform_output/eel_0000.dat');
time = data(:,1);
x_pos = data(:,5);
y_pos = data(:,6);

figure;
plot(x_pos, y_pos, 'b-', 'LineWidth', 2);
xlabel('X position (m)');
ylabel('Y position (m)');
title('Eel Swimming Trajectory');
axis equal; grid on;
```

### Compute Swimming Efficiency

```matlab
data = load('anguilliform_output/eel_0000.dat');
time = data(:,1);
drag = data(:,2);
power = data(:,4);
velocity = data(:,7);

% Average values (after reaching steady state, t > 2s)
idx = time > 2.0;
avg_velocity = mean(velocity(idx));
avg_power = mean(power(idx));
avg_drag = mean(abs(drag(idx)));

% Froude efficiency
thrust = avg_drag;  % At steady state, thrust ≈ drag
efficiency = (thrust * avg_velocity) / avg_power;

fprintf('Average velocity: %.3f m/s (%.2f BL/s)\n', ...
        avg_velocity, avg_velocity);
fprintf('Average power: %.3f W\n', avg_power);
fprintf('Froude efficiency: %.1f%%\n', efficiency * 100);
```

## References

1. **Tytell, E. D., & Lauder, G. V. (2004).** "The hydrodynamics of eel swimming I. Wake structure." *Journal of Experimental Biology*, 207(11), 1825-1841.

2. **Kern, S., & Koumoutsakos, P. (2006).** "Simulations of optimized anguilliform swimming." *Journal of Experimental Biology*, 209(24), 4841-4857.

3. **Borazjani, I., & Sotiropoulos, F. (2010).** "On the role of form and kinematics on the hydrodynamics of self-propelled body/caudal fin swimming." *Journal of Experimental Biology*, 213(1), 89-107.

## Customization Ideas

### 1. Variable Frequency Swimming
Modify `IBEELKinematics.cpp` to vary frequency over time:
```cpp
double current_frequency = d_frequency * (1.0 + 0.2 * sin(0.5 * d_current_time));
double omega = 2.0 * M_PI * current_frequency;
```

### 2. Burst-and-Coast Swimming
Alternate between active swimming and coasting:
```cpp
// Active for 2s, coast for 1s
double cycle_time = fmod(d_current_time, 3.0);
double active = (cycle_time < 2.0) ? 1.0 : 0.0;
double velocity_y = active * (-amplitude * d_omega * cos(phase));
```

### 3. Different Body Shapes
Modify `generate_eel_geometry.m`:
```matlab
% More tapered tail:
width = max_width * (sin(pi * X).^0.5);

% Bulbous head:
width = max_width * (sin(pi * X).^0.7 + 0.2 * exp(-10*(X-0.1).^2));
```

### 4. C-start Escape Response
Implement rapid body bending for escape:
```cpp
// Rapid C-bend for first 0.1 seconds
if (d_current_time < 0.1) {
    amplitude = 0.3 * d_amplitude;  // 3x normal amplitude
    omega = 4.0 * d_omega;          // 4x frequency
}
```

---

**For more information, see:** `../../docs/Fish_Swimming_Kinematics_Guide.md`

**Next steps:** Try the carangiform (tuna) example in `../carangiform/`
