# Carangiform (Tuna-like) Swimming Simulation

## Overview

This example demonstrates **carangiform swimming** using IBAMR's ConstraintIB method. Carangiform swimming is characterized by:

- **Posterior body undulation**: Only the rear portion of the body participates in propulsion
- **Quadratic amplitude envelope**: Amplitude increases quadratically with most motion concentrated at the tail
- **Less than one wave**: Wavelength exceeds body length (λ > L)
- **High efficiency**: Superior propulsive efficiency (75-85%) compared to anguilliform mode
- **Streamlined fusiform body**: Football-shaped body with narrow caudal peduncle and enlarged tail
- **Examples in nature**: Tuna, mackerel, jacks (Carangidae family), sharks, dolphins

Carangiform swimmers are among the fastest and most efficient swimmers in nature, capable of sustained high-speed cruising and trans-oceanic migrations.

## Kinematics

### Amplitude Envelope

The carangiform amplitude envelope uses a **quadratic function** that creates minimal anterior motion and concentrates undulation at the posterior body:

```
A(X) = a₀ + a₁X + a₂X²
```

Where:
- `X = x/L`: Normalized position along body (0 = head, 1 = tail)
- `a₀ = 0.01`: Constant term (small baseline amplitude)
- `a₁ = -0.05`: Linear term (NEGATIVE - reduces anterior body motion)
- `a₂ = 0.14`: Quadratic term (creates strong posterior concentration)

The negative linear coefficient `a₁` is critical: it counteracts the constant term in the anterior region, resulting in a nearly rigid head and trunk with amplitude increasing rapidly toward the tail.

### Complete Swimming Equation

```
y(x,t) = A(x) · sin(kx - ωt)
      = [a₀ + a₁(x/L) + a₂(x/L)²] · sin(2π·x/λ - 2π·f·t)
```

### Default Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Frequency | f | 3.0 Hz | Tail beat frequency (higher than anguilliform) |
| Wavelength | λ | 1.5 L | 150% of body length (longer than body) |
| Constant coeff. | a₀ | 0.01 | Baseline amplitude (1% of length) |
| Linear coeff. | a₁ | -0.05 | Negative to suppress anterior motion |
| Quadratic coeff. | a₂ | 0.14 | Concentrates motion at tail (14%) |
| Body length | L | 1.0 m | Fish length |
| Max body width | w_max | 0.12 m | 12% of length (robust fusiform) |
| Caudal peduncle | - | 0.75 L | Narrow region at 75% body length |

### Amplitude Distribution

```
Position (X):    0.0   0.2   0.4   0.6   0.8   1.0
Amplitude (m): 0.010 0.000 0.006 0.028 0.066 0.120
```

Notice how the amplitude is minimal in the anterior half (even zero at X=0.2) and increases dramatically in the posterior third of the body. This is characteristic of carangiform swimming.

### Comparison with Anguilliform Mode

| Feature | Anguilliform | Carangiform |
|---------|-------------|-------------|
| Amplitude function | Exponential | Quadratic |
| Body participation | Whole body | Posterior 1/3 |
| Wavelength | 0.8 L (< 1) | 1.5 L (> 1) |
| Frequency | 2.0 Hz | 3.0 Hz |
| Swimming speed | 0.8-1.2 BL/s | 1.5-3.0 BL/s |
| Efficiency | 60-65% | 75-85% |
| Head motion | Moderate recoil | Minimal recoil |

## File Structure

```
carangiform/
├── IBEELKinematics_Carangiform.cpp   # C++ kinematics implementation
├── IBEELKinematics.h                   # Header file (shared with anguilliform)
├── generate_tuna_geometry.m             # MATLAB geometry generator
├── input2d_carangiform                  # IBAMR input configuration
└── README.md                           # This file
```

## Usage Instructions

### Step 1: Generate Geometry

Run MATLAB to create the tuna body geometry:

```bash
matlab -batch "generate_tuna_geometry"
```

This creates:
- `eel2d.vertex`: Lagrangian mesh file (256 points, closed contour)
- `tuna_geometry.png`: Comprehensive visualization showing:
  - Body outline with Lagrangian points
  - Width profile with peduncle and caudal fin regions
  - Amplitude envelope comparison
  - Local fineness ratio (streamlining)
  - Cross-sectional area distribution

**Geometry features:**
- Fusiform (football-shaped) body with maximum width at mid-body
- Narrow caudal peduncle at 75% length (characteristic of fast swimmers)
- Enlarged caudal fin region at 90% length
- Aspect ratio ~8.3:1 (length to width)
- Optimized for low drag and high thrust

### Step 2: Set Up IBAMR Project

Copy files to your IBAMR project:

```bash
# Create project directory
mkdir -p ~/ibamr-projects/carangiform
cd ~/ibamr-projects/carangiform

# Copy source files
cp IBEELKinematics_Carangiform.cpp IBEELKinematics.cpp
cp IBEELKinematics.h .
cp input2d_carangiform input2d
cp eel2d.vertex .

# Copy example main program (you need example.cpp from IBAMR examples)
cp $IBAMR_ROOT/examples/ConstraintIB/eel2d/example.cpp .
```

### Step 3: Create CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.10)
project(carangiform)

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

**Note:** If you encounter linking errors related to IBAMR libraries, ensure:
- IBAMR is properly installed and IBAMR_DIR is set
- All dependencies (PETSc, HDF5, SAMRAI) are available
- Use `make VERBOSE=1` to debug compilation issues

### Step 5: Run Simulation

```bash
cd ..  # Back to project root
mpirun -np 4 ./build/main2d input2d
```

**Simulation details:**
- Duration: 10 seconds (30 swimming cycles at 3 Hz)
- Time step: 0.8 ms (smaller than anguilliform due to higher frequency)
- Output interval: Every 100 steps (0.08 seconds)
- Grid: 1000×200 base grid with 3 levels of adaptive refinement
- Domain: 5m × 1m (larger for faster swimming)
- Computational time: ~3-6 hours on 4 cores
- Memory: ~4-8 GB RAM recommended

**Expected console output:**
```
+++++++++++++++++++++++++++++++++++++++++++++++++++
At beginning of timestep # 0
Simulation time is 0
IBHierarchyIntegrator::advanceHierarchy(): time interval = [0,0.0008], dt = 0.0008
...
Writing visualization files...
+++++++++++++++++++++++++++++++++++++++++++++++++++
At beginning of timestep # 100
Simulation time is 0.08
Average forward velocity: 0.234 m/s
Thrust: 0.85 N
Power: 2.3 W
```

### Step 6: Visualize Results

Use VisIt to view results:

```bash
visit -o viz_carangiform/dumps.visit
```

**Recommended visualizations:**

1. **Velocity Magnitude with Streamlines**
   - Add plot: Pseudocolor → velocity_magnitude
   - Add plot: Vector → velocity (subsample 4x for clarity)
   - Observe: High-speed wake jets from tail region

2. **Vorticity (Wake Structure)**
   - Add plot: Pseudocolor → vorticity
   - Use diverging colormap (blue-white-red)
   - Observe: Reverse Kármán vortex street (thrust-producing)
   - Look for: Single or double row of vortices depending on Strouhal number

3. **Pressure Distribution**
   - Add plot: Pseudocolor → pressure
   - Observe: Pressure waves along body
   - Look for: High pressure on concave side, low on convex side

4. **Lagrangian Structure Motion**
   - Add plot: Mesh → eel2d
   - Add plot: Pseudocolor → velocity_magnitude (background)
   - Create animation to see swimming motion
   - Observe: Rigid anterior body, undulating tail

**Animation tips:**
- Frame rate: 10-15 fps for smooth visualization
- Use "Save movie" feature in VisIt
- Focus camera on fish (follow center of mass)

## Parameter Tuning

### Frequency (f)

Controls tail beat frequency and swimming speed:

```
f = 2.0 Hz  →  Slow cruising   →  v ≈ 1.2 m/s  →  More efficient
f = 3.0 Hz  →  Moderate (default) →  v ≈ 2.0 m/s
f = 4.0 Hz  →  Fast burst      →  v ≈ 2.8 m/s  →  Higher power
```

**In input2d:** Modify `kinematics_velocity_function { frequency = 3.0 }`

**Effect on performance:**
- Higher f → faster swimming but lower efficiency
- Optimal frequency depends on body size and Reynolds number
- Real tuna: 1-5 Hz depending on species and swimming mode

**Warning:** For f > 4 Hz, reduce time step to dt_max = 0.0005 for stability.

### Wavelength (λ)

Controls the spatial extent of the traveling wave:

```
λ = 1.2 L  →  0.83 waves on body  →  More undulation
λ = 1.5 L  →  0.67 waves (default) →  Optimal efficiency
λ = 2.0 L  →  0.50 waves          →  Less undulation, more rigid
```

**In input2d:** Modify `kinematics_velocity_function { wavelength = 1.5 }`

**Effect on swimming:**
- Shorter λ → more body bending, lower speed, anguilliform-like
- Longer λ → stiffer body, higher speed, thunniform-like
- Carangiform range: λ = 1.2-2.0 L
- Thunniform (pure tail) mode: λ > 2.5 L

### Quadratic Coefficients (a₀, a₁, a₂)

Control the shape of the amplitude envelope:

#### Constant term (a₀):
```
a₀ = 0.005  →  Lower baseline  →  More concentrated tail motion
a₀ = 0.010  →  Default
a₀ = 0.020  →  Higher baseline →  More whole-body motion
```

#### Linear term (a₁):
```
a₁ = -0.08  →  Very rigid head/trunk
a₁ = -0.05  →  Default (rigid anterior)
a₁ = -0.02  →  Moderate anterior motion
a₁ = 0.00   →  No suppression (more anguilliform)
```

**Critical:** a₁ must be negative for true carangiform motion! Positive values create anguilliform-like whole-body undulation.

#### Quadratic term (a₂):
```
a₂ = 0.10   →  Less tail amplitude  →  Lower thrust
a₂ = 0.14   →  Default
a₂ = 0.18   →  More tail amplitude  →  Higher thrust, possible instability
```

**In input2d:** Modify these lines in `kinematics_velocity_function` block:
```
a0 = 0.01
a1 = -0.05
a2 = 0.14
```

### Amplitude Profiles for Different Parameters

**Effect of a₁ (with a₀=0.01, a₂=0.14):**
```
a₁ = -0.08 (very rigid):
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 -0.006 0.006 0.034 0.078 0.138

a₁ = -0.05 (default):
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 0.000 0.006 0.028 0.066 0.120

a₁ = 0.00 (no suppression):
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 0.016 0.032 0.060 0.099 0.150
```

**Effect of a₂ (with a₀=0.01, a₁=-0.05):**
```
a₂ = 0.10:
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 0.000 0.002 0.016 0.046 0.092

a₂ = 0.14 (default):
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 0.000 0.006 0.028 0.066 0.120

a₂ = 0.18:
X:   0.0   0.2   0.4   0.6   0.8   1.0
A: 0.010 0.002 0.013 0.042 0.087 0.148
```

### Body Geometry Parameters

To modify body shape, edit `generate_tuna_geometry.m`:

```matlab
max_width = 0.12;       % Body robustness (8-15% typical)
peduncle_factor = 0.6;  % Peduncle narrowing (0.5-0.7)
tail_factor = 1.3;      % Caudal fin size (1.2-1.5)
```

**Effect of body width:**
- Wider bodies (15%) → more drag, more inertia, slower acceleration
- Narrower bodies (8%) → less drag, faster swimming, may need stiffer motion

**Reynolds number effects:**
- Re = UL/ν ≈ 2×10⁶ for a 1m tuna at 2 m/s
- High Re → turbulent boundary layer, skin friction drag dominant
- CFD simulation at Re = 10³-10⁴ captures qualitative features

## Output Files

### Force and Power Data

`carangiform_output/tuna_*.dat`: Contains time-series data (ASCII format)

**Columns:**
1. Time (s)
2. Drag force (N) - positive = resistance
3. Lateral force (N) - oscillates with swimming motion
4. Mechanical power (W) - instantaneous power input
5. Center of mass X (m) - forward displacement
6. Center of mass Y (m) - lateral recoil
7. Forward velocity (m/s) - swimming speed
8. Lateral velocity (m/s) - lateral oscillation speed

**File naming:**
- `tuna_0000.dat` - processor 0 (master)
- `tuna_0001.dat` - processor 1
- etc.

For single-body simulations, use `tuna_0000.dat` which contains complete data.

### Visualization Files

`viz_carangiform/`: VisIt database files (HDF5/Silo format)

**Contents:**
- `dumps.visit` - Main database file (open this in VisIt)
- `lag_data.cycle_NNNNNN` - Lagrangian structure data
- `viz_ib2d/` - Eulerian field data (velocity, pressure, vorticity)
- `summary.samrai` - SAMRAI patch hierarchy data

**Fields available:**
- `velocity_magnitude` - Speed field |u| (m/s)
- `velocity_0`, `velocity_1` - Velocity components u, v (m/s)
- `pressure` - Pressure field p (Pa)
- `vorticity` - Vorticity ω = ∂v/∂x - ∂u/∂y (1/s)
- `mesh_eel2d` - Lagrangian structure position

### Restart Files

`restart_carangiform/`: Restart checkpoint files for resuming simulation

To restart from checkpoint at t=5.0s:
```bash
mpirun -np 4 ./build/main2d input2d restart_carangiform/restore.5000
```

## Expected Results

### Performance Metrics

| Metric | Expected Value | Typical Range | Description |
|--------|----------------|---------------|-------------|
| Forward velocity | 2.0 m/s | 1.5-3.0 m/s | 1.5-3.0 BL/s (body lengths per second) |
| Lateral recoil | 0.01-0.02 m | 1-2% of L | Much smaller than anguilliform (3-5%) |
| Swimming speed | 2.0 BL/s | 1.5-3.0 BL/s | Faster than anguilliform (0.8-1.2 BL/s) |
| Thrust force | 1.2 N | 0.8-2.0 N | Higher than anguilliform |
| Mechanical power | 2.5 W | 1.5-4.0 W | Power = Force × Velocity |
| Froude efficiency | 78% | 75-85% | η = Thrust×U / Power |
| Strouhal number | 0.28 | 0.25-0.35 | St = fA/U (optimal: 0.25-0.35) |
| Reynolds number | 2×10³ | 10³-10⁴ | Re = UL/ν (CFD regime) |
| Cost of transport | 0.15 J/m | 0.10-0.20 J/m | Energy per distance |

### Swimming Performance vs. Time

**Startup phase (0-2 seconds):**
- Fish accelerates from rest
- Velocity increases gradually to steady state
- Initial transients in force and power
- Vortex wake develops

**Steady swimming (2-10 seconds):**
- Constant average velocity (oscillates ±5% about mean)
- Periodic force and power (frequency = 2f due to two tail beats per cycle)
- Established wake structure
- Use this region for performance analysis

### Flow Features

1. **Reverse Kármán Vortex Street**
   - Vortices arranged to produce thrust (opposite of bluff body wake)
   - Single or double row depending on tail amplitude and Strouhal number
   - Vortices shed primarily from tail region, not whole body
   - Spacing: λ_vortex ≈ U/f ≈ 0.67 m

2. **Jet Wake**
   - High-velocity jet behind fish (wake velocity > swimming velocity)
   - Indicates efficient momentum transfer to fluid
   - Narrow jet (focused from caudal region)

3. **Pressure Distribution**
   - High pressure on concave (accelerating) side of tail
   - Low pressure on convex (decelerating) side
   - Pressure difference creates thrust force
   - Traveling pressure wave along body

4. **Boundary Layer**
   - Thin attached boundary layer on most of body
   - Possible separation at tail during high-amplitude motion
   - Turbulent at high Re, laminar at low Re

### Motion Characteristics

1. **Posterior Body Undulation**
   - Head and anterior body (~50%) remain nearly rigid
   - Amplitude increases rapidly from mid-body to tail
   - Clear traveling wave visible in posterior 50%
   - Smooth, continuous motion

2. **Minimal Lateral Recoil**
   - Head moves laterally only 1-2% of body length
   - Much less than anguilliform mode (3-5%)
   - More efficient: less energy wasted in lateral motion
   - Nearly straight-line trajectory

3. **High Forward Speed**
   - Typical: 1.5-3.0 body lengths per second
   - 50-150% faster than anguilliform mode
   - Real tuna: 2-15 BL/s depending on activity (cruising vs. burst)

4. **Smooth Acceleration**
   - Gradual ramp-up from rest to steady swimming
   - No jerky motions or instabilities
   - Periodic oscillations superimposed on steady motion

### Comparison with Experimental Data

**Published values for real tuna (Scombridae):**
- Swimming speed: 2-10 BL/s (cruising), up to 20 BL/s (burst)
- Tail beat frequency: 1-5 Hz
- Stride length: U/f ≈ 0.5-1.5 L
- Strouhal number: 0.25-0.35
- Efficiency: 75-90% (among highest in nature)

**CFD simulation typically shows:**
- 10-30% lower speeds (due to simplified geometry and low Re)
- Correct qualitative trends with parameter variations
- Accurate capture of wake structure and vortex patterns
- Good agreement with PIV flow visualization experiments

## Troubleshooting

### Problem: Simulation crashes early (t < 0.5s)

**Possible causes:**
- Time step too large for high-frequency motion
- Insufficient grid resolution near tail
- Amplitude too large causing extreme deformations

**Solutions:**
```
// In input2d:
dt_max = 0.0005              // Reduce from 0.0008
cfl = 0.2                    // Reduce from 0.3
max_levels = 4               // Increase from 3
vorticity_rel_thresh = 0.1, 0.25, 0.5  // Finer threshold

// Or reduce motion amplitude:
a2 = 0.10                    // Reduce from 0.14
frequency = 2.5              // Reduce from 3.0
```

### Problem: Tuna doesn't move forward / swims backward

**Possible causes:**
- Wrong sign in velocity calculation (most common!)
- Wavelength too short (λ < L creates anguilliform mode)
- Missing translational momentum calculation
- Amplitude too small to generate thrust

**Diagnostics:**
```bash
# Check force data:
tail -100 carangiform_output/tuna_0000.dat
# Column 2 should be negative (thrust) on average after t > 2s
# Column 7 should be positive (forward velocity)
```

**Solutions:**
- Check `velocity[1] = -amplitude * omega * cos(phase);` has negative sign in .cpp file
- Ensure `calculate_translational_momentum = TRUE` in input2d
- Increase wavelength: `wavelength = 1.5` (must be > 1.0 for carangiform)
- Increase tail amplitude: `a2 = 0.16`
- Make domain larger to allow flow development: `x_lo = -3.0, x_up = 3.0`

### Problem: Excessive lateral oscillations / unstable trajectory

**Possible causes:**
- Amplitude too large (a₂ > 0.18)
- Frequency too high for body stiffness
- Grid resolution insufficient
- Wrong a₁ coefficient (should be negative!)

**Solutions:**
```
// Reduce motion amplitude:
a2 = 0.12                    // Reduce quadratic term
a1 = -0.08                   // More negative (stiffer anterior)

// Or reduce frequency:
frequency = 2.5              // Reduce from 3.0

// Improve resolution:
max_levels = 4
domain_boxes = [(0,0), (1999,399)]  // Finer base grid
```

### Problem: Very slow swimming (v < 1.0 m/s)

**Possible causes:**
- Wavelength too short (creating inefficient anguilliform mode)
- Tail amplitude too small
- Frequency too low
- Positive a₁ coefficient (wrong sign!)

**Solutions:**
```
// Increase propulsive motion:
wavelength = 1.8             // Increase from 1.5
a2 = 0.16                    // Increase tail amplitude
frequency = 3.5              // Increase frequency

// Check sign:
a1 = -0.05                   // Must be NEGATIVE!
```

### Problem: No visible tail motion in visualization

**Possible causes:**
- Quadratic term a₂ too small
- Viewing wrong time range
- Need to zoom in on tail region

**Solutions:**
- Increase tail amplitude: `a2 = 0.16` or higher
- In VisIt, zoom to region x ∈ [0.3, 0.6], y ∈ [-0.15, 0.15]
- Animate over several cycles (at least 1 second of simulation time)
- Verify motion by plotting amplitude envelope (see Analysis Scripts below)

### Problem: Slow computational performance

**Solutions:**

**Reduce simulation cost:**
```
end_time = 5.0               // Reduce from 10.0 (still ~15 cycles)
viz_dump_interval = 200      // Reduce output frequency
max_levels = 2               // Reduce refinement
domain_boxes = [(0,0), (699,139)]  // Smaller domain
```

**Optimize parallel performance:**
```bash
# Try different processor counts:
mpirun -np 8 ./build/main2d input2d   # More processors
mpirun -np 16 ./build/main2d input2d  # Even more

# Enable load balancing:
LoadBalancer {
   bin_pack_method = "SPATIAL"
   max_workload_factor = 0.8  # Tighter balance
}
```

**Profile performance:**
```
// In input2d, enable detailed timing:
timer_dump_interval = 10
```

Then examine the output log to find bottlenecks.

### Problem: Vorticity field looks noisy

**Causes:**
- Insufficient grid resolution
- Numerical oscillations
- Need more levels of refinement

**Solutions:**
```
max_levels = 4               // More refinement
vorticity_rel_thresh = 0.1, 0.25, 0.5  // Tag more regions

// Use higher-order delta function:
delta_fcn = "IB_6"           // Instead of "IB_4"
```

### Problem: Simulation runs but results seem unphysical

**Debugging checklist:**

1. **Check geometry file:**
   ```bash
   head -20 eel2d.vertex
   # Should show 256 (or similar) and then coordinate pairs
   ```

2. **Verify parameter values:**
   ```bash
   grep frequency input2d
   grep wavelength input2d
   grep "a0\|a1\|a2" input2d
   ```

3. **Examine force data:**
   ```matlab
   data = load('carangiform_output/tuna_0000.dat');
   plot(data(:,1), data(:,7));  % Velocity vs time
   xlabel('Time (s)'); ylabel('Velocity (m/s)');
   ```

4. **Check for NaN or Inf:**
   ```bash
   grep -i "nan\|inf" output_carangiform
   ```

5. **Visualize in VisIt:**
   - Check that fish geometry appears correct
   - Verify motion is smooth (no jumps or gaps)
   - Confirm flow field looks reasonable (no extreme values)

## Analysis Scripts

### Plot Swimming Trajectory and Speed

```matlab
%% Load and analyze carangiform swimming data
data = load('carangiform_output/tuna_0000.dat');
time = data(:,1);
drag = data(:,2);
lateral_force = data(:,3);
power = data(:,4);
x_pos = data(:,5);
y_pos = data(:,6);
u_vel = data(:,7);
v_vel = data(:,8);

% Plot trajectory
figure('Position', [100 100 1200 400]);
subplot(1,3,1);
plot(x_pos, y_pos, 'b-', 'LineWidth', 2);
xlabel('X position (m)', 'FontSize', 12);
ylabel('Y position (m)', 'FontSize', 12);
title('Tuna Swimming Trajectory', 'FontSize', 14, 'FontWeight', 'bold');
axis equal; grid on;

% Plot velocity vs time
subplot(1,3,2);
plot(time, u_vel, 'b-', 'LineWidth', 2);
hold on;
yline(mean(u_vel(time > 2)), 'r--', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Forward velocity (m/s)', 'FontSize', 12);
title('Swimming Speed', 'FontSize', 14, 'FontWeight', 'bold');
legend('Instantaneous', 'Mean (t>2s)', 'Location', 'best');
grid on;

% Plot power vs time
subplot(1,3,3);
plot(time, power, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
title('Mechanical Power', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Print summary statistics
fprintf('\n========================================\n');
fprintf('CARANGIFORM SWIMMING ANALYSIS\n');
fprintf('========================================\n');
idx = time > 2.0;  % Steady-state region
fprintf('Average forward velocity: %.3f m/s (%.2f BL/s)\n', ...
        mean(u_vel(idx)), mean(u_vel(idx)));
fprintf('Lateral recoil amplitude: %.4f m (%.2f%% of L)\n', ...
        std(y_pos(idx)), std(y_pos(idx))*100);
fprintf('Average power: %.3f W\n', mean(power(idx)));
fprintf('Distance traveled: %.2f m (%.1f BL)\n', ...
        x_pos(end) - x_pos(1), x_pos(end) - x_pos(1));
fprintf('========================================\n');
```

### Compute Swimming Efficiency and Strouhal Number

```matlab
%% Swimming efficiency analysis
data = load('carangiform_output/tuna_0000.dat');
time = data(:,1);
drag = -data(:,2);  % Negative of drag = thrust
power = data(:,4);
velocity = data(:,7);

% Steady-state analysis (after t > 2s)
idx = time > 2.0 & time < 9.0;
time_ss = time(idx);
drag_ss = drag(idx);
power_ss = power(idx);
velocity_ss = velocity(idx);

% Time-averaged values
avg_velocity = mean(velocity_ss);
avg_thrust = mean(drag_ss);
avg_power = mean(power_ss);

% Froude efficiency (propulsive efficiency)
% η = (Thrust × Velocity) / Power
froude_efficiency = (avg_thrust * avg_velocity) / avg_power;

% Strouhal number
% St = f × A / U
% where A = peak-to-peak tail amplitude
frequency = 3.0;  % Hz (from input file)
tail_amplitude = 0.12;  % m (a₂ coefficient for X=1)
strouhal = frequency * tail_amplitude / avg_velocity;

% Cost of transport (COT)
% Energy per unit distance
simulation_time = time_ss(end) - time_ss(1);
distance = avg_velocity * simulation_time;
total_energy = avg_power * simulation_time;
COT = total_energy / distance;  % J/m

% Display results
fprintf('\n========================================\n');
fprintf('SWIMMING EFFICIENCY ANALYSIS\n');
fprintf('========================================\n');
fprintf('Average velocity:      %.3f m/s (%.2f BL/s)\n', ...
        avg_velocity, avg_velocity);
fprintf('Average thrust force:  %.3f N\n', avg_thrust);
fprintf('Average power:         %.3f W\n', avg_power);
fprintf('Froude efficiency:     %.1f%%\n', froude_efficiency * 100);
fprintf('Strouhal number:       %.3f (optimal: 0.25-0.35)\n', strouhal);
fprintf('Cost of transport:     %.3f J/m\n', COT);
fprintf('Reynolds number:       %.0f\n', avg_velocity * 1.0 / 1e-6);
fprintf('========================================\n');

% Optimal performance check
if strouhal >= 0.25 && strouhal <= 0.35
    fprintf('✓ Strouhal number in optimal range!\n');
else
    fprintf('⚠ Strouhal number outside optimal range\n');
end

if froude_efficiency > 0.75
    fprintf('✓ High efficiency achieved!\n');
elseif froude_efficiency > 0.65
    fprintf('○ Moderate efficiency\n');
else
    fprintf('⚠ Low efficiency - check parameters\n');
end
```

### Analyze Wake Structure and Vorticity

```matlab
%% Frequency analysis of forces and motion
data = load('carangiform_output/tuna_0000.dat');
time = data(:,1);
lateral_force = data(:,3);
power = data(:,4);

% Focus on steady-state
idx = time > 3.0 & time < 9.0;
time_fft = time(idx);
lateral_fft = lateral_force(idx);
power_fft = power(idx);

% Compute FFT
dt = mean(diff(time_fft));
Fs = 1/dt;  % Sampling frequency
N = length(lateral_fft);

% FFT of lateral force
Y_lat = fft(lateral_fft - mean(lateral_fft));
P_lat = abs(Y_lat/N);
P_lat = P_lat(1:N/2+1);
P_lat(2:end-1) = 2*P_lat(2:end-1);
f = Fs*(0:(N/2))/N;

% FFT of power
Y_pow = fft(power_fft - mean(power_fft));
P_pow = abs(Y_pow/N);
P_pow = P_pow(1:N/2+1);
P_pow(2:end-1) = 2*P_pow(2:end-1);

% Plot frequency spectra
figure('Position', [100 100 1200 500]);

subplot(1,2,1);
plot(f, P_lat, 'b-', 'LineWidth', 2);
xlim([0 10]);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Lateral Force Spectrum', 'FontSize', 14, 'FontWeight', 'bold');
xline(3.0, 'r--', 'f', 'LineWidth', 2);  % Swimming frequency
grid on;

subplot(1,2,2);
plot(f, P_pow, 'r-', 'LineWidth', 2);
xlim([0 10]);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Power Spectrum', 'FontSize', 14, 'FontWeight', 'bold');
xline(6.0, 'r--', '2f', 'LineWidth', 2);  % Expected power frequency = 2f
grid on;

% Find dominant frequencies
[~, idx_lat] = max(P_lat(f > 0.5));
dominant_freq_lat = f(f > 0.5);
dominant_freq_lat = dominant_freq_lat(idx_lat);

fprintf('\nFrequency Analysis:\n');
fprintf('Dominant lateral force frequency: %.2f Hz\n', dominant_freq_lat);
fprintf('Expected swimming frequency: 3.00 Hz\n');
```

### Plot Amplitude Envelope

```matlab
%% Visualize amplitude envelope
X = linspace(0, 1, 100);
a0 = 0.01;
a1 = -0.05;
a2 = 0.14;

% Calculate amplitude
A = a0 + a1*X + a2*X.^2;
A(A < 0) = 0;  % Non-negative

% Plot
figure('Position', [100 100 800 600]);
plot(X, A, 'b-', 'LineWidth', 3);
hold on;
plot(X, a0*ones(size(X)), 'r--', 'LineWidth', 1.5);
plot(X, a0 + a1*X, 'g--', 'LineWidth', 1.5);
plot(X, a2*X.^2, 'm--', 'LineWidth', 1.5);

xlabel('Normalized body position (X)', 'FontSize', 14);
ylabel('Amplitude (m)', 'FontSize', 14);
title('Carangiform Amplitude Envelope: A(X) = a₀ + a₁X + a₂X²', ...
      'FontSize', 16, 'FontWeight', 'bold');
legend('Total A(X)', 'Constant a₀', 'Linear a₁X', 'Quadratic a₂X²', ...
       'Location', 'northwest');
grid on;

% Add annotations
text(0.05, 0.005, 'Rigid head', 'FontSize', 12, 'Color', 'red');
text(0.75, 0.09, 'Undulating tail', 'FontSize', 12, 'Color', 'red');

% Mark key positions
xline(0.5, 'k--', 'Mid-body', 'LineWidth', 1);
xline(0.75, 'k--', 'Peduncle', 'LineWidth', 1);

fprintf('\nAmplitude at key positions:\n');
fprintf('X = 0.0 (head):      A = %.4f m (%.1f%% of L)\n', ...
        a0, a0*100);
fprintf('X = 0.2 (anterior):  A = %.4f m (%.1f%% of L)\n', ...
        a0 + a1*0.2 + a2*0.2^2, (a0 + a1*0.2 + a2*0.2^2)*100);
fprintf('X = 0.5 (mid-body):  A = %.4f m (%.1f%% of L)\n', ...
        a0 + a1*0.5 + a2*0.5^2, (a0 + a1*0.5 + a2*0.5^2)*100);
fprintf('X = 0.75 (peduncle): A = %.4f m (%.1f%% of L)\n', ...
        a0 + a1*0.75 + a2*0.75^2, (a0 + a1*0.75 + a2*0.75^2)*100);
fprintf('X = 1.0 (tail):      A = %.4f m (%.1f%% of L)\n', ...
        a0 + a1 + a2, (a0 + a1 + a2)*100);
```

### Compare Multiple Simulations

```matlab
%% Compare different parameter sets
% Assumes you've run multiple simulations with different parameters

% Load data from different runs
data1 = load('run1_f2.5/carangiform_output/tuna_0000.dat');
data2 = load('run2_f3.0/carangiform_output/tuna_0000.dat');
data3 = load('run3_f3.5/carangiform_output/tuna_0000.dat');

% Extract steady-state velocities
idx1 = data1(:,1) > 2.0;
idx2 = data2(:,1) > 2.0;
idx3 = data3(:,1) > 2.0;

vel1 = mean(data1(idx1, 7));
vel2 = mean(data2(idx2, 7));
vel3 = mean(data3(idx3, 7));

pow1 = mean(data1(idx1, 4));
pow2 = mean(data2(idx2, 4));
pow3 = mean(data3(idx3, 4));

freq = [2.5, 3.0, 3.5];
velocity = [vel1, vel2, vel3];
power = [pow1, pow2, pow3];
efficiency = (velocity .* 1.0) ./ power * 100;  % Assuming thrust ≈ 1N

% Plot comparison
figure('Position', [100 100 1200 400]);

subplot(1,3,1);
plot(freq, velocity, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Swimming speed (m/s)', 'FontSize', 12);
title('Speed vs Frequency', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(1,3,2);
plot(freq, power, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
title('Power vs Frequency', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(1,3,3);
plot(freq, efficiency, 'go-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Efficiency (%)', 'FontSize', 12);
title('Efficiency vs Frequency', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
```

## References

### Primary Literature on Carangiform Swimming

1. **Videler, J. J., & Hess, F. (1984).** "Fast continuous swimming of two pelagic predators, saithe (Pollachius virens) and mackerel (Scomber scombrus): a kinematic analysis." *Journal of Experimental Biology*, 109(1), 209-228.
   - Classic paper establishing carangiform kinematics
   - Detailed measurements of body wave parameters

2. **Donley, J. M., & Dickson, K. A. (2000).** "Swimming kinematics of juvenile kawakawa tuna (Euthynnus affinis) and chub mackerel (Scomber japonicus)." *Journal of Experimental Biology*, 203(20), 3103-3116.
   - Quantitative data on tuna swimming modes
   - Body wave amplitude distributions

3. **Lauder, G. V., & Tytell, E. D. (2006).** "Hydrodynamics of undulatory propulsion." *Fish Physiology*, 23, 425-468.
   - Comprehensive review of swimming modes
   - Comparison of anguilliform vs. carangiform efficiency

### Computational Studies

4. **Borazjani, I., & Sotiropoulos, F. (2008).** "Numerical investigation of the hydrodynamics of carangiform swimming in the transitional and inertial flow regimes." *Journal of Experimental Biology*, 211(10), 1541-1558.
   - CFD simulations of carangiform swimming
   - Reynolds number effects on performance
   - Validation against experimental data

5. **Zhu, Q., Wolfgang, M. J., Yue, D. K. P., & Triantafyllou, M. S. (2002).** "Three-dimensional flow structures and vorticity control in fish-like swimming." *Journal of Fluid Mechanics*, 468, 1-28.
   - 3D computational study of fish swimming
   - Vortex wake analysis
   - Optimal Strouhal number

6. **Liu, H., & Kawachi, K. (1999).** "A numerical study of undulatory swimming." *Journal of Computational Physics*, 155(2), 223-247.
   - Immersed boundary method for fish swimming
   - Parameter studies of amplitude envelope

### Efficiency and Energetics

7. **Anderson, J. M., Streitlien, K., Barrett, D. S., & Triantafyllou, M. S. (1998).** "Oscillating foils of high propulsive efficiency." *Journal of Fluid Mechanics*, 360, 41-72.
   - Optimal Strouhal number (0.25-0.35) for thrust production
   - Relationship between kinematics and efficiency

8. **Eloy, C. (2012).** "Optimal Strouhal number for swimming animals." *Journal of Fluids and Structures*, 30, 205-218.
   - Theoretical analysis of optimal swimming parameters
   - Universal scaling laws

### Biological Context

9. **Sfakiotakis, M., Lane, D. M., & Davies, J. B. C. (1999).** "Review of fish swimming modes for aquatic locomotion." *IEEE Journal of Oceanic Engineering*, 24(2), 237-252.
   - Comprehensive classification of swimming modes
   - Engineering applications

10. **Webb, P. W. (1975).** "Hydrodynamics and energetics of fish propulsion." *Bulletin of the Fisheries Research Board of Canada*, 190, 1-159.
    - Foundational work on fish swimming mechanics
    - Body and caudal fin (BCF) propulsion modes

### IBAMR and Immersed Boundary Method

11. **Griffith, B. E., & Patankar, N. A. (2020).** "Immersed methods for fluid-structure interaction." *Annual Review of Fluid Mechanics*, 52, 421-448.
    - Review of IB methods including ConstraintIB
    - Applications to biofluid mechanics

12. **Bhalla, A. P. S., Bale, R., Griffith, B. E., & Patankar, N. A. (2013).** "A unified mathematical framework and an adaptive numerical method for fluid-structure interaction with rigid, deforming, and elastic bodies." *Journal of Computational Physics*, 250, 446-476.
    - ConstraintIB method formulation
    - Rigid body kinematics

## Customization Ideas

### 1. Variable Frequency Swimming (Gait Transition)

Modify `IBEELKinematics_Carangiform.cpp` to vary frequency over time:

```cpp
// In getVelocity() function, replace fixed frequency with time-varying:
// Gradual acceleration: increase frequency from 2 Hz to 4 Hz
double current_frequency = 2.0 + 2.0 * (d_current_time / 5.0);
if (current_frequency > 4.0) current_frequency = 4.0;
double omega = 2.0 * M_PI * current_frequency;
```

**Application:** Simulate fish accelerating from cruising to burst swimming.

### 2. Transition Between Swimming Modes

Blend carangiform and anguilliform kinematics:

```cpp
// Hybrid amplitude envelope
// α = 0: pure carangiform (quadratic)
// α = 1: pure anguilliform (exponential)
double alpha_blend = 0.5;  // 50-50 blend
double A_carang = d_a0 + d_a1 * X + d_a2 * X * X;
double A_anguil = 0.1 * exp(1.0 * (X - 1.0));
double amplitude = (1.0 - alpha_blend) * A_carang + alpha_blend * A_anguil;
```

**Application:** Explore intermediate swimming modes (sub-carangiform).

### 3. Burst-and-Glide Swimming

Alternate between active swimming and passive coasting:

```cpp
// In getVelocity():
double cycle_period = 4.0;  // 4-second cycle
double active_fraction = 0.6;  // 60% active, 40% glide
double cycle_time = fmod(d_current_time, cycle_period);
double active = (cycle_time < cycle_period * active_fraction) ? 1.0 : 0.0;

// Apply to velocity
double velocity_y = active * (-amplitude * d_omega * cos(phase));
```

**Application:** Energy-saving swimming strategy, reduce cost of transport by 20-30%.

### 4. C-Start Escape Response

Implement rapid body bending for escape maneuvers:

```cpp
// In getVelocity():
if (d_current_time < 0.15) {
    // Stage 1 (0-0.05s): Rapid C-bend
    amplitude = 0.25;  // 25% of body length!
    omega = 8.0 * d_omega;  // 8x frequency

} else if (d_current_time < 0.30) {
    // Stage 2 (0.15-0.30s): Propulsive stroke
    amplitude = 0.18;
    omega = 6.0 * d_omega;

} else {
    // Normal swimming resumes
    // ... standard kinematics
}
```

**Application:** Predator avoidance, maximum acceleration (up to 40 m/s²).

### 5. Thunniform Swimming (Pure Tail Mode)

Modify for even more concentrated motion (tuna at high speed):

```cpp
// In input2d:
wavelength = 2.5      // Much longer than body
a0 = 0.005            // Minimal baseline
a1 = -0.10            // Very rigid anterior
a2 = 0.20             // Extreme tail concentration
frequency = 4.0       // High frequency
```

**Application:** Simulate high-speed cruising mode of tuna and billfish.

### 6. Different Body Shapes

Modify `generate_tuna_geometry.m` for different species:

```matlab
%% Mackerel (slender carangiform)
max_width = 0.08;       % Slender (8%)
peduncle_factor = 0.7;  % Less pronounced
tail_factor = 1.2;      % Smaller tail

%% Trevally/Jack (robust carangiform)
max_width = 0.15;       % Robust (15%)
peduncle_factor = 0.5;  % Very narrow peduncle
tail_factor = 1.5;      % Large forked tail

%% Shark (sub-carangiform to thunniform)
max_width = 0.10;       % Moderate (10%)
peduncle_factor = 0.6;
tail_factor = 1.4;
% Add heterocercal tail asymmetry
```

### 7. Schooling Behavior

Run multiple fish with phase offsets:

```cpp
// In input2d, define multiple structures:
ConstraintIBKinematics {
   fish1 {
      structure_names = "fish1"
      kinematics_velocity_function {
         frequency = 3.0
         phase_offset = 0.0  // Leader
         initCenterOfMass = 0.0, 0.0
      }
   }

   fish2 {
      structure_names = "fish2"
      kinematics_velocity_function {
         frequency = 3.0
         phase_offset = 0.5  // Half cycle behind
         initCenterOfMass = -0.5, 0.3  // Offset position
      }
   }
}
```

**Application:** Study hydrodynamic benefits of schooling (energy savings up to 20%).

### 8. Optimize for Maximum Efficiency

Use parameter sweeps to find optimal configuration:

```bash
# Run parameter sweep script
for freq in 2.0 2.5 3.0 3.5 4.0; do
  for a2 in 0.10 0.12 0.14 0.16 0.18; do
    # Modify input2d
    sed -i "s/frequency = .*/frequency = $freq/" input2d
    sed -i "s/a2 = .*/a2 = $a2/" input2d

    # Run simulation
    mpirun -np 4 ./build/main2d input2d

    # Save results
    mv carangiform_output results_f${freq}_a${a2}/
  done
done

# Analyze all results
matlab -batch "analyze_parameter_sweep"
```

### 9. Add Pectoral Fins

Extend kinematics to include pectoral fin motion:

```cpp
// Add pectoral fin oscillation (for maneuvering)
if (isPectoralFin(point)) {
    double fin_frequency = 2.0 * d_frequency;  // Higher frequency
    double fin_amplitude = 0.05;  // Small amplitude
    double fin_phase = d_omega * d_current_time;

    velocity[0] = fin_amplitude * omega * sin(fin_phase);  // Fore-aft
    velocity[1] = -fin_amplitude * omega * cos(fin_phase); // Dorso-ventral
}
```

**Application:** Study maneuvering and station-keeping.

### 10. Adjust for Different Reynolds Numbers

Scale simulation for different fish sizes:

```cpp
// In input2d:
// Small mackerel (L=0.3m, U=0.6m/s, Re=1.8×10⁵)
mu = 1.0e-6 * 1000 / (0.6 * 0.3)  // Adjust viscosity

// Large tuna (L=2.0m, U=4.0m/s, Re=8×10⁶)
mu = 1.0e-6 * 1000 / (4.0 * 2.0)
```

Or use dimensionless formulation for direct Re control.

---

## Additional Resources

**IBAMR Documentation:**
- Main website: https://ibamr.github.io/
- ConstraintIB examples: `$IBAMR_ROOT/examples/ConstraintIB/`
- User guide and tutorials

**VisIt Visualization:**
- VisIt website: https://visit-dav.github.io/visit-website/
- Python scripting for automated analysis
- Movie generation and post-processing

**Project Documentation:**
- Main README: `../../README.md`
- Anguilliform example: `../anguilliform/README.md`
- Kinematics guide: `../../docs/Fish_Swimming_Kinematics_Guide.md`

---

**Questions or Issues?**

If you encounter problems or have questions about this example:
1. Check the Troubleshooting section above
2. Review IBAMR documentation and examples
3. Examine the source code comments in `.cpp` files
4. Compare with the anguilliform example for guidance

**Next Steps:**
- Experiment with different parameter combinations
- Analyze efficiency vs. speed trade-offs
- Compare anguilliform and carangiform modes
- Extend to 3D simulations
- Implement custom swimming gaits

**Good luck with your carangiform swimming simulations!**
