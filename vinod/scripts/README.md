# Utility Scripts

Collection of scripts to automate common IBAMR tasks.

## Available Scripts

### build.sh
Builds IBAMR examples with CMake.

**Usage**:
```bash
./scripts/build.sh [example_name] [2d|3d] [Debug|Release]
```

**Examples**:
```bash
# Build simple_cylinder in 2D with Release mode (default)
./scripts/build.sh simple_cylinder

# Build in 3D with Debug mode
./scripts/build.sh simple_cylinder 3d Debug

# Build different example
./scripts/build.sh custom_force 2d Release
```

**Features**:
- Automatic dimension detection and configuration
- Parallel build with all available cores
- Color-coded output
- Error handling

---

### run.sh
Runs IBAMR simulations with MPI.

**Usage**:
```bash
./scripts/run.sh [example_name] [num_procs] [input_file]
```

**Examples**:
```bash
# Run with 4 processors (default)
./scripts/run.sh simple_cylinder

# Run with 8 processors
./scripts/run.sh simple_cylinder 8

# Run with custom input file
./scripts/run.sh simple_cylinder 4 my_input.input
```

**Features**:
- Automatic executable detection
- MPI process management
- Output logging to file
- Visualization instructions

---

### analyze.py
Analyzes IBAMR simulation logs and creates plots.

**Usage**:
```bash
python3 scripts/analyze.py <log_file>
```

**Examples**:
```bash
# Analyze default log
python3 scripts/analyze.py IB2d.log

# Analyze custom log
python3 scripts/analyze.py examples/simple_cylinder/simulation.log
```

**Features**:
- Extracts timestep information
- Plots dt and CFL evolution
- Computes statistics (min, max, mean, std)
- Saves plots as PNG files

**Requirements**:
```bash
pip install matplotlib numpy
```

**Output**:
- Console: Summary statistics
- File: `timestep_analysis.png`

---

### clean.sh
Removes build artifacts and simulation output.

**Usage**:
```bash
./scripts/clean.sh [option]
```

**Options**:
- `build` - Clean only build directories and executables
- `output` - Clean only simulation output (viz, restart, logs)
- `all` - Clean everything (default)

**Examples**:
```bash
# Clean everything
./scripts/clean.sh

# Clean only build files
./scripts/clean.sh build

# Clean only output data
./scripts/clean.sh output
```

**What gets cleaned**:
- Build directories: `build_*/`
- Executables: `main2d`, `main3d`
- Object files: `*.o`
- Visualization: `viz_*/`
- Restart files: `restart_*/`
- Post-processing: `postproc_*/`
- Logs: `*.log`
- Plots: `*.png`

---

## Complete Workflow Example

```bash
# 1. Build the example
./scripts/build.sh simple_cylinder 2d Release

# 2. Run the simulation
./scripts/run.sh simple_cylinder 4 input2d

# 3. Analyze results
python3 scripts/analyze.py examples/simple_cylinder/IB2d.log

# 4. Visualize with VisIt
visit -o examples/simple_cylinder/viz_IB2d/dumps.visit

# 5. Clean up when done
./scripts/clean.sh output
```

## Making Scripts Executable

```bash
chmod +x scripts/*.sh
```

## Advanced Usage

### Batch Runs (Parameter Sweep)
Create a script to run multiple simulations:

```bash
#!/bin/bash
for mu in 0.001 0.01 0.1; do
    # Modify input file with new viscosity
    sed "s/MU = .*/MU = $mu/" input2d > input2d_mu$mu

    # Run simulation
    ./scripts/run.sh simple_cylinder 4 input2d_mu$mu

    # Rename output
    mv viz_IB2d viz_IB2d_mu$mu
done
```

### Parallel Builds
Build multiple examples simultaneously:

```bash
#!/bin/bash
for example in simple_cylinder custom_force; do
    ./scripts/build.sh $example 2d Release &
done
wait
echo "All builds complete"
```

### Automated Analysis
Analyze all log files in examples:

```bash
#!/bin/bash
find examples -name "*.log" -type f | while read log; do
    echo "Analyzing: $log"
    python3 scripts/analyze.py "$log"
done
```

## Customization

### Adding New Scripts

1. Create script in `scripts/` directory
2. Add execute permission: `chmod +x scripts/new_script.sh`
3. Follow naming convention: use `.sh` for bash, `.py` for python
4. Add usage instructions as comments at top
5. Document in this README

### Modifying Existing Scripts

All scripts are designed to be easily modified:
- Color definitions at top for changing output colors
- Configuration variables for default values
- Clear function definitions for extending functionality

## Troubleshooting

**Script not found**:
```bash
# Make sure you're in the IBAMR root directory or use full path
cd /path/to/IBAMR
./vinod/scripts/build.sh
```

**Permission denied**:
```bash
chmod +x scripts/*.sh
```

**Python script fails**:
```bash
# Install required packages
pip3 install matplotlib numpy

# Or use conda
conda install matplotlib numpy
```

**Build fails**:
- Check that IBAMR dependencies are installed
- Verify CMake can find IBAMR (set CMAKE_PREFIX_PATH if needed)
- Try Debug build for more information: `./scripts/build.sh example 2d Debug`

## Tips

1. **Use tab completion**: Bash will auto-complete example names
2. **Check logs**: Always review `simulation_output.log` after runs
3. **Clean regularly**: Use `./scripts/clean.sh output` between runs
4. **Save results**: Copy important output before cleaning
5. **Version control**: Keep input files in git, ignore output in `.gitignore`

## References

- [IBAMR Documentation](http://ibamr.github.io)
- [CMake Documentation](https://cmake.org/documentation/)
- [MPI Documentation](https://www.open-mpi.org/doc/)
