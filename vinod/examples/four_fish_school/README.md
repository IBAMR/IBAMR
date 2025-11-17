# Four Fish School Simulation

This example demonstrates fish schooling behavior using IBAMR's immersed boundary method.

## Setup Instructions

### Step 1: Copy Your Code Here

On your local machine, copy files from your Four_fish_school repository:

```bash
# Navigate to your local Four_fish_school repository
cd /path/to/Four_fish_school

# Copy all relevant files to the IBAMR vinod workspace
cp *.cpp /path/to/IBAMR/vinod/examples/four_fish_school/
cp *.h /path/to/IBAMR/vinod/examples/four_fish_school/
cp *.input /path/to/IBAMR/vinod/examples/four_fish_school/
cp *.vertex /path/to/IBAMR/vinod/examples/four_fish_school/
cp Makefile* /path/to/IBAMR/vinod/examples/four_fish_school/
cp CMakeLists.txt /path/to/IBAMR/vinod/examples/four_fish_school/ 2>/dev/null

# Copy any additional geometry or data files
cp *.spring /path/to/IBAMR/vinod/examples/four_fish_school/ 2>/dev/null
cp *.beam /path/to/IBAMR/vinod/examples/four_fish_school/ 2>/dev/null
```

### Step 2: Review and Organize

Once files are copied, organize them:

```
four_fish_school/
├── main.cpp              # Main simulation code
├── FishKinematics.h      # Custom fish motion classes (if any)
├── FishKinematics.cpp
├── input2d               # 2D simulation parameters
├── input3d               # 3D simulation parameters (if applicable)
├── fish.vertex           # Fish geometry
├── fish.spring           # Spring connections (if applicable)
├── CMakeLists.txt        # Build configuration
└── README.md             # This file
```

### Step 3: Integration Checklist

After copying files, verify:

- [ ] All source files (.cpp, .h) are present
- [ ] Input files (.input) are complete
- [ ] Geometry files (.vertex, .spring, etc.) are included
- [ ] Build files (CMakeLists.txt or Makefile) are present
- [ ] No absolute paths in code (should use relative paths)
- [ ] Dependencies are documented

## Expected Files from Four_fish_school

Based on typical IBAMR fish schooling simulations, you should have:

### Core Files
- **main.cpp** - Main simulation driver
- **input2d** or **input3d** - Parameter configuration
- **fish geometry files** - .vertex, .spring, .beam files defining fish shape

### Possible Custom Classes
- **FishKinematics** - Swimming motion prescription
- **FishForces** - Force generation for swimming
- **SchoolingController** - Coordination between fish
- **CustomOutput** - Specialized data output

### Configuration
- Physical parameters (fluid viscosity, density)
- Fish geometry parameters
- Swimming kinematics parameters
- Grid and AMR settings

## Building

Once files are in place:

```bash
# Using the build script
cd /path/to/IBAMR
./vinod/scripts/build.sh four_fish_school 2d Release

# Or manually
cd vinod/examples/four_fish_school
mkdir build && cd build
cmake ..
make -j
```

## Running

```bash
# Using the run script
./vinod/scripts/run.sh four_fish_school 4 input2d

# Or manually
cd vinod/examples/four_fish_school
mpirun -np 4 ./main2d input2d
```

## Analysis

After running, analyze results:

```bash
python3 vinod/scripts/analyze.py vinod/examples/four_fish_school/*.log
```

## Next Steps

After copying your code here:
1. I can review it against IBAMR best practices
2. Optimize performance and code structure
3. Add documentation and comments
4. Create additional analysis tools
5. Extend functionality as needed

---

**Status**: Awaiting code files from Four_fish_school repository
