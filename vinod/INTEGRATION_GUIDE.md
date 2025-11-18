# Integration Guide: Adding Your Code to Vinod Workspace

This guide explains how to integrate your existing IBAMR projects into the vinod workspace.

## Quick Integration Steps

### 1. Copy Your Repository Files

On your **local machine** (not here), run these commands:

```bash
# Navigate to your IBAMR fork
cd /path/to/your/IBAMR

# Clone your Four_fish_school repository locally
cd /tmp
git clone https://github.com/vinodthale/Four_fish_school.git

# Copy files to vinod workspace
cp -r Four_fish_school/* /path/to/your/IBAMR/vinod/examples/four_fish_school/

# Return to IBAMR directory
cd /path/to/your/IBAMR

# Add and commit
git add vinod/examples/four_fish_school/
git commit -m "Add four fish school simulation to vinod workspace"
git push origin claude/understand-codebase-01J178x5nJhA26M5sbtXqxea
```

### 2. Alternative: Use Git Subtree

If you want to maintain connection to the original Four_fish_school repo:

```bash
cd /path/to/your/IBAMR

# Add Four_fish_school as a subtree
git subtree add --prefix=vinod/examples/four_fish_school \
    https://github.com/vinodthale/Four_fish_school.git main --squash

# Later, to pull updates from Four_fish_school
git subtree pull --prefix=vinod/examples/four_fish_school \
    https://github.com/vinodthale/Four_fish_school.git main --squash
```

### 3. Alternative: Use Git Submodule

For a lightweight reference:

```bash
cd /path/to/your/IBAMR

# Add as submodule
git submodule add https://github.com/vinodthale/Four_fish_school.git \
    vinod/examples/four_fish_school

# Initialize and update
git submodule init
git submodule update

# Commit the submodule reference
git commit -m "Add Four_fish_school as submodule"
git push
```

## Integration for Different Project Types

### Scenario A: Complete IBAMR Application

If Four_fish_school is a complete IBAMR application:

```
vinod/examples/four_fish_school/
├── main.cpp
├── input2d
├── geometry files
└── CMakeLists.txt
```

**Action**: Copy as-is to `vinod/examples/four_fish_school/`

### Scenario B: Custom Classes Only

If it contains reusable custom classes:

```
Four_fish_school/
├── FishKinematics.h
├── FishKinematics.cpp
└── SchoolingController.h
```

**Action**:
1. Copy `.h` files to `vinod/src/`
2. Copy `.cpp` files to `vinod/src/`
3. Update `vinod/src/README.md` with documentation

### Scenario C: Mix of Examples and Classes

If it has multiple examples and shared code:

```
vinod/
├── src/
│   ├── FishKinematics.h       # Shared classes
│   └── FishKinematics.cpp
└── examples/
    ├── single_fish/           # Example 1
    ├── two_fish/              # Example 2
    └── four_fish_school/      # Example 3
```

## Post-Integration Checklist

After copying files, verify:

### Code Quality
- [ ] Code compiles without errors
- [ ] No hardcoded absolute paths
- [ ] Proper IBAMR header includes
- [ ] Follows IBAMR coding conventions

### Documentation
- [ ] README.md in example directory
- [ ] Comments explain key algorithms
- [ ] Input file parameters documented
- [ ] Usage instructions provided

### Build System
- [ ] CMakeLists.txt or Makefile present
- [ ] Can build with vinod/scripts/build.sh
- [ ] Dependencies clearly specified

### Testing
- [ ] Example runs without crashing
- [ ] Output is reasonable
- [ ] Can visualize results with VisIt

## Common Integration Issues

### Issue 1: Absolute Paths

**Problem**: Code has hardcoded paths like `/home/old_user/project/`

**Solution**: Replace with relative paths
```cpp
// Bad
string vertex_file = "/home/old_user/Four_fish_school/fish.vertex";

// Good
string vertex_file = "fish.vertex";
```

### Issue 2: Missing Dependencies

**Problem**: Code relies on custom libraries not in IBAMR

**Solution**:
1. Document in README
2. Add to CMakeLists.txt
3. Consider providing fallback implementation

### Issue 3: Different IBAMR Version

**Problem**: Code written for older/newer IBAMR version

**Solution**:
1. Check IBAMR version compatibility
2. Update deprecated function calls
3. Test with current IBAMR version (0.19.0-pre)

### Issue 4: Build System Differences

**Problem**: Original uses autotools, vinod uses CMake

**Solution**: Create new CMakeLists.txt:
```cmake
cmake_minimum_required(VERSION 3.15.0)
project(four_fish_school CXX)

find_package(IBAMR REQUIRED)

add_executable(main2d main.cpp)
target_link_libraries(main2d PRIVATE IBAMR::IBAMR2d)
```

## Review and Optimization

Once integrated, I can help with:

### Code Review
- IBAMR best practices compliance
- Memory management and smart pointers
- MPI safety and parallel efficiency
- Error handling

### Performance Optimization
- Profiling and bottleneck identification
- AMR tuning for fish schooling
- Load balancing optimization
- I/O optimization

### Feature Extensions
- Additional fish behaviors
- Different schooling patterns
- Interactive control
- Advanced analysis tools

### Documentation
- Code comments and explanations
- Theory and algorithm documentation
- Usage tutorials
- Parameter tuning guides

## Example Integration Session

Here's what we'll do once you copy the files:

```bash
# 1. Copy your code locally
cp -r /path/to/Four_fish_school/* vinod/examples/four_fish_school/

# 2. Push to this session (if using Claude Code locally)
# Files will appear here automatically

# 3. I will:
#    - Review all files
#    - Check IBAMR compatibility
#    - Suggest improvements
#    - Add documentation
#    - Optimize as needed
#    - Commit changes

# 4. You can then:
#    - Test the integrated code
#    - Use vinod scripts for automation
#    - Extend functionality
#    - Share improvements back
```

## Next Steps

**Choose your integration method:**

1. **Manual Copy** (simplest)
   - Copy files directly
   - No git connection to original repo
   - Full control

2. **Git Subtree** (recommended)
   - Maintains history
   - Can pull updates
   - Clean integration

3. **Git Submodule** (for reference)
   - Lightweight
   - Points to original repo
   - Requires submodule commands

**Then:**
1. Copy your Four_fish_school code using chosen method
2. Let me know when ready
3. I'll review and integrate properly
4. We'll test and optimize together

---

**Ready to integrate?** Just copy your Four_fish_school files to the vinod workspace and let me know!
