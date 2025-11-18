# Development Notes

Research notes and ideas for IBAMR development.

## Current Research Focus

### Fluid-Structure Interaction Scenarios
- Flexible filaments in flow
- Multiple interacting structures
- Complex geometries with IB method

### Performance Optimization
- Load balancing strategies for irregular geometries
- Efficient Lagrangian-Eulerian coupling
- AMR refinement criteria optimization

## Ideas and TODO

### Short Term
- [ ] Implement parameter sweep utility
- [ ] Create visualization scripts for common quantities
- [ ] Add restart file analysis tools
- [ ] Document best practices for mesh convergence studies

### Long Term
- [ ] Investigate alternative IB kernels
- [ ] Explore machine learning for adaptive refinement
- [ ] Develop domain-specific examples (biomedical, aerospace)
- [ ] Create automated testing framework

## Lessons Learned

### Timestep Selection
- CFL condition is critical for stability
- IB method adds additional stability constraints
- Start with smaller `DT` and gradually increase

### Grid Resolution
- Rule of thumb: ~4-6 grid points per IB kernel width
- Lagrangian spacing should match or be finer than Eulerian grid
- Use AMR to concentrate resolution near structures

### Parallel Scaling
- Load balancing becomes important for complex geometries
- Consider workload when partitioning Lagrangian points
- Profile to identify bottlenecks

## Useful Commands

### Build and Run
```bash
# Configure with CMake
cmake -DCMAKE_BUILD_TYPE=Release -DIBAMR_DIMENSIONS=2 ..
make -j8

# Run with profiling
mpirun -np 4 ./main2d input2d 2>&1 | tee output.log
```

### Analysis
```bash
# Extract timing data
grep "Time to" output.log

# Check conservation properties
grep "Volume" output.log | awk '{print $3}' | python3 -m matplotlib.pylab
```

### Debugging
```bash
# Run with debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..
mpirun -np 1 gdb --args ./main2d input2d

# Check memory leaks
mpirun -np 1 valgrind --leak-check=full ./main2d input2d
```

## References

### Key Papers
1. Peskin (2002) - "The immersed boundary method" - Foundational IB method paper
2. Griffith & Hornung (2012) - Original IBAMR paper
3. Griffith et al. (2007) - Adaptive IB method with structured AMR

### Useful Resources
- IBAMR Documentation: http://ibamr.github.io
- IBAMR Users Group: ibamr-users@googlegroups.com
- SAMRAI Documentation: https://computing.llnl.gov/projects/samrai
- PETSc Documentation: https://petsc.org/release/

## Meeting Notes

### [Date: Add dates as you work]

**Topics Discussed:**
-

**Action Items:**
-

**Results:**
-

---

*Last Updated: 2025-11-17*
