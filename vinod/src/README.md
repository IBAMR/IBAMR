# Custom Source Files

This directory contains custom C++ implementations extending IBAMR functionality.

## Files

### CustomForceFunction
**Files**: `CustomForceFunction.h`, `CustomForceFunction.cpp`

**Purpose**: Applies custom forces to immersed boundary structures

**Features**:
- Time-dependent force ramping
- Position-dependent forcing (customizable)
- Smooth startup with cosine ramp function

**Usage**:
```cpp
#include "vinod/src/CustomForceFunction.h"

// Create force function
Pointer<CustomForceFunction> force_fcn = new CustomForceFunction(
    "custom_force",
    input_db->getDatabase("CustomForce")
);

// Configure
force_fcn->setForceMagnitude(10.0);
force_fcn->setForceDirection(1.0, 0.0);  // Force in +x direction

// Register with IB method
ib_method->registerLagrangianForceFunction(force_fcn);
```

**Input Database Parameters**:
```
CustomForce {
   force_magnitude = 10.0      // Magnitude of applied force
   force_direction = 1.0, 0.0  // Direction vector
   start_time = 0.0            // When to start applying force
   ramp_time = 1.0             // Time to ramp up to full force
}
```

## Building Custom Components

### With CMake
Add to your `CMakeLists.txt`:
```cmake
# Include vinod source directory
include_directories(${PROJECT_SOURCE_DIR}/vinod/src)

# Add source files
set(SOURCES
    main.cpp
    ${PROJECT_SOURCE_DIR}/vinod/src/CustomForceFunction.cpp
)

add_executable(myapp ${SOURCES})
target_link_libraries(myapp IBAMR::IBAMR2d)
```

### With Makefile
Add to your `Makefile`:
```make
VINOD_DIR = ../../vinod
INCLUDES += -I$(VINOD_DIR)/src
SOURCES += $(VINOD_DIR)/src/CustomForceFunction.cpp
```

## Extending Custom Components

### Adding New Custom Classes

1. **Create header file**: Define class inheriting from appropriate IBAMR base
2. **Create implementation**: Implement required virtual methods
3. **Document**: Add doxygen-style comments
4. **Test**: Create example in `examples/` demonstrating usage
5. **Update**: Add documentation to `docs/CUSTOM_COMPONENTS.md`

### Common Base Classes to Extend

| Base Class | Purpose | Example Use Cases |
|------------|---------|-------------------|
| `IBLagrangianForceStrategy` | Custom IB forces | Actuation, external forces |
| `IBLagrangianSourceStrategy` | Source/sink terms | Mass addition/removal |
| `CartGridFunction` | Eulerian functions | Initial conditions, boundaries |
| `LInitStrategy` | Lagrangian initialization | Complex geometries |
| `IBMethodPostProcessStrategy` | Post-processing | Custom output, analysis |

### Best Practices

1. **Namespace**: Consider using a namespace to avoid conflicts
2. **Smart Pointers**: Use `SAMRAI::tbox::Pointer` for SAMRAI objects
3. **Const Correctness**: Mark methods const where appropriate
4. **MPI Safety**: Ensure all operations are parallel-safe
5. **Memory**: Properly restore PETSc arrays after use
6. **Documentation**: Document all parameters and usage

## Testing

Create unit tests or examples for each custom component:
```bash
cd examples/custom_force
make
mpirun -np 4 ./main2d input2d
```

## Future Components

Planned additions:
- [ ] `CustomBoundaryCondition` - Specialized boundary conditions
- [ ] `CustomOutputWriter` - Domain-specific output formats
- [ ] `AdaptiveTimestepper` - Custom time-stepping logic
- [ ] `CustomRefinementCriteria` - Specialized AMR tagging
