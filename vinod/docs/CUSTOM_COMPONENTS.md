# Custom IBAMR Components

This document describes custom C++ classes and components developed for extending IBAMR functionality.

## Custom Classes

### CustomForceFunction
**Location**: `src/CustomForceFunction.h`, `src/CustomForceFunction.cpp`

**Purpose**: Implements custom force application to IB structures

**Key Methods:**
- `setDataOnPatch()` - Applies forces to Lagrangian points
- `computeForce()` - Calculates force based on position and time

**Usage Example:**
```cpp
#include "CustomForceFunction.h"

Pointer<CustomForceFunction> force_fcn = new CustomForceFunction(
    "custom_force",
    app_initializer->getComponentDatabase("CustomForce")
);

ib_method->registerLagrangianForceFunction(force_fcn);
```

### CustomInitializer
**Location**: `src/CustomInitializer.h`, `src/CustomInitializer.cpp`

**Purpose**: Custom initial conditions for fluid and structure

**Key Methods:**
- `initializeFluidVelocity()` - Set initial velocity field
- `initializePressure()` - Set initial pressure field
- `initializeStructure()` - Set initial IB structure configuration

**Usage Example:**
```cpp
#include "CustomInitializer.h"

Pointer<CustomInitializer> init = new CustomInitializer(
    "custom_init",
    grid_geometry,
    input_db->getDatabase("CustomInit")
);
```

## Design Patterns

### Strategy Pattern
Custom components use IBAMR's strategy pattern:
- Inherit from IBAMR strategy base classes
- Override virtual methods for custom behavior
- Register with hierarchy integrator

### Common Base Classes
- `IBLagrangianForceStrategy` - For custom IB forces
- `IBLagrangianSourceStrategy` - For source/sink terms
- `CartGridFunction` - For grid-based initial/boundary conditions
- `LInitStrategy` - For Lagrangian initialization

## Integration with IBAMR

### Registering Custom Components

```cpp
// In main application code
#include "src/CustomForceFunction.h"

// Create and register
Pointer<CustomForceFunction> custom_force = new CustomForceFunction(...);
ib_method->registerLagrangianForceFunction(custom_force);
```

### Building Custom Code

Add to your `Makefile` or `CMakeLists.txt`:
```make
SOURCES += $(VINOD_DIR)/src/CustomForceFunction.cpp
SOURCES += $(VINOD_DIR)/src/CustomInitializer.cpp
```

## Best Practices

1. **Naming**: Use descriptive names with "Custom" prefix to avoid conflicts
2. **Documentation**: Document all public methods with doxygen comments
3. **Testing**: Create simple test cases in `examples/` for each component
4. **Thread Safety**: Ensure MPI safety for parallel execution
5. **Memory Management**: Use SAMRAI smart pointers (`tbox::Pointer`)

## Future Extensions

Planned custom components:
- [ ] Adaptive time-stepping controller
- [ ] Custom refinement criteria for AMR
- [ ] Specialized output writers
- [ ] Performance profiling utilities
