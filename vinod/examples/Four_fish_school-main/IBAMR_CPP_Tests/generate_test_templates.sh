#!/bin/bash
# Script to generate template files for all IBAMR C++ tests

# Test definitions: TEST_NUM|TEST_DIR|TEST_NAME|DESCRIPTION
declare -a TESTS=(
"02|Test02_Diffusion_Analytic|Pure Diffusion Analytic|Validates diffusion operator against analytical Gaussian solution"
"03|Test03_Advection_Analytic|Pure Advection|Tests advection operator with moving Gaussian"
"04|Test04_MMS|Method of Manufactured Solutions|Verifies combined advection-diffusion with manufactured solution"
"05|Test05_Discontinuous|Discontinuous Initial Condition|Tests stability with top-hat function"
"06|Test06_MassConservation|Mass Conservation|Validates global mass budget over time"
"07|Test07_BCs|Boundary Conditions|Tests Dirichlet, Neumann, and Robin boundary conditions"
"08|Test08_SphereSource|Sphere Source|Compares to Lei et al. sphere source data"
"09|Test09_HighSc|High Schmidt Number|Tests stability at Sc = 100-1000"
"10|Test10_MovingIB|Moving Immersed Boundary|Validates IB-scalar coupling with moving ellipsoid"
"11|Test11_AMR|AMR Sensitivity|Tests for adaptive mesh refinement artifacts"
"12|Test12_TimeStep|Time-step Convergence|Validates temporal convergence and CFL effects"
"13|Test13_LongRun|Long Run Stability|Tests long-term stability and conservation"
"14|Test14_Benchmarks|Benchmark Comparison|Compares to Lei et al. and Kamran et al. data"
)

for test_info in "${TESTS[@]}"; do
    IFS='|' read -r NUM DIR NAME DESC <<< "$test_info"

    echo "Creating test ${NUM}: ${NAME}..."

    # Create main.cpp template
    cat > ${DIR}/main.cpp << 'EOFMAIN'
// =============================================================================
// TEST __NUM__: __NAME__
// =============================================================================
// __DESC__
// =============================================================================

#include <SAMRAI_config.h>
#include <petscsys.h>

// SAMRAI headers
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// IBAMR headers
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

// IBTK headers
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Common test utilities
#include "TestUtilities.h"
#include "ErrorCalculator.h"
#include "AnalyticalSolutions.h"

#include <ibamr/app_namespaces.h>

using namespace TestUtilities;

int main(int argc, char* argv[])
{
    // Initialize IBAMR/SAMRAI/PETSc
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TestUtils::printTestHeader("__NAME__", __NUM__);

    Timer total_timer;
    total_timer.start();

    bool test_passed = true;
    ResultLogger logger("test__NUM2___results.txt");

    {
        // Initialize application
        TestUtils::printProgress("Initializing application...");
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "test__NUM2__.log");

        // TODO: Implement test-specific logic here
        // See Test01 for reference implementation

        TestUtils::printProgress("Test implementation pending...");
        test_passed = false;  // Mark as incomplete

    }

    total_timer.stop();
    logger.logTestResult("Test __NUM__: __NAME__", test_passed);

    TestUtils::printTestFooter("Test __NUM__: __NAME__", test_passed,
                              test_passed ? "Test completed" : "Test not yet implemented");

    return test_passed ? 0 : 1;
}
EOFMAIN

    # Replace placeholders
    sed -i "s/__NUM__/${NUM}/g" ${DIR}/main.cpp
    sed -i "s/__NUM2__/$(printf "%02d" $NUM)/g" ${DIR}/main.cpp
    sed -i "s/__NAME__/${NAME}/g" ${DIR}/main.cpp
    sed -i "s/__DESC__/${DESC}/g" ${DIR}/main.cpp

    # Create README template
    cat > ${DIR}/README.md << EOFREADME
# Test ${NUM}: ${NAME}

## Purpose

${DESC}

## Test Description

**Status**: Template created - implementation pending

This test validates: (TODO: Add specific validation criteria)

## Pass/Fail Criteria

The test PASSES if:
- [ ] (TODO: Add specific criteria)

## Building and Running

\`\`\`bash
cd ${DIR}
mkdir build && cd build
cmake ..
make
./test${NUM}_* ../input2d
\`\`\`

## Expected Output

(TODO: Describe expected output)

## Implementation Notes

(TODO: Add implementation details based on Test01 example)

---
**Test Status**: Template
**Last Updated**: 2025-11-17
EOFREADME

    # Create basic input2d file
    cat > ${DIR}/input2d << EOFINPUT
// Test ${NUM}: ${NAME}
// ${DESC}

diffusion_coefficient = 0.001

CartesianGeometry {
   domain_boxes = [ (0,0) , (63,63) ]
   x_lo = -1.0, -1.0
   x_up = 1.0, 1.0
   periodic_dimension = 0, 0
}

Main {
   log_file_name = "test$(printf "%02d" $NUM).log"
   log_all_nodes = FALSE
   viz_writer = "VisIt"
   viz_dump_interval = 10
   viz_dump_dirname = "viz_test$(printf "%02d" $NUM)"
   timer_enabled = TRUE

   DT = 0.01
   START_TIME = 0.0
   END_TIME = 1.0
}

// TODO: Add test-specific configuration

INSStaggeredHierarchyIntegrator {
   mu = 0.01
   rho = 1.0
   start_time = 0.0
   end_time = 1.0
   num_cycles = 1
}

AdvDiffHierarchyIntegrator {
   start_time = 0.0
   end_time = 1.0
   num_cycles = 1
   diffusion_time_stepping_type = "CRANK_NICOLSON"
}

GriddingAlgorithm {
   max_levels = 1
}

StandardTagAndInitialize {
   tagging_method = "REFINE_BOXES"
}

LoadBalancer {
   bin_pack_method = "SPATIAL"
}

OdorInitialConditions {
   function = "exp(-10.0*(X_0*X_0 + X_1*X_1))"
}

OdorBcCoefs {
   acoef_function_0 = "1.0"
   acoef_function_1 = "1.0"
   acoef_function_2 = "1.0"
   acoef_function_3 = "1.0"

   bcoef_function_0 = "0.0"
   bcoef_function_1 = "0.0"
   bcoef_function_2 = "0.0"
   bcoef_function_3 = "0.0"

   gcoef_function_0 = "0.0"
   gcoef_function_1 = "0.0"
   gcoef_function_2 = "0.0"
   gcoef_function_3 = "0.0"
}
EOFINPUT

    echo "Created ${DIR}/main.cpp, README.md, and input2d"
done

echo ""
echo "All test templates created successfully!"
