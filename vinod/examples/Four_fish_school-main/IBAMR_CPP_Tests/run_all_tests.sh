#!/bin/bash
# Script to run all IBAMR C++ tests sequentially

set +e  # Don't exit on error (some tests may fail)

echo "=========================================="
echo "Running IBAMR C++ Test Suite"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Check if build directory exists
if [ ! -d "build" ]; then
    echo "ERROR: Build directory not found"
    echo "Please run build_all_tests.sh first"
    exit 1
fi

# Test definitions
declare -a TESTS=(
"01:test01_smoke:Test01_SmokeTest:Smoke Test"
"02:test02_diffusion:Test02_Diffusion_Analytic:Pure Diffusion"
"03:test03_advection:Test03_Advection_Analytic:Pure Advection"
"04:test04_mms:Test04_MMS:MMS"
"05:test05_discontinuous:Test05_Discontinuous:Discontinuous"
"06:test06_mass_conservation:Test06_MassConservation:Mass Conservation"
"07:test07_bcs:Test07_BCs:Boundary Conditions"
"08:test08_sphere_source:Test08_SphereSource:Sphere Source"
"09:test09_high_sc:Test09_HighSc:High Schmidt"
"10:test10_moving_ib:Test10_MovingIB:Moving IB"
"11:test11_amr:Test11_AMR:AMR"
"12:test12_timestep:Test12_TimeStep:Time-step"
"13:test13_long_run:Test13_LongRun:Long Run"
"14:test14_benchmarks:Test14_Benchmarks:Benchmarks"
)

# Results tracking
declare -a RESULTS
PASSED=0
FAILED=0
SKIPPED=0

# Run each test
for test_info in "${TESTS[@]}"; do
    IFS=':' read -r NUM EXEC DIR NAME <<< "$test_info"

    echo "=========================================="
    echo "Test ${NUM}: ${NAME}"
    echo "=========================================="

    # Check if executable exists
    if [ ! -f "build/${EXEC}" ]; then
        echo "SKIPPED: Executable not found (build/${EXEC})"
        RESULTS+=("${NUM}:SKIP:${NAME}")
        ((SKIPPED++))
        echo ""
        continue
    fi

    # Check if directory exists
    if [ ! -d "${DIR}" ]; then
        echo "SKIPPED: Test directory not found (${DIR})"
        RESULTS+=("${NUM}:SKIP:${NAME}")
        ((SKIPPED++))
        echo ""
        continue
    fi

    # Run test
    cd ${DIR}
    START_TIME=$(date +%s)

    ../build/${EXEC} input2d > test${NUM}_output.log 2>&1
    EXIT_CODE=$?

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    cd ..

    # Check result
    if [ $EXIT_CODE -eq 0 ]; then
        echo "PASSED (${ELAPSED}s)"
        RESULTS+=("${NUM}:PASS:${NAME}:${ELAPSED}s")
        ((PASSED++))
    else
        echo "FAILED (exit code ${EXIT_CODE})"
        RESULTS+=("${NUM}:FAIL:${NAME}:Exit code ${EXIT_CODE}")
        ((FAILED++))
    fi

    echo ""
done

# Print summary
echo "=========================================="
echo "TEST SUITE SUMMARY"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Results:"
printf "%-6s %-8s %-30s %s\n" "Test" "Status" "Name" "Notes"
echo "------------------------------------------------------------------"

for result in "${RESULTS[@]}"; do
    IFS=':' read -r NUM STATUS NAME NOTES <<< "$result"
    printf "%-6s %-8s %-30s %s\n" "#${NUM}" "${STATUS}" "${NAME}" "${NOTES}"
done

echo "------------------------------------------------------------------"
echo "Total: $((PASSED + FAILED + SKIPPED)) tests"
echo "Passed: ${PASSED}"
echo "Failed: ${FAILED}"
echo "Skipped: ${SKIPPED}"
echo "=========================================="

# Exit with error if any tests failed
if [ $FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
