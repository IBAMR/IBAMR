#ifndef ERROR_CALCULATOR_H
#define ERROR_CALCULATOR_H

#include <vector>
#include <string>
#include <functional>

// SAMRAI headers for patch hierarchy operations
#include <PatchHierarchy.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>

using namespace SAMRAI;

namespace TestUtilities {

/**
 * @brief Utilities for computing errors and convergence rates
 */
class ErrorCalculator {
public:
    /**
     * @brief Compute L2 norm of error on patch hierarchy
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param computed_idx Index of computed solution variable
     * @param exact_idx Index of exact solution variable
     * @param time Current simulation time
     * @return L2 norm of error
     */
    static double computeL2Error(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int computed_idx,
        int exact_idx,
        double time = 0.0);

    /**
     * @brief Compute Linf (max) norm of error on patch hierarchy
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param computed_idx Index of computed solution variable
     * @param exact_idx Index of exact solution variable
     * @return Linf norm of error
     */
    static double computeLinfError(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int computed_idx,
        int exact_idx);

    /**
     * @brief Compute L1 norm of error on patch hierarchy
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param computed_idx Index of computed solution variable
     * @param exact_idx Index of exact solution variable
     * @return L1 norm of error
     */
    static double computeL1Error(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int computed_idx,
        int exact_idx);

    /**
     * @brief Compute convergence rate from error data
     * Uses log-log fit: error = C * h^p, where p is convergence rate
     *
     * @param grid_spacings Vector of grid spacings (h)
     * @param errors Vector of corresponding errors
     * @return Convergence rate p
     */
    static double computeConvergenceRate(
        const std::vector<double>& grid_spacings,
        const std::vector<double>& errors);

    /**
     * @brief Compute convergence rate between two consecutive refinements
     *
     * @param h_coarse Coarse grid spacing
     * @param h_fine Fine grid spacing
     * @param error_coarse Error on coarse grid
     * @param error_fine Error on fine grid
     * @return Convergence rate
     */
    static double computeConvergenceRate(
        double h_coarse, double h_fine,
        double error_coarse, double error_fine);

    /**
     * @brief Compute relative error
     *
     * @param computed Computed value
     * @param exact Exact value
     * @return Relative error: |computed - exact| / |exact|
     */
    static double computeRelativeError(double computed, double exact);

    /**
     * @brief Compute total mass on patch hierarchy
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param var_idx Index of variable
     * @return Total integrated mass
     */
    static double computeTotalMass(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int var_idx);

    /**
     * @brief Check for negative values in field
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param var_idx Index of variable
     * @param tolerance Tolerance for negative values (default -1e-12)
     * @return Number of negative cells found
     */
    static int checkForNegatives(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int var_idx,
        double tolerance = -1.0e-12);

    /**
     * @brief Check for NaN or Inf values in field
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param var_idx Index of variable
     * @return True if NaN/Inf found, false otherwise
     */
    static bool checkForNaNInf(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int var_idx);

    /**
     * @brief Get min and max values in field
     *
     * @param patch_hierarchy Pointer to patch hierarchy
     * @param var_idx Index of variable
     * @param min_val Output: minimum value
     * @param max_val Output: maximum value
     */
    static void getMinMax(
        tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
        int var_idx,
        double& min_val,
        double& max_val);

    /**
     * @brief Print error summary
     *
     * @param test_name Name of test
     * @param L1_error L1 error
     * @param L2_error L2 error
     * @param Linf_error Linf error
     * @param relative_error Relative error (optional)
     */
    static void printErrorSummary(
        const std::string& test_name,
        double L1_error,
        double L2_error,
        double Linf_error,
        double relative_error = -1.0);

    /**
     * @brief Print convergence study results
     *
     * @param test_name Name of test
     * @param grid_spacings Vector of grid spacings
     * @param errors Vector of errors
     * @param expected_rate Expected convergence rate
     * @param tolerance Tolerance for pass/fail
     */
    static void printConvergenceStudy(
        const std::string& test_name,
        const std::vector<double>& grid_spacings,
        const std::vector<double>& errors,
        double expected_rate = 2.0,
        double tolerance = 0.3);

    /**
     * @brief Determine if test passed based on error criteria
     *
     * @param error_value Computed error
     * @param threshold Threshold for passing
     * @return True if test passed, false otherwise
     */
    static bool testPassed(double error_value, double threshold);

    /**
     * @brief Determine if convergence test passed
     *
     * @param computed_rate Computed convergence rate
     * @param expected_rate Expected convergence rate
     * @param tolerance Tolerance (default ±0.3)
     * @return True if test passed, false otherwise
     */
    static bool convergenceTestPassed(
        double computed_rate,
        double expected_rate = 2.0,
        double tolerance = 0.3);

private:
    /**
     * @brief Helper function to perform linear regression
     * Used for convergence rate calculation
     */
    static void linearRegression(
        const std::vector<double>& x,
        const std::vector<double>& y,
        double& slope,
        double& intercept);
};

} // namespace TestUtilities

#endif // ERROR_CALCULATOR_H
