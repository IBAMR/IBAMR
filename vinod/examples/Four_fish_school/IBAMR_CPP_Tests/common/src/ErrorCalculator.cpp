#include "ErrorCalculator.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

// SAMRAI headers
#include <CellData.h>
#include <CellIterator.h>
#include <Patch.h>
#include <PatchLevel.h>

using namespace SAMRAI;

namespace TestUtilities {

double ErrorCalculator::computeL2Error(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int computed_idx,
    int exact_idx,
    double time)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);

    // Create temporary variable for difference
    tbox::Pointer<hier::PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    tbox::Pointer<pdat::CellVariable<NDIM,double>> error_var = new pdat::CellVariable<NDIM,double>("error_temp");

    // Compute difference: error = computed - exact
    hier_cc_data_ops.subtract(computed_idx, computed_idx, exact_idx);

    // Compute L2 norm
    double L2_norm = hier_cc_data_ops.L2Norm(computed_idx);

    return L2_norm;
}

double ErrorCalculator::computeLinfError(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int computed_idx,
    int exact_idx)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);

    // Compute difference
    hier_cc_data_ops.subtract(computed_idx, computed_idx, exact_idx);

    // Compute max norm
    double Linf_norm = hier_cc_data_ops.maxNorm(computed_idx);

    return Linf_norm;
}

double ErrorCalculator::computeL1Error(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int computed_idx,
    int exact_idx)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);

    // Compute difference
    hier_cc_data_ops.subtract(computed_idx, computed_idx, exact_idx);

    // Compute L1 norm
    double L1_norm = hier_cc_data_ops.L1Norm(computed_idx);

    return L1_norm;
}

double ErrorCalculator::computeConvergenceRate(
    const std::vector<double>& grid_spacings,
    const std::vector<double>& errors)
{
    if (grid_spacings.size() != errors.size() || grid_spacings.size() < 2) {
        std::cerr << "ERROR: Need at least 2 data points for convergence rate calculation\n";
        return 0.0;
    }

    // Use log-log fit: log(error) = log(C) + p * log(h)
    std::vector<double> log_h(grid_spacings.size());
    std::vector<double> log_error(errors.size());

    for (size_t i = 0; i < grid_spacings.size(); ++i) {
        log_h[i] = std::log(grid_spacings[i]);
        log_error[i] = std::log(errors[i]);
    }

    double slope, intercept;
    linearRegression(log_h, log_error, slope, intercept);

    return slope;  // This is the convergence rate p
}

double ErrorCalculator::computeConvergenceRate(
    double h_coarse, double h_fine,
    double error_coarse, double error_fine)
{
    if (error_coarse <= 0.0 || error_fine <= 0.0) {
        std::cerr << "ERROR: Errors must be positive for convergence rate calculation\n";
        return 0.0;
    }

    double log_ratio_h = std::log(h_coarse / h_fine);
    double log_ratio_error = std::log(error_coarse / error_fine);

    return log_ratio_error / log_ratio_h;
}

double ErrorCalculator::computeRelativeError(double computed, double exact)
{
    if (std::abs(exact) < 1.0e-14) {
        return std::abs(computed - exact);  // Absolute error if exact is near zero
    }
    return std::abs(computed - exact) / std::abs(exact);
}

double ErrorCalculator::computeTotalMass(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int var_idx)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    const double weight = 1.0;  // Weight for integration

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);

    // Integrate over domain (sum of C * dV)
    double total_mass = hier_cc_data_ops.integral(var_idx, weight);

    return total_mass;
}

int ErrorCalculator::checkForNegatives(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int var_idx,
    double tolerance)
{
    int negative_count = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln) {
        tbox::Pointer<hier::PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);

        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            tbox::Pointer<pdat::CellData<NDIM,double>> data = patch->getPatchData(var_idx);

            const hier::Box<NDIM>& patch_box = patch->getBox();
            for (pdat::CellIterator<NDIM> ci(patch_box); ci; ci++) {
                const pdat::CellIndex<NDIM>& idx = ci();
                if ((*data)(idx) < tolerance) {
                    negative_count++;
                }
            }
        }
    }

    return negative_count;
}

bool ErrorCalculator::checkForNaNInf(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int var_idx)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln) {
        tbox::Pointer<hier::PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);

        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            tbox::Pointer<pdat::CellData<NDIM,double>> data = patch->getPatchData(var_idx);

            const hier::Box<NDIM>& patch_box = patch->getBox();
            for (pdat::CellIterator<NDIM> ci(patch_box); ci; ci++) {
                const pdat::CellIndex<NDIM>& idx = ci();
                double value = (*data)(idx);
                if (std::isnan(value) || std::isinf(value)) {
                    return true;
                }
            }
        }
    }

    return false;
}

void ErrorCalculator::getMinMax(
    tbox::Pointer<hier::PatchHierarchy<NDIM>> patch_hierarchy,
    int var_idx,
    double& min_val,
    double& max_val)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);

    min_val = hier_cc_data_ops.min(var_idx);
    max_val = hier_cc_data_ops.max(var_idx);
}

void ErrorCalculator::printErrorSummary(
    const std::string& test_name,
    double L1_error,
    double L2_error,
    double Linf_error,
    double relative_error)
{
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "  ERROR SUMMARY: " << test_name << "\n";
    std::cout << "========================================\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  L1 error:    " << L1_error << "\n";
    std::cout << "  L2 error:    " << L2_error << "\n";
    std::cout << "  Linf error:  " << Linf_error << "\n";
    if (relative_error >= 0.0) {
        std::cout << "  Relative:    " << relative_error << "\n";
    }
    std::cout << "========================================\n\n";
}

void ErrorCalculator::printConvergenceStudy(
    const std::string& test_name,
    const std::vector<double>& grid_spacings,
    const std::vector<double>& errors,
    double expected_rate,
    double tolerance)
{
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "  CONVERGENCE STUDY: " << test_name << "\n";
    std::cout << "========================================\n";

    std::cout << std::setw(15) << "Grid Spacing"
              << std::setw(15) << "Error"
              << std::setw(15) << "Rate\n";
    std::cout << "----------------------------------------\n";

    for (size_t i = 0; i < grid_spacings.size(); ++i) {
        std::cout << std::scientific << std::setprecision(4);
        std::cout << std::setw(15) << grid_spacings[i]
                  << std::setw(15) << errors[i];

        if (i > 0) {
            double rate = computeConvergenceRate(
                grid_spacings[i-1], grid_spacings[i],
                errors[i-1], errors[i]);
            std::cout << std::setw(15) << rate;
        }
        std::cout << "\n";
    }

    double overall_rate = computeConvergenceRate(grid_spacings, errors);
    std::cout << "----------------------------------------\n";
    std::cout << "Overall convergence rate: " << overall_rate << "\n";
    std::cout << "Expected rate: " << expected_rate << " +/- " << tolerance << "\n";

    bool passed = convergenceTestPassed(overall_rate, expected_rate, tolerance);
    std::cout << "Result: " << (passed ? "PASSED" : "FAILED") << "\n";
    std::cout << "========================================\n\n";
}

bool ErrorCalculator::testPassed(double error_value, double threshold)
{
    return error_value < threshold;
}

bool ErrorCalculator::convergenceTestPassed(
    double computed_rate,
    double expected_rate,
    double tolerance)
{
    double difference = std::abs(computed_rate - expected_rate);
    return difference <= tolerance;
}

void ErrorCalculator::linearRegression(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double& slope,
    double& intercept)
{
    size_t n = x.size();
    if (n != y.size() || n < 2) {
        slope = 0.0;
        intercept = 0.0;
        return;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;

    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_xx += x[i] * x[i];
    }

    double n_double = static_cast<double>(n);
    double denominator = n_double * sum_xx - sum_x * sum_x;

    if (std::abs(denominator) < 1.0e-14) {
        slope = 0.0;
        intercept = sum_y / n_double;
        return;
    }

    slope = (n_double * sum_xy - sum_x * sum_y) / denominator;
    intercept = (sum_y - slope * sum_x) / n_double;
}

} // namespace TestUtilities
