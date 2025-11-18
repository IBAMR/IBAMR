#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <string>
#include <fstream>
#include <iomanip>
#include <vector>

// SAMRAI/IBAMR headers
#include <PatchHierarchy.h>
#include <tbox/Pointer.h>

using namespace SAMRAI;

namespace TestUtilities {

/**
 * @brief General utilities for test infrastructure
 */
class TestUtils {
public:
    /**
     * @brief Print test header
     *
     * @param test_name Name of the test
     * @param test_number Test number (1-14)
     */
    static void printTestHeader(const std::string& test_name, int test_number);

    /**
     * @brief Print test footer with verdict
     *
     * @param test_name Name of the test
     * @param passed True if test passed, false otherwise
     * @param message Optional message to display
     */
    static void printTestFooter(const std::string& test_name, bool passed,
                               const std::string& message = "");

    /**
     * @brief Print section separator
     *
     * @param title Section title
     */
    static void printSectionSeparator(const std::string& title);

    /**
     * @brief Print progress message
     *
     * @param message Progress message
     */
    static void printProgress(const std::string& message);

    /**
     * @brief Print warning message
     *
     * @param message Warning message
     */
    static void printWarning(const std::string& message);

    /**
     * @brief Print error message
     *
     * @param message Error message
     */
    static void printError(const std::string& message);

    /**
     * @brief Write results to file
     *
     * @param filename Output filename
     * @param content Content to write
     * @param append If true, append to file; otherwise overwrite
     */
    static void writeToFile(const std::string& filename,
                           const std::string& content,
                           bool append = false);

    /**
     * @brief Create results summary
     *
     * @param test_name Name of test
     * @param passed Test result
     * @param details Additional details
     * @return Formatted summary string
     */
    static std::string createResultsSummary(
        const std::string& test_name,
        bool passed,
        const std::string& details);

    /**
     * @brief Check mass conservation
     * Computes relative change in total mass
     *
     * @param initial_mass Initial total mass
     * @param current_mass Current total mass
     * @param tolerance Tolerance for drift
     * @return True if mass conserved within tolerance
     */
    static bool checkMassConservation(
        double initial_mass,
        double current_mass,
        double tolerance = 1.0e-10);

    /**
     * @brief Compute relative mass drift
     *
     * @param initial_mass Initial total mass
     * @param current_mass Current total mass
     * @return Relative mass drift: |M_current - M_initial| / M_initial
     */
    static double computeMassDrift(double initial_mass, double current_mass);

    /**
     * @brief Validate simulation parameters
     * Checks for reasonable CFL, Schmidt number, etc.
     *
     * @param dt Time step
     * @param dx Grid spacing
     * @param velocity Characteristic velocity
     * @param diffusivity Diffusion coefficient
     * @return True if parameters are reasonable
     */
    static bool validateParameters(
        double dt,
        double dx,
        double velocity,
        double diffusivity);

    /**
     * @brief Compute CFL number for advection
     * CFL_adv = u * dt / dx
     *
     * @param velocity Characteristic velocity
     * @param dt Time step
     * @param dx Grid spacing
     * @return CFL number
     */
    static double computeCFL_advection(double velocity, double dt, double dx);

    /**
     * @brief Compute CFL number for diffusion
     * CFL_diff = kappa * dt / dx^2
     *
     * @param diffusivity Diffusion coefficient
     * @param dt Time step
     * @param dx Grid spacing
     * @return Diffusion CFL number
     */
    static double computeCFL_diffusion(double diffusivity, double dt, double dx);

    /**
     * @brief Compute Schmidt number
     * Sc = nu / kappa (kinematic viscosity / diffusivity)
     *
     * @param viscosity Kinematic viscosity
     * @param diffusivity Diffusion coefficient
     * @return Schmidt number
     */
    static double computeSchmidtNumber(double viscosity, double diffusivity);

    /**
     * @brief Compute Peclet number
     * Pe = U * L / kappa
     *
     * @param velocity Characteristic velocity
     * @param length Characteristic length
     * @param diffusivity Diffusion coefficient
     * @return Peclet number
     */
    static double computePecletNumber(
        double velocity,
        double length,
        double diffusivity);

    /**
     * @brief Format time string
     * Converts seconds to human-readable format
     *
     * @param seconds Time in seconds
     * @return Formatted string (e.g., "2h 34m 12s")
     */
    static std::string formatTime(double seconds);

    /**
     * @brief Get timestamp string
     *
     * @return Current timestamp in format "YYYY-MM-DD HH:MM:SS"
     */
    static std::string getTimestamp();

    /**
     * @brief Create directory if it doesn't exist
     *
     * @param dir_path Directory path
     * @return True if successful or already exists
     */
    static bool createDirectory(const std::string& dir_path);

    /**
     * @brief Check if file exists
     *
     * @param filename File path
     * @return True if file exists
     */
    static bool fileExists(const std::string& filename);

    /**
     * @brief Format double with fixed precision
     *
     * @param value Value to format
     * @param precision Number of decimal places
     * @return Formatted string
     */
    static std::string formatDouble(double value, int precision = 6);

    /**
     * @brief Format double in scientific notation
     *
     * @param value Value to format
     * @param precision Number of significant figures
     * @return Formatted string
     */
    static std::string formatScientific(double value, int precision = 3);

    /**
     * @brief Create a table row for results
     *
     * @param columns Vector of column values
     * @param widths Vector of column widths
     * @return Formatted table row string
     */
    static std::string createTableRow(
        const std::vector<std::string>& columns,
        const std::vector<int>& widths);

    /**
     * @brief Create table header
     *
     * @param headers Vector of header names
     * @param widths Vector of column widths
     * @return Formatted header string with separator
     */
    static std::string createTableHeader(
        const std::vector<std::string>& headers,
        const std::vector<int>& widths);

private:
    static constexpr int LINE_WIDTH = 80;
    static constexpr char SEPARATOR_CHAR = '=';
    static constexpr char SUBSEPARATOR_CHAR = '-';
};

/**
 * @brief Timer class for performance measurement
 */
class Timer {
public:
    Timer();

    void start();
    void stop();
    void reset();

    double getElapsedTime() const;  // Returns time in seconds
    std::string getFormattedTime() const;

private:
    double start_time_;
    double elapsed_time_;
    bool running_;

    static double getCurrentTime();  // Platform-independent time
};

/**
 * @brief Result logger for systematic result recording
 */
class ResultLogger {
public:
    ResultLogger(const std::string& filename);
    ~ResultLogger();

    void log(const std::string& message);
    void logError(const std::string& error_type, double error_value);
    void logParameter(const std::string& param_name, double param_value);
    void logParameter(const std::string& param_name, const std::string& param_value);
    void logTestResult(const std::string& test_name, bool passed);
    void logConvergenceRate(double rate, double expected, bool passed);
    void close();

private:
    std::ofstream file_;
    std::string filename_;
    bool is_open_;
};

} // namespace TestUtilities

#endif // TEST_UTILITIES_H
