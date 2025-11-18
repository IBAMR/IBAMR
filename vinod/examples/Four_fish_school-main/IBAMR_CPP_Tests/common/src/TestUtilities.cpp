#include "TestUtilities.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <sys/stat.h>

namespace TestUtilities {

// ============================================================================
// TestUtils Implementation
// ============================================================================

void TestUtils::printTestHeader(const std::string& test_name, int test_number)
{
    std::cout << "\n";
    std::cout << std::string(LINE_WIDTH, SEPARATOR_CHAR) << "\n";
    std::cout << "  TEST " << test_number << ": " << test_name << "\n";
    std::cout << std::string(LINE_WIDTH, SEPARATOR_CHAR) << "\n";
    std::cout << "  Start time: " << getTimestamp() << "\n";
    std::cout << std::string(LINE_WIDTH, SEPARATOR_CHAR) << "\n\n";
}

void TestUtils::printTestFooter(const std::string& test_name, bool passed,
                                const std::string& message)
{
    std::cout << "\n";
    std::cout << std::string(LINE_WIDTH, SEPARATOR_CHAR) << "\n";
    std::cout << "  TEST RESULT: " << (passed ? "PASSED" : "FAILED") << "\n";
    if (!message.empty()) {
        std::cout << "  " << message << "\n";
    }
    std::cout << "  End time: " << getTimestamp() << "\n";
    std::cout << std::string(LINE_WIDTH, SEPARATOR_CHAR) << "\n\n";
}

void TestUtils::printSectionSeparator(const std::string& title)
{
    std::cout << "\n";
    std::cout << std::string(LINE_WIDTH, SUBSEPARATOR_CHAR) << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(LINE_WIDTH, SUBSEPARATOR_CHAR) << "\n\n";
}

void TestUtils::printProgress(const std::string& message)
{
    std::cout << "[INFO] " << message << std::endl;
}

void TestUtils::printWarning(const std::string& message)
{
    std::cout << "[WARNING] " << message << std::endl;
}

void TestUtils::printError(const std::string& message)
{
    std::cerr << "[ERROR] " << message << std::endl;
}

void TestUtils::writeToFile(const std::string& filename,
                            const std::string& content,
                            bool append)
{
    std::ios_base::openmode mode = append ? std::ios::app : std::ios::out;
    std::ofstream file(filename, mode);

    if (!file.is_open()) {
        printError("Could not open file: " + filename);
        return;
    }

    file << content;
    file.close();
}

std::string TestUtils::createResultsSummary(
    const std::string& test_name,
    bool passed,
    const std::string& details)
{
    std::ostringstream oss;

    oss << "Test: " << test_name << "\n";
    oss << "Result: " << (passed ? "PASSED" : "FAILED") << "\n";
    oss << "Timestamp: " << getTimestamp() << "\n";
    if (!details.empty()) {
        oss << "Details:\n" << details << "\n";
    }
    oss << std::string(60, '-') << "\n";

    return oss.str();
}

bool TestUtils::checkMassConservation(
    double initial_mass,
    double current_mass,
    double tolerance)
{
    double drift = computeMassDrift(initial_mass, current_mass);
    return drift < tolerance;
}

double TestUtils::computeMassDrift(double initial_mass, double current_mass)
{
    if (std::abs(initial_mass) < 1.0e-14) {
        return std::abs(current_mass - initial_mass);
    }
    return std::abs(current_mass - initial_mass) / std::abs(initial_mass);
}

bool TestUtils::validateParameters(
    double dt,
    double dx,
    double velocity,
    double diffusivity)
{
    bool valid = true;

    double CFL_adv = computeCFL_advection(velocity, dt, dx);
    double CFL_diff = computeCFL_diffusion(diffusivity, dt, dx);

    if (CFL_adv > 1.0) {
        printWarning("Advection CFL > 1.0: " + formatDouble(CFL_adv, 3));
        valid = false;
    }

    if (CFL_diff > 0.5) {
        printWarning("Diffusion CFL > 0.5: " + formatDouble(CFL_diff, 3));
        valid = false;
    }

    return valid;
}

double TestUtils::computeCFL_advection(double velocity, double dt, double dx)
{
    return std::abs(velocity) * dt / dx;
}

double TestUtils::computeCFL_diffusion(double diffusivity, double dt, double dx)
{
    return diffusivity * dt / (dx * dx);
}

double TestUtils::computeSchmidtNumber(double viscosity, double diffusivity)
{
    if (std::abs(diffusivity) < 1.0e-14) {
        printWarning("Diffusivity near zero in Schmidt number calculation");
        return 1.0e10;
    }
    return viscosity / diffusivity;
}

double TestUtils::computePecletNumber(
    double velocity,
    double length,
    double diffusivity)
{
    if (std::abs(diffusivity) < 1.0e-14) {
        printWarning("Diffusivity near zero in Peclet number calculation");
        return 1.0e10;
    }
    return velocity * length / diffusivity;
}

std::string TestUtils::formatTime(double seconds)
{
    int hours = static_cast<int>(seconds) / 3600;
    int minutes = (static_cast<int>(seconds) % 3600) / 60;
    int secs = static_cast<int>(seconds) % 60;

    std::ostringstream oss;
    if (hours > 0) {
        oss << hours << "h " << minutes << "m " << secs << "s";
    } else if (minutes > 0) {
        oss << minutes << "m " << secs << "s";
    } else {
        oss << secs << "s";
    }

    return oss.str();
}

std::string TestUtils::getTimestamp()
{
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    auto now_tm = *std::localtime(&now_c);

    std::ostringstream oss;
    oss << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

bool TestUtils::createDirectory(const std::string& dir_path)
{
    struct stat st;
    if (stat(dir_path.c_str(), &st) == 0) {
        return true;  // Directory already exists
    }

    // Try to create directory
    int result = mkdir(dir_path.c_str(), 0755);
    return (result == 0);
}

bool TestUtils::fileExists(const std::string& filename)
{
    std::ifstream file(filename);
    return file.good();
}

std::string TestUtils::formatDouble(double value, int precision)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

std::string TestUtils::formatScientific(double value, int precision)
{
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << value;
    return oss.str();
}

std::string TestUtils::createTableRow(
    const std::vector<std::string>& columns,
    const std::vector<int>& widths)
{
    if (columns.size() != widths.size()) {
        return "";
    }

    std::ostringstream oss;
    for (size_t i = 0; i < columns.size(); ++i) {
        oss << std::setw(widths[i]) << columns[i];
    }
    oss << "\n";

    return oss.str();
}

std::string TestUtils::createTableHeader(
    const std::vector<std::string>& headers,
    const std::vector<int>& widths)
{
    std::ostringstream oss;

    // Header row
    oss << createTableRow(headers, widths);

    // Separator
    int total_width = 0;
    for (int w : widths) total_width += w;
    oss << std::string(total_width, SUBSEPARATOR_CHAR) << "\n";

    return oss.str();
}

// ============================================================================
// Timer Implementation
// ============================================================================

Timer::Timer() : start_time_(0.0), elapsed_time_(0.0), running_(false)
{
}

void Timer::start()
{
    start_time_ = getCurrentTime();
    running_ = true;
}

void Timer::stop()
{
    if (running_) {
        elapsed_time_ += getCurrentTime() - start_time_;
        running_ = false;
    }
}

void Timer::reset()
{
    start_time_ = 0.0;
    elapsed_time_ = 0.0;
    running_ = false;
}

double Timer::getElapsedTime() const
{
    double total = elapsed_time_;
    if (running_) {
        total += getCurrentTime() - start_time_;
    }
    return total;
}

std::string Timer::getFormattedTime() const
{
    return TestUtils::formatTime(getElapsedTime());
}

double Timer::getCurrentTime()
{
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = now.time_since_epoch();
    auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(duration);
    return seconds.count();
}

// ============================================================================
// ResultLogger Implementation
// ============================================================================

ResultLogger::ResultLogger(const std::string& filename)
    : filename_(filename), is_open_(false)
{
    file_.open(filename_);
    if (file_.is_open()) {
        is_open_ = true;
        file_ << "Test Results Log\n";
        file_ << "Generated: " << TestUtils::getTimestamp() << "\n";
        file_ << std::string(60, '=') << "\n\n";
    } else {
        TestUtils::printError("Could not open log file: " + filename_);
    }
}

ResultLogger::~ResultLogger()
{
    close();
}

void ResultLogger::log(const std::string& message)
{
    if (is_open_) {
        file_ << message << "\n";
        file_.flush();
    }
}

void ResultLogger::logError(const std::string& error_type, double error_value)
{
    if (is_open_) {
        file_ << std::scientific << std::setprecision(6);
        file_ << error_type << ": " << error_value << "\n";
        file_.flush();
    }
}

void ResultLogger::logParameter(const std::string& param_name, double param_value)
{
    if (is_open_) {
        file_ << param_name << ": " << param_value << "\n";
        file_.flush();
    }
}

void ResultLogger::logParameter(const std::string& param_name, const std::string& param_value)
{
    if (is_open_) {
        file_ << param_name << ": " << param_value << "\n";
        file_.flush();
    }
}

void ResultLogger::logTestResult(const std::string& test_name, bool passed)
{
    if (is_open_) {
        file_ << "\nTest: " << test_name << "\n";
        file_ << "Result: " << (passed ? "PASSED" : "FAILED") << "\n";
        file_ << "Time: " << TestUtils::getTimestamp() << "\n";
        file_ << std::string(60, '-') << "\n\n";
        file_.flush();
    }
}

void ResultLogger::logConvergenceRate(double rate, double expected, bool passed)
{
    if (is_open_) {
        file_ << std::fixed << std::setprecision(3);
        file_ << "Convergence rate: " << rate << "\n";
        file_ << "Expected rate: " << expected << "\n";
        file_ << "Result: " << (passed ? "PASSED" : "FAILED") << "\n";
        file_.flush();
    }
}

void ResultLogger::close()
{
    if (is_open_) {
        file_ << "\nLog closed: " << TestUtils::getTimestamp() << "\n";
        file_.close();
        is_open_ = false;
    }
}

} // namespace TestUtilities
