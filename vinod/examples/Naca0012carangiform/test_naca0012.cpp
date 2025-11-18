// Simple test to verify NACA0012 thickness distribution
#include <iostream>
#include <cmath>
#include <iomanip>

static const double CHORD_LENGTH = 1.0;
static const double MAX_THICKNESS = 0.12;

// NACA 4-digit airfoil thickness distribution coefficients
static const double NACA_A0 =  0.2969;
static const double NACA_A1 = -0.1260;
static const double NACA_A2 = -0.3516;
static const double NACA_A3 =  0.2843;
static const double NACA_A4 = -0.1015;

inline double naca0012_thickness(const double x_c)
{
    if (x_c < 0.0 || x_c > 1.0) return 0.0;

    const double sqrt_x_c = std::sqrt(x_c);
    const double x_c2 = x_c * x_c;
    const double x_c3 = x_c2 * x_c;
    const double x_c4 = x_c3 * x_c;

    // NACA 4-digit thickness distribution
    const double y_t = (MAX_THICKNESS / 0.2) * (
        NACA_A0 * sqrt_x_c +
        NACA_A1 * x_c +
        NACA_A2 * x_c2 +
        NACA_A3 * x_c3 +
        NACA_A4 * x_c4
    );

    return y_t * CHORD_LENGTH;
}

int main()
{
    std::cout << "NACA0012 Thickness Distribution Test\n";
    std::cout << "=====================================\n";
    std::cout << "Chord Length: " << CHORD_LENGTH << "\n";
    std::cout << "Max Thickness: " << MAX_THICKNESS << " (12% of chord)\n\n";
    std::cout << std::setw(10) << "x/c" << std::setw(15) << "Half-thickness" << "\n";
    std::cout << std::string(25, '-') << "\n";

    double max_t = 0.0;
    double x_at_max_t = 0.0;

    for (double x_c = 0.0; x_c <= 1.0; x_c += 0.1)
    {
        double t = naca0012_thickness(x_c);
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << x_c
                  << std::setw(15) << std::setprecision(6) << t << "\n";

        if (t > max_t)
        {
            max_t = t;
            x_at_max_t = x_c;
        }
    }

    std::cout << "\nMaximum half-thickness: " << max_t
              << " at x/c = " << x_at_max_t << "\n";
    std::cout << "Full thickness: " << 2.0 * max_t << "\n";
    std::cout << "Expected max thickness: " << MAX_THICKNESS * CHORD_LENGTH << "\n";

    return 0;
}
