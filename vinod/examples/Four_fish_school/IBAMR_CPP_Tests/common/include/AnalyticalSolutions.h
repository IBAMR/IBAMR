#ifndef ANALYTICAL_SOLUTIONS_H
#define ANALYTICAL_SOLUTIONS_H

#include <cmath>
#include <vector>

namespace TestUtilities {

/**
 * @brief Collection of analytical solutions for validation tests
 */
class AnalyticalSolutions {
public:
    /**
     * @brief Gaussian diffusion solution in 1D
     * C(x,t) = (1/sqrt(4*pi*kappa*t)) * exp(-x^2 / (4*kappa*t))
     *
     * @param x Position
     * @param t Time
     * @param kappa Diffusion coefficient
     * @param x0 Initial center position
     * @param C0 Initial amplitude
     * @return Concentration at (x,t)
     */
    static double gaussianDiffusion1D(double x, double t, double kappa,
                                     double x0 = 0.0, double C0 = 1.0);

    /**
     * @brief Gaussian diffusion solution in 2D
     * C(x,y,t) = (C0/(4*pi*kappa*t)) * exp(-(r^2)/(4*kappa*t))
     * where r^2 = (x-x0)^2 + (y-y0)^2
     *
     * @param x X position
     * @param y Y position
     * @param t Time
     * @param kappa Diffusion coefficient
     * @param x0 Initial center X
     * @param y0 Initial center Y
     * @param C0 Initial amplitude (total mass)
     * @return Concentration at (x,y,t)
     */
    static double gaussianDiffusion2D(double x, double y, double t, double kappa,
                                     double x0 = 0.0, double y0 = 0.0, double C0 = 1.0);

    /**
     * @brief Gaussian diffusion solution in 3D
     *
     * @param x X position
     * @param y Y position
     * @param z Z position
     * @param t Time
     * @param kappa Diffusion coefficient
     * @param x0 Initial center X
     * @param y0 Initial center Y
     * @param z0 Initial center Z
     * @param C0 Initial amplitude
     * @return Concentration at (x,y,z,t)
     */
    static double gaussianDiffusion3D(double x, double y, double z, double t, double kappa,
                                     double x0 = 0.0, double y0 = 0.0, double z0 = 0.0,
                                     double C0 = 1.0);

    /**
     * @brief Pure advection of Gaussian profile
     * C(x,t) = C0 * exp(-((x - x0 - u*t)^2) / (2*sigma^2))
     *
     * @param x Position
     * @param t Time
     * @param u Advection velocity
     * @param x0 Initial center
     * @param sigma Standard deviation
     * @param C0 Amplitude
     * @return Concentration at (x,t)
     */
    static double advectedGaussian1D(double x, double t, double u,
                                    double x0 = 0.0, double sigma = 0.1, double C0 = 1.0);

    /**
     * @brief Pure advection of Gaussian profile in 2D
     *
     * @param x X position
     * @param y Y position
     * @param t Time
     * @param u X velocity
     * @param v Y velocity
     * @param x0 Initial center X
     * @param y0 Initial center Y
     * @param sigma Standard deviation
     * @param C0 Amplitude
     * @return Concentration at (x,y,t)
     */
    static double advectedGaussian2D(double x, double y, double t, double u, double v,
                                    double x0 = 0.0, double y0 = 0.0,
                                    double sigma = 0.1, double C0 = 1.0);

    /**
     * @brief Manufactured solution for MMS test
     * C(x,y,t) = exp(-t) * sin(pi*x) * sin(pi*y)
     *
     * @param x X position
     * @param y Y position
     * @param t Time
     * @return Concentration
     */
    static double manufacturedSolution2D(double x, double y, double t);

    /**
     * @brief Source term for manufactured solution
     * S = dC/dt - kappa*Laplacian(C) + u*dC/dx + v*dC/dy
     *
     * @param x X position
     * @param y Y position
     * @param t Time
     * @param kappa Diffusion coefficient
     * @param u X velocity
     * @param v Y velocity
     * @return Source term value
     */
    static double manufacturedSource2D(double x, double y, double t,
                                      double kappa, double u, double v);

    /**
     * @brief Top-hat (step) function
     * C(x) = C0 if |x - x0| < width/2, else 0
     *
     * @param x Position
     * @param x0 Center position
     * @param width Width of top-hat
     * @param C0 Amplitude
     * @return Concentration
     */
    static double topHat1D(double x, double x0, double width, double C0 = 1.0);

    /**
     * @brief Top-hat function in 2D (cylinder)
     *
     * @param x X position
     * @param y Y position
     * @param x0 Center X
     * @param y0 Center Y
     * @param radius Cylinder radius
     * @param C0 Amplitude
     * @return Concentration
     */
    static double topHat2D(double x, double y, double x0, double y0,
                          double radius, double C0 = 1.0);

    /**
     * @brief Steady-state diffusion around a sphere source
     * For low Reynolds number, steady-state solution:
     * C(r) = (Q / (4*pi*kappa)) * (1/r) for r > R
     * where Q is source strength, R is sphere radius
     *
     * @param x X position
     * @param y Y position
     * @param z Z position
     * @param x0 Sphere center X
     * @param y0 Sphere center Y
     * @param z0 Sphere center Z
     * @param R Sphere radius
     * @param Q Source strength
     * @param kappa Diffusion coefficient
     * @return Concentration
     */
    static double sphereSource3D(double x, double y, double z,
                                double x0, double y0, double z0,
                                double R, double Q, double kappa);

    /**
     * @brief Steady-state diffusion around a cylinder source (2D)
     *
     * @param x X position
     * @param y Y position
     * @param x0 Cylinder center X
     * @param y0 Cylinder center Y
     * @param R Cylinder radius
     * @param Q Source strength per unit length
     * @param kappa Diffusion coefficient
     * @return Concentration
     */
    static double cylinderSource2D(double x, double y,
                                  double x0, double y0,
                                  double R, double Q, double kappa);

private:
    static constexpr double PI = 3.141592653589793238463;
    static constexpr double EPSILON = 1.0e-14;  // Small number to avoid division by zero
};

} // namespace TestUtilities

#endif // ANALYTICAL_SOLUTIONS_H
