#include "AnalyticalSolutions.h"
#include <cmath>
#include <algorithm>

namespace TestUtilities {

double AnalyticalSolutions::gaussianDiffusion1D(
    double x, double t, double kappa, double x0, double C0)
{
    if (t < EPSILON) {
        // At t=0, return delta function approximation
        return (std::abs(x - x0) < EPSILON) ? C0 : 0.0;
    }

    double dx = x - x0;
    double denom = std::sqrt(4.0 * PI * kappa * t);
    double exponent = -(dx * dx) / (4.0 * kappa * t);

    return (C0 / denom) * std::exp(exponent);
}

double AnalyticalSolutions::gaussianDiffusion2D(
    double x, double y, double t, double kappa,
    double x0, double y0, double C0)
{
    if (t < EPSILON) {
        // At t=0, return delta function approximation
        double dx = x - x0;
        double dy = y - y0;
        double r = std::sqrt(dx*dx + dy*dy);
        return (r < EPSILON) ? C0 : 0.0;
    }

    double dx = x - x0;
    double dy = y - y0;
    double r_squared = dx*dx + dy*dy;

    double denom = 4.0 * PI * kappa * t;
    double exponent = -r_squared / (4.0 * kappa * t);

    return (C0 / denom) * std::exp(exponent);
}

double AnalyticalSolutions::gaussianDiffusion3D(
    double x, double y, double z, double t, double kappa,
    double x0, double y0, double z0, double C0)
{
    if (t < EPSILON) {
        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        return (r < EPSILON) ? C0 : 0.0;
    }

    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;
    double r_squared = dx*dx + dy*dy + dz*dz;

    double denom = std::pow(4.0 * PI * kappa * t, 1.5);
    double exponent = -r_squared / (4.0 * kappa * t);

    return (C0 / denom) * std::exp(exponent);
}

double AnalyticalSolutions::advectedGaussian1D(
    double x, double t, double u, double x0, double sigma, double C0)
{
    double x_center = x0 + u * t;
    double dx = x - x_center;
    double exponent = -(dx * dx) / (2.0 * sigma * sigma);

    return C0 * std::exp(exponent);
}

double AnalyticalSolutions::advectedGaussian2D(
    double x, double y, double t, double u, double v,
    double x0, double y0, double sigma, double C0)
{
    double x_center = x0 + u * t;
    double y_center = y0 + v * t;

    double dx = x - x_center;
    double dy = y - y_center;
    double r_squared = dx*dx + dy*dy;

    double exponent = -r_squared / (2.0 * sigma * sigma);

    return C0 * std::exp(exponent);
}

double AnalyticalSolutions::manufacturedSolution2D(double x, double y, double t)
{
    return std::exp(-t) * std::sin(PI * x) * std::sin(PI * y);
}

double AnalyticalSolutions::manufacturedSource2D(
    double x, double y, double t, double kappa, double u, double v)
{
    // C = exp(-t) * sin(pi*x) * sin(pi*y)
    // dC/dt = -exp(-t) * sin(pi*x) * sin(pi*y) = -C
    // dC/dx = exp(-t) * pi * cos(pi*x) * sin(pi*y)
    // dC/dy = exp(-t) * pi * sin(pi*x) * cos(pi*y)
    // d2C/dx2 = -exp(-t) * pi^2 * sin(pi*x) * sin(pi*y) = -pi^2 * C
    // d2C/dy2 = -exp(-t) * pi^2 * sin(pi*x) * sin(pi*y) = -pi^2 * C
    // Laplacian(C) = -2*pi^2 * C
    //
    // Source: S = dC/dt - kappa*Laplacian(C) + u*dC/dx + v*dC/dy

    double C = manufacturedSolution2D(x, y, t);

    double dC_dt = -C;
    double dC_dx = std::exp(-t) * PI * std::cos(PI * x) * std::sin(PI * y);
    double dC_dy = std::exp(-t) * PI * std::sin(PI * x) * std::cos(PI * y);
    double laplacian_C = -2.0 * PI * PI * C;

    double source = dC_dt - kappa * laplacian_C + u * dC_dx + v * dC_dy;

    return source;
}

double AnalyticalSolutions::topHat1D(double x, double x0, double width, double C0)
{
    double distance = std::abs(x - x0);
    return (distance < width / 2.0) ? C0 : 0.0;
}

double AnalyticalSolutions::topHat2D(
    double x, double y, double x0, double y0, double radius, double C0)
{
    double dx = x - x0;
    double dy = y - y0;
    double r = std::sqrt(dx*dx + dy*dy);

    return (r < radius) ? C0 : 0.0;
}

double AnalyticalSolutions::sphereSource3D(
    double x, double y, double z,
    double x0, double y0, double z0,
    double R, double Q, double kappa)
{
    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Avoid singularity at sphere surface
    if (r < R + EPSILON) {
        return Q / (4.0 * PI * kappa * R);
    }

    return Q / (4.0 * PI * kappa * r);
}

double AnalyticalSolutions::cylinderSource2D(
    double x, double y, double x0, double y0,
    double R, double Q, double kappa)
{
    double dx = x - x0;
    double dy = y - y0;
    double r = std::sqrt(dx*dx + dy*dy);

    // Avoid singularity at cylinder surface
    if (r < R + EPSILON) {
        // At surface: C = -(Q/(2*pi*kappa)) * ln(R)
        return -(Q / (2.0 * PI * kappa)) * std::log(R);
    }

    // For r > R: C = -(Q/(2*pi*kappa)) * ln(r)
    // Note: This is relative to C->0 as r->infinity
    return -(Q / (2.0 * PI * kappa)) * std::log(r);
}

} // namespace TestUtilities
