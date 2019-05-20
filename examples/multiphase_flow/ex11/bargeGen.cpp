#include <cmath>
#include <fstream>
#include <iomanip>

int
main()
{
    const double Lx = 5.0;
    const double Ly = 2.5;
    const int Nx = 500 * 2 * 2;
    const int Ny = 250 * 2 * 2;

    // Dimensional parameters
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    const double Length = 0.3;
    const double Width = Length / 3.0;
    const double Xcom = 0.0;
    const double Ycom = 0.0;
    const double theta = 15;             // in degrees
    const double t = theta * M_PI / 180; // in radians

    const int NumPtsX = std::ceil(Length / dx) + 1;
    const int NumPtsY = std::ceil(Width / dy) + 1;

    std::ofstream barge_coord_stream("barge2d.vertex");
    barge_coord_stream << NumPtsX * NumPtsY << "\n";
    // Generate Rectangular Prism
    for (int j = 1; j <= NumPtsY; ++j)
    {
        const double y = Ycom - Width / 2 + (j - 1) * dy;
        for (int i = 1; i <= NumPtsX; ++i)
        {
            const double x = Xcom - Length / 2 + (i - 1) * dy;
            const double r_x = x * std::cos(t) - y * std::sin(t);
            const double r_y = x * std::sin(t) + y * std::cos(t);
            barge_coord_stream << std::scientific << std::setprecision(7) << std::setfill('0') << r_x << "\t" << r_y
                               << "\n";
        }
    }
}
