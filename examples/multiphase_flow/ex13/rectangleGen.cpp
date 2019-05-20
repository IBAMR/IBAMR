#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <vector>

int main()
{
    const double Lx = 3.22;
    const double Ly = 1.0;
    const double Lz = 1.0;
    const double Nx = 161;
    const double Ny = 50;
    const double Nz = 50;

    const double dx = Lx/Nx;
    const double dy = Ly/Ny;
    const double dz = Lz/Nz;

    const double Length = 0.16;
    const double Width = 0.4;
    const double Height = 0.16;
    const double Xcom = 0.74;
    const double Ycom = 0.3 + Width/2;
    const double Zcom = Height/2.0;

    const std::size_t NumPtsX = std::ceil(Length/dx) + 1;
    const std::size_t NumPtsY = std::ceil(Width/dy) + 1;
    const std::size_t NumPtsZ = std::ceil(Height/dz) + 1;

    std::vector<double> LagX;
    std::vector<double> LagY;
    std::vector<double> LagZ;
    // Generate rectangular prism
    for (std::size_t k = 0; k < NumPtsZ; ++k)
    {
        const double z = Zcom - Height/2.0 + k*dz;
        for (std::size_t j = 0; j < NumPtsY; ++j)
        {
            const double y = Ycom - Width/2.0 + j*dy;
            for (std::size_t i = 0; i < NumPtsX; ++i)
            {
                const double x = Xcom - Length/2.0 + i*dy;
                LagX.push_back(x);
                LagY.push_back(y);
                LagZ.push_back(z);
            }
        }
    }

    std::ofstream out("rectangle3d.vertex");
    out << LagZ.size() << '\n';
    out << std::scientific << std::uppercase;
    for (std::size_t i = 0; i < LagZ.size(); ++i)
    {
        out << std::setprecision(7) << LagX[i] << "\t\t"
            << std::setprecision(7) << LagY[i] << "\t\t"
            << std::setprecision(7) << LagZ[i] << '\n';
    }
}
