#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <vector>

int main()
{
    // mesh parameters:
    const double Lx = 5.0;
    const double Ly = 2.5;
    const double Nx = 500*2*2;
    const double Ny = 250*2*2;
    const double dx = Lx/Nx;
    const double dy = Ly/Ny;

    const double Length = 0.3;
    const double Width = Length/3.0;
    const double theta = 15.0;
    const double t = theta*M_PI/180.0;

    const int NumPtsX = std::ceil(Length/dx) + 1;
    const int NumPtsY = std::ceil(Width/dy) + 1;

    const double X_com = 0.0;
    const double Y_com = 0.0;

    std::vector<double> LagX;
    std::vector<double> LagY;

    // Generate rectangular prism:
    for (int j = 0; j < NumPtsY; ++j)
    {
        const double y = Y_com - Width/2.0 + j*dy;
        for (int i = 0; i < NumPtsX; ++i)
        {
            const double x = X_com - Length/2.0 + i*dy;
            LagX.push_back(x * std::cos(t) - y* std::sin(t));
            LagY.push_back(x * std::sin(t) + y* std::cos(t));
        }
    }

    std::ofstream out("barge2d.vertex");
    out << LagX.size() << '\n';
    out << std::scientific << std::uppercase;
    for (unsigned int index = 0; index < LagX.size(); ++index)
    {
        out << std::setprecision(7) << LagX[index] << "\t\t"
            << std::setprecision(7) << LagY[index] << '\n';
    }

}
