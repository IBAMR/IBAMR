#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

int
main()
{
    const double Lx = 12.0;
    const double Ly = 2.0;
    const double Lz = 3.0;

    // domain parameters (z is y in 2D)
    const std::size_t N = 50 * 2;
    const std::size_t Nx = 6 * N;
    const std::size_t Ny = N;
    const std::size_t Nz = 3 * N / 2;

    // dimensional parameters
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    const double dz = Lz / Nz;

    const double Length = 1.2;
    const double Width = 1.2;
    const double Length_From_Left = Lx / 2.0;
    const double Distance_From_Edge = Ly / 2.0;
    const double Height_Above_Bottom = 2.3;
    const double angle = 25.0; // degrees

    const double x0 = Length_From_Left;
    const double z0 = Height_Above_Bottom;

    const double Height = Length / 2.0 * std::tan(angle * M_PI / 180.0);
    const double m = Height / (Length / 2.0); // slope of xz line
    const double br = z0 - m * x0;
    const double bl = z0 + m * x0;

    const std::size_t NumPtsX = std::ceil(Length / dx) + 1;
    const std::size_t NumPtsY = std::ceil(Width / dy) + 1;
    const std::size_t NumPtsZ = std::ceil(Height / dz) + 1;

    // 2D wedge
    {
        std::ofstream out("wedge2d.vertex");
        std::vector<double> LagX;
        std::vector<double> LagZ;
        for (std::size_t k = 0; k < NumPtsZ; ++k)
        {
            const double z = Height_Above_Bottom + k * dz;
            for (std::size_t i = 0; i < NumPtsX; ++i)
            {
                const double x = Length_From_Left - Length / 2.0 + i * dx;
                const double Zr = m * x + br;
                const double Zl = -m * x + bl;
                if (Zr <= z && Zl <= z)
                {
                    LagX.push_back(x);
                    LagZ.push_back(z);
                }
            }
        }

        out << LagX.size() << '\n';
        for (std::size_t k = 0; k < LagX.size(); ++k)
        {
            out << std::scientific << std::setprecision(7) << std::setfill('0') << LagX[k] << "\t\t" << LagZ[k] << "\n";
        }
    }

    // 3D wedge
    {
        std::ofstream out("wedge3d.vertex");
        std::vector<double> LagX;
        std::vector<double> LagY;
        std::vector<double> LagZ;
        for (std::size_t k = 0; k < NumPtsZ; ++k)
        {
            const double z = Height_Above_Bottom + k * dz;
            for (std::size_t j = 0; j < NumPtsY; ++j)
            {
                const double y = Distance_From_Edge - Width / 2 + j * dy;
                for (std::size_t i = 0; i < NumPtsX; ++i)
                {
                    const double x = Length_From_Left - Length / 2.0 + i * dx;
                    const double Zr = m * x + br;
                    const double Zl = -m * x + bl;
                    if (Zr <= z && Zl <= z)
                    {
                        LagX.push_back(x);
                        LagY.push_back(y);
                        LagZ.push_back(z);
                    }
                }
            }
        }

        out << LagX.size() << '\n';
        for (std::size_t k = 0; k < LagX.size(); ++k)
        {
            out << std::scientific << std::setprecision(7) << std::setfill('0') << LagX[k] << "\t\t" << LagY[k]
                << "\t\t" << LagZ[k] << "\n";
        }
    }
}
