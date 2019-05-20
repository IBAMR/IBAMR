#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <vector>

int main()
{
    // mesh parameters:
    const double radius = 0.5;
    const double Lx = 5.0;
    const double Ly = 5.0;
    const double Nx = 200;
    const double dx = Lx/Nx;

    // disk parameters:
    const double h = dx;
    const double X_shift = Lx/2.0;
    const double Y_shift = Ly/2.0;
    // h/16 to put marker away from grid lines on second level
    // The markers are symmetric on single level (by taking odd #cells)
    const double X_com = 0.0;
    const double Y_com = 0.0;

    // write out the points.
    const int Npts = std::floor(radius/h);

    std::vector<double> X_first;
    std::vector<double> Y_first;
    // first quadrant:
    for (int j = 0; j < Npts; ++j)
    {
        const double y = Y_com + (j + 1)*h;
        for (int i = 0; i < Npts; ++i)
        {
            const double x = X_com + (i + 1)*h;
            const double dist2 = std::pow(x - X_com, 2) + std::pow(y - Y_com, 2);
            if (dist2 <= radius*radius)
            {
                X_first.push_back(x);
                Y_first.push_back(y);
            }
        }
    }

    // second quadrant:
    std::vector<double> X_second = X_first;
    std::vector<double> Y_second = Y_first;
    std::transform(X_second.begin(), X_second.end(), X_second.begin(), std::negate<double>{});

    // third quadrant:
    std::vector<double> X_third = X_first;
    std::vector<double> Y_third = Y_first;
    std::transform(X_third.begin(), X_third.end(), X_third.begin(), std::negate<double>{});
    std::transform(Y_third.begin(), Y_third.end(), Y_third.begin(), std::negate<double>{});

    // fourth quadrant:
    std::vector<double> X_fourth = X_first;
    std::vector<double> Y_fourth = Y_first;
    std::transform(Y_fourth.begin(), Y_fourth.end(), Y_fourth.begin(), std::negate<double>{});

    // X axis:
    std::vector<double> X_xaxis;
    std::vector<double> Y_xaxis;
    for (int k = -Npts + 1; k < Npts; ++k)
    {
        X_xaxis.push_back(X_com + k*h);
        Y_xaxis.push_back(X_com);
    }

    // Y axis:
    std::vector<double> X_yaxis;
    std::vector<double> Y_yaxis;
    for (int k = -Npts + 1; k < Npts; ++k)
    {
        X_yaxis.push_back(X_com);
        Y_yaxis.push_back(Y_com + k*h);
    }

    std::vector<std::vector<double> *> xs {&X_first, &X_second, &X_third, &X_fourth, &X_xaxis, &X_yaxis};
    std::vector<std::vector<double> *> ys {&Y_first, &Y_second, &Y_third, &Y_fourth, &Y_xaxis, &Y_yaxis};

    std::ofstream out("cylinder2d.vertex");
    std::size_t total_n_pts = 0;
    for (std::size_t outer = 0; outer < xs.size(); ++outer)
        total_n_pts += xs[outer]->size();
    out << total_n_pts << '\n';
    for (std::size_t outer = 0; outer < xs.size(); ++outer)
    {
        const std::vector<double> &x = *xs[outer];
        const std::vector<double> &y = *ys[outer];
        for (std::size_t inner = 0; inner < x.size(); ++inner)
            out << std::setfill('0') << std::setw(8) << std::left << x[inner] + X_shift << '\t'
                << std::setfill('0') << std::setw(8) << std::left << y[inner] + Y_shift << '\n';
    }

}
