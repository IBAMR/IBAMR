#ifndef included_IBTK_interpolation_utilities_inc
#define included_IBTK_interpolation_utilities_inc
#include "ibtk/interpolation_utilities.h"

namespace IBTK
{

namespace Interpolation
{
template <typename VectorArray>
IBTK::MatrixXd
formMonomialBasis(const std::vector<VectorArray>& pts, int deg, double ds, const VectorArray& shft)
{
    // Number of polynomials
    size_t num_pts = pts.size();
    int num_poly = getNumberOfPolynomials(deg);
    IBTK::MatrixXd mat = IBTK::MatrixXd::Zero(num_pts, num_poly);
    // Shift and scale monomials
    for (size_t row = 0; row < pts.size(); ++row)
    {
        Point pt = pts[row];
        for (unsigned int d = 0; d < NDIM; ++d) pt(d) = (pt(d) - shft(d)) / ds;
        int col = 0;
#if (NDIM == 2)
        for (int i = 0; i <= deg; ++i)
        {
            for (int j = 0; j <= deg; ++j)
            {
                if ((i + j) <= deg)
                {
                    mat(row, col++) = Interpolation::pow(pt(0), i) * Interpolation::pow(pt(1), j);
                }
            }
        }
#endif
#if (NDIM == 3)
        for (int i = 0; i <= deg; ++i)
        {
            for (int j = i; j >= 0; --j)
            {
                for (int k = i - j; k >= 0; --k)
                {
                    mat(row, col++) = Interpolation::pow(pt(0), j) * Interpolation::pow(pt(1), k) *
                                      Interpolation::pow(pt(2), i - k - j);
                }
            }
        }
#endif
    }
    return mat;
}

inline int
getNumberOfPolynomials(int deg)
{
#if (NDIM == 2)
    return (deg + 1) * (deg + 2) / 2;
#endif
#if (NDIM == 3)
    return (deg + 1) * (deg + 2) * (deg + 3) / 6;
#endif
}

inline double
pow(const double q, const int i)
{
    if (i == 0) return 1.0;
    if (q == 0.0)
        return 0.0;
    else
        return std::pow(q, i);
}
} // namespace Interpolation
} // namespace IBTK

#endif
