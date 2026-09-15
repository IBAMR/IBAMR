// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_matrix_free_fixture
#define included_matrix_free_fixture

#include <ibtk/config.h>

#include <Patch.h>
#include <SideData.h>

#include <string>
#include <vector>

namespace MatrixFreeTest
{
/*! \brief Make a single patch with shifted indices and anisotropic spacing. */
SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> make_patch(int cells);

/*! \brief Fill every side including ghosts with a smooth nonconstant field. */
void fill_field(SAMRAI::pdat::SideData<NDIM, double>& field);

/*! \brief Make deterministic interior marker positions in cell or shuffled order. */
std::vector<double> make_positions(const SAMRAI::hier::Patch<NDIM>& patch, int count, bool shuffle);

/*! \brief Evaluate a scalar reference at an arbitrary grid-to-marker distance. */
double reference_weight(const std::string& kernel, int axis, int direction, double distance);

using FortranInterpolate = void(const double*,
                                const double*,
                                const double*,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
#if (NDIM == 3)
                                const int&,
                                const int&,
#endif
                                const int&,
                                const int&
#if (NDIM == 3)
                                ,
                                const int&
#endif
                                ,
                                const double*,
                                const int*,
                                const double*,
                                const int&,
                                const double*,
                                double*);
using FortranSpread = void(const double*,
                           const double*,
                           const double*,
                           const int&,
                           const int*,
                           const double*,
                           const int&,
                           const double*,
                           const double*,
                           const int&,
                           const int&,
                           const int&,
                           const int&,
#if (NDIM == 3)
                           const int&,
                           const int&,
#endif
                           const int&,
                           const int&
#if (NDIM == 3)
                           ,
                           const int&
#endif
                           ,
                           double*);

/*! \brief Select a numerical loop before starting the timer. */
FortranInterpolate* get_fortran_interpolate(const std::string& kernel);

/*! \brief Select a numerical loop before starting the timer. */
FortranSpread* get_fortran_spread(const std::string& kernel);
} // namespace MatrixFreeTest
#endif
