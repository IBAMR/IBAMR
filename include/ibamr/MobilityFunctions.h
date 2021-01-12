// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_MobilityFunctions
#define included_IBAMR_MobilityFunctions

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

namespace IBAMR
{
/*!
 * \brief Class MobilityFunctions provides empirical functions/fits for
 * forming mobility matrix needed for fully constraint based IB method.
 */
class MobilityFunctions
{
public:
    /*!
     * \brief Construct the mobility matrix from empirical fits.
     *
     * \param kernel_name IB kernel function.
     * \note Supported IB kernels are "IB_3", "IB_4" and "IB_6".
     *
     * \param mu Fluid density.
     *
     * \param rho Fluid density.
     *
     * \param dt The time step size.
     *
     * \param dx Cartesian grid spacing.
     *
     * \param X Array of IB markers' location.
     *
     * \param num_nodes Number of Lagrangian markers.
     *
     * \param reset_constants Boolean indicating whether all constants
     * are to be reset if beta (the viscous CFL number) is changed (otherwise
     * will use the previous beta for fitting formula).
     *
     * \param periodic_correction Input parameter for incorporating
     * periodic domain correction. Set it to zero if not needed.
     *
     * \param l_domain Length of the domain.
     * \note This is used only only for 2 dimensional steady stokes.
     *
     * \param mm Pointer to mobility matrix.
     */
    static void constructEmpiricalMobilityMatrix(const char* kernel_name,
                                                 const double mu,
                                                 const double rho,
                                                 const double dt,
                                                 const double dx,
                                                 const double* X,
                                                 const int num_nodes,
                                                 const int reset_constants,
                                                 const double periodic_correction,
                                                 const double l_domain,
                                                 double* mm);

    /*!
     * \brief Construct the mobility matrix from Rotne-Pragner-Yamakawa tensor.
     *
     * \param kernel_name IB kernel function.
     * \note Supported IB kernels are "IB_3", "IB_4" and "IB_6".
     *
     * \param mu Fluid density.
     *
     * \param dx Cartesian grid spacing.
     *
     * \param X Array of IB markers' location.
     *
     * \param num_nodes Number of Lagrangian markers.
     *
     * \param periodic_correction Input parameter for incorporating
     * periodic domain correction. Set it to zero if not known.
     *
     * \param l_domain Length of the domain.
     *
     * \param mm Pointer to mobility matrix.
     */
    static void constructRPYMobilityMatrix(const char* kernel_name,
                                           const double mu,
                                           const double dx,
                                           const double* X,
                                           const int num_nodes,
                                           const double periodic_correction,
                                           double* mm);
}; // MobilityFunctions

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // ifndef included_IBAMR_MobilityFunctions
