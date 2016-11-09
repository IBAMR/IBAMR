// Filename: MobilityFunctions.h
// Created on 17 Feb 2016 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_MobilityFunctions
#define included_IBAMR_MobilityFunctions

/////////////////////////////// INCLUDES /////////////////////////////////////

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
