// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_AcousticStreamingPETScMatUtilities
#define included_IBAMR_AcousticStreamingPETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "IntVector.h"
#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include "petscao.h"
#include "petscmat.h"

#include <set>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
template <int DIM>
class CoarseFineBoundary;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AcousticStreamingPETScMatUtilities provides utility functions for
 * <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat objects used in solving acoustic
 * streaming equations.
 *
 * see IBAMR::AcousticStreamingHierarchyIntegrator and IBAMR::FOAcousticStreamingPETScLevelSolver
 */
class AcousticStreamingPETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to a MAC
     * discretization of the first order acoustic streaming equations on a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelFOAcousticStreamingOp(
        Mat& mat,
        double omega,
        double sound_speed,
        int rho_idx,
        int mu_idx,
        int lambda_idx,
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
        double data_time,
        const std::vector<int>& num_dofs_per_proc,
        int u_dof_index_idx,
        int p_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AcousticStreamingPETScMatUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AcousticStreamingPETScMatUtilities(const AcousticStreamingPETScMatUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AcousticStreamingPETScMatUtilities& operator=(const AcousticStreamingPETScMatUtilities& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AcousticStreamingPETScMatUtilities
