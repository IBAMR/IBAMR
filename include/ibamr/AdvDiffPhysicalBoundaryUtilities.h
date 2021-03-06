// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_AdvDiffPhysicalBoundaryUtilities
#define included_IBAMR_AdvDiffPhysicalBoundaryUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/Pointer.h"

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
template <int DIM, class TYPE>
class FaceData;
} // namespace pdat
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
 * \brief Class AdvDiffPhysicalBoundaryUtilities is a utility class that
 * provides functions useful for dealing with physical boundary conditions for
 * advection-diffusion solvers.
 */
class AdvDiffPhysicalBoundaryUtilities
{
public:
    /*!
     * \brief Set physical boundary conditions at physical boundaries.
     */
    static void setPhysicalBoundaryConditions(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > Q_data,
                                              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > u_ADV_data,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                              const double fill_time,
                                              const bool inflow_boundaries_only,
                                              const bool homogeneous_bc);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPhysicalBoundaryUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPhysicalBoundaryUtilities(const AdvDiffPhysicalBoundaryUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPhysicalBoundaryUtilities& operator=(const AdvDiffPhysicalBoundaryUtilities& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_AdvDiffPhysicalBoundaryUtilities
