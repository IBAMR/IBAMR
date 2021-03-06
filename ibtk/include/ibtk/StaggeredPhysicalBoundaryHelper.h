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

#ifndef included_IBTK_StaggeredPhysicalBoundaryHelper
#define included_IBTK_StaggeredPhysicalBoundaryHelper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <map>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Box;
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class ArrayData;
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
template <class TYPE>
class Array;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class StaggeredPhysicalBoundaryHelper provides helper functions to
 * handle physical boundary conditions for a staggered grid discretizations.
 */
class StaggeredPhysicalBoundaryHelper : SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    StaggeredPhysicalBoundaryHelper() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~StaggeredPhysicalBoundaryHelper() = default;

    /*!
     * \brief Copy data to u_data_out_idx from u_data_in_idx at Dirichlet
     * boundaries over the specified range of levels in the patch hierarchy.
     */
    void copyDataAtDirichletBoundaries(int u_out_data_idx,
                                       int u_in_data_idx,
                                       int coarsest_ln = -1,
                                       int finest_ln = -1) const;

    /*!
     * \brief Copy data to u_data_out_idx from u_data_in_idx at Dirichlet
     * boundaries on a single patch.
     */
    void copyDataAtDirichletBoundaries(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > u_out_data,
                                       SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > u_in_data,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Setup a masking function over the specified range of levels in the
     * patch hierarchy.
     */
    void setupMaskingFunction(int mask_data_idx, int coarsest_ln = -1, int finest_ln = -1) const;

    /*!
     * \brief Setup a masking function on a single patch.
     */
    void setupMaskingFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, int> > u_data,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Return a boolean value indicating whether a patch has Dirichlet
     * boundaries.
     */
    bool patchTouchesDirichletBoundary(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Return a boolean value indicating whether a patch has Dirichlet
     * boundaries in the specified coordinate axis.
     */
    bool patchTouchesDirichletBoundaryAxis(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                           const unsigned int axis) const;

    /*!
     * \brief Cache boundary coefficient data.
     */
    void cacheBcCoefData(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                         double fill_time,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Clear cached boundary coefficient data.
     */
    virtual void clearBcCoefData();

protected:
    /*!
     * \brief Setup boundary boxes used for setting boundary condition
     * coefficients.
     */
    static void setupBcCoefBoxes(SAMRAI::hier::Box<NDIM>& bc_coef_box,
                                 SAMRAI::hier::BoundaryBox<NDIM>& trimmed_bdry_box,
                                 const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch);

    /*!
     * Cached hierarchy-related information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    std::vector<std::map<int, SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> > > > d_physical_codim1_boxes;
    std::vector<std::map<int, std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, bool> > > > >
        d_dirichlet_bdry_locs;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredPhysicalBoundaryHelper(const StaggeredPhysicalBoundaryHelper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredPhysicalBoundaryHelper& operator=(const StaggeredPhysicalBoundaryHelper& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_StaggeredPhysicalBoundaryHelper
