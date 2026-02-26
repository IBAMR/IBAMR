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

#ifndef included_IBTK_StaggeredPhysicalBoundaryHelper
#define included_IBTK_StaggeredPhysicalBoundaryHelper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>
#include <ibtk/samrai_compatibility_names.h>

#include <tbox/DescribedClass.h>

#include <SAMRAIArray.h>
#include <SAMRAIArrayData.h>
#include <SAMRAIBoundaryBox.h>
#include <SAMRAIBox.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRobinBcCoefStrategy.h>
#include <SAMRAISideData.h>

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
                                       int coarsest_ln = invalid_level_number,
                                       int finest_ln = invalid_level_number) const;

    /*!
     * \brief Copy data to u_data_out_idx from u_data_in_idx at Dirichlet
     * boundaries on a single patch.
     */
    void copyDataAtDirichletBoundaries(SAMRAIPointer<SAMRAISideData<double>> u_out_data,
                                       SAMRAIPointer<SAMRAISideData<double>> u_in_data,
                                       SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Setup a masking function over the specified range of levels in the
     * patch hierarchy.
     */
    void setupMaskingFunction(int mask_data_idx,
                              int coarsest_ln = invalid_level_number,
                              int finest_ln = invalid_level_number) const;

    /*!
     * \brief Setup a masking function on a single patch.
     */
    void setupMaskingFunction(SAMRAIPointer<SAMRAISideData<int>> u_data, SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Return a boolean value indicating whether a patch has Dirichlet
     * boundaries.
     */
    bool patchTouchesDirichletBoundary(SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Return a boolean value indicating whether a patch has Dirichlet
     * boundaries in the specified coordinate axis.
     */
    bool patchTouchesDirichletBoundaryAxis(SAMRAIPointer<SAMRAIPatch> patch, const unsigned int axis) const;

    /*!
     * \brief Cache boundary coefficient data.
     */
    void cacheBcCoefData(const std::vector<SAMRAIRobinBcCoefStrategy*>& u_bc_coefs,
                         double fill_time,
                         SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy);

    /*!
     * \brief Clear cached boundary coefficient data.
     */
    virtual void clearBcCoefData();

protected:
    /*!
     * \brief Setup boundary boxes used for setting boundary condition
     * coefficients.
     */
    static void setupBcCoefBoxes(SAMRAIBox& bc_coef_box,
                                 SAMRAIBoundaryBox& trimmed_bdry_box,
                                 const SAMRAIBoundaryBox& bdry_box,
                                 SAMRAIPointer<SAMRAIPatch> patch);

    /*!
     * Cached hierarchy-related information.
     */
    SAMRAIPointer<SAMRAIPatchHierarchy> d_hierarchy;
    std::vector<std::map<int, SAMRAIArray<SAMRAIBoundaryBox>>> d_physical_codim1_boxes;
    std::vector<std::map<int, std::vector<SAMRAIPointer<SAMRAIArrayData<bool>>>>> d_dirichlet_bdry_locs;

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

#endif // #ifndef included_IBTK_StaggeredPhysicalBoundaryHelper
