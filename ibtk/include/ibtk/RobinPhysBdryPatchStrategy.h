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

#ifndef included_IBTK_RobinPhysBdryPatchStrategy
#define included_IBTK_RobinPhysBdryPatchStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"

#include <set>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class RobinPhysBdryPatchStrategy is an abstract strategy class that
 * extends the SAMRAI::xfer::RefinePatchStrategy interface to provide routines
 * specific for setting Robin boundary conditions at physical boundaries.  This
 * class also provides default implementations of some methods defined in
 * SAMRAI::xfer::RefinePatchStrategy that are generally not needed for filling
 * ghost cell values at physical boundaries.
 */
class RobinPhysBdryPatchStrategy : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    RobinPhysBdryPatchStrategy() = default;

    /*!
     * \brief Destructor.
     */
    ~RobinPhysBdryPatchStrategy() = default;

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    void setPatchDataIndex(int patch_data_index);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const std::set<int>& patch_data_indices);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices);

    /*!
     * \brief Reset the Robin boundary condition specification object employed
     * by this class to set physical boundary conditions.
     *
     * \note \a bc_coef cannot be NULL.
     */
    void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

    /*!
     * \brief Reset the Robin boundary condition specification object employed
     * by this class to set physical boundary conditions.
     *
     * \note None of the elements of \a bc_coefs can be NULL.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Set whether boundary filling should employ homogeneous boundary
     * conditions.
     *
     * \note By default, inhomogeneous boundary conditions are assumed.
     */
    void setHomogeneousBc(bool homogeneous_bc);

    /*!
     * \return Whether boundary filling employs homogeneous boundary conditions.
     */
    bool getHomogeneousBc() const;

    /*!
     * \name Partial implementation of SAMRAI::xfer::RefinePatchStrategy
     * interface.
     */
    //\{

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The preprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * The default implementation does nothing.  This behavior can be overridden by
     * subclasses.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void preprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                          const SAMRAI::hier::Patch<NDIM>& coarse,
                          const SAMRAI::hier::Box<NDIM>& fine_box,
                          const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    /*!
     * Function to perform user-defined postprocess data refine operations.
     * This member function is called after standard refine operations
     * (expressed using concrete subclasses of the SAMRAI::xfer::RefineOperator
     * base class).  The postprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * The implementation does nothing.  This behavior can be overridden by
     * subclasses.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void postprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio) override;

    //\}

    /*!
     * Function to accumulate data near physical boundaries from values set in
     * the ghost cell region using the adjoint of the operator used to
     * extrapolate the ghost cell values.  This function can be used to
     * construct the adjoint of linear operators that use ghost cell data.
     *
     * \note A default implementation is provided that emits an error message.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    virtual void accumulateFromPhysicalBoundaryData(SAMRAI::hier::Patch<NDIM>& patch,
                                                    double fill_time,
                                                    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

protected:
    /*
     * The patch data indices corresponding to the "scratch" patch data that
     * requires extrapolation of ghost cell values at physical boundaries.
     */
    std::set<int> d_patch_data_indices;

    /*
     * The RobinBcCoefStrategy objects used to specify Robin boundary conditions
     * for each data depth.
     *
     * The boolean value indicates whether homogeneous boundary conditions
     * should be used.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    bool d_homogeneous_bc = false;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RobinPhysBdryPatchStrategy(const RobinPhysBdryPatchStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RobinPhysBdryPatchStrategy& operator=(const RobinPhysBdryPatchStrategy& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_RobinPhysBdryPatchStrategy
