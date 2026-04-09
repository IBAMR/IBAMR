// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_INSSGSKinematics
#define included_IBAMR_INSSGSKinematics

#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <EdgeVariable.h>
#include <NodeVariable.h>
#include <SideVariable.h>

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Helper object for reconstructing resolved SGS-model kinematic data
 * from a ghost-filled staggered velocity field.
 */
class INSSGSKinematics : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    explicit INSSGSKinematics(std::string object_name);

    /*!
     * \brief Fill the internal ghosted staggered velocity field.
     */
    void fillGhostedVelocity(int U_idx,
                             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                             double data_time);

    /*!
     * \brief Compute the cell-centered symmetric strain-rate tensor from the
     * internal ghosted staggered velocity field.
     */
    void computeCellCenteredStrainRate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                       double data_time);

#if (NDIM == 2)
    /*!
     * \brief Compute the node-centered symmetric strain-rate tensor directly
     * from the ghost-filled staggered velocity field using node-appropriate
     * finite-difference stencils.
     */
    void computeNodeCenteredStrainRate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                       double data_time);
#endif

#if (NDIM == 3)
    /*!
     * \brief Compute the edge-centered symmetric strain-rate tensor directly
     * from the ghost-filled staggered velocity field using edge-appropriate
     * finite-difference stencils.
     */
    void computeEdgeCenteredStrainRate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                       double data_time);
#endif

    /*!
     * \brief Return the internal ghost-filled staggered velocity index.
     */
    int getGhostedVelocityPatchDataIndex() const;

    /*!
     * \brief Return the internal cell-centered symmetric strain-rate index.
     */
    int getCellCenteredStrainRatePatchDataIndex() const;

#if (NDIM == 2)
    /*!
     * \brief Return the internal node-centered symmetric strain-rate index.
     */
    int getNodeCenteredStrainRatePatchDataIndex() const;
#endif

#if (NDIM == 3)
    /*!
     * \brief Return the internal edge-centered symmetric strain-rate index.
     */
    int getEdgeCenteredStrainRatePatchDataIndex() const;
#endif

private:
    /*!
     * Object name used in subordinate variable and context names.
     */
    std::string d_object_name;

    /*!
     * Scratch side-centered velocity with valid ghost data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_U_scratch_var;

    /*!
     * Cell-centered symmetric strain-rate tensor with depth
     * `NDIM * (NDIM + 1) / 2`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_S_cc_var;

#if (NDIM == 2)
    /*!
     * Node-centered symmetric strain-rate tensor with depth
     * `NDIM * (NDIM + 1) / 2`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_S_nc_var;
#endif

#if (NDIM == 3)
    /*!
     * Edge-centered symmetric strain-rate tensor with depth
     * `NDIM * (NDIM + 1) / 2`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double>> d_S_ec_var;
#endif

    /*!
     * Patch-data descriptor index for d_U_scratch_var.
     */
    int d_U_scratch_idx = IBTK::invalid_index;

    /*!
     * Patch-data descriptor index for d_S_cc_var.
     */
    int d_S_cc_idx = IBTK::invalid_index;

#if (NDIM == 2)
    /*!
     * Patch-data descriptor index for d_S_nc_var.
     */
    int d_S_nc_idx = IBTK::invalid_index;
#endif

#if (NDIM == 3)
    /*!
     * Patch-data descriptor index for d_S_ec_var.
     */
    int d_S_ec_idx = IBTK::invalid_index;
#endif
};
} // namespace IBAMR

#endif
