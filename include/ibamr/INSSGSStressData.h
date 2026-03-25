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

#ifndef included_IBAMR_INSSGSStressData
#define included_IBAMR_INSSGSStressData

#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <EdgeVariable.h>
#include <NodeVariable.h>
#include <SideVariable.h>

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Stores an SGS stress tensor on a staggered layout suitable for the
 * side-centered momentum equation.
 *
 * Diagonal components are cell-centered. In 2D, the single shear component is
 * node-centered. In 3D, the shear components are edge-centered, with edge axis
 * `0/1/2` storing `tau_yz/tau_xz/tau_xy`, respectively.
 */
class INSSGSStressData : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    explicit INSSGSStressData(std::string object_name);

    /*!
     * \brief Allocate the internal staggered stress patch data.
     */
    void allocatePatchData(double data_time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * \brief Deallocate the internal staggered stress patch data.
     */
    void deallocatePatchData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * \brief Zero the internal staggered stress data.
     */
    void setToZero(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * \brief Fill the staggered stress representation from a cell-centered
     * symmetric stress tensor with depth `NDIM * (NDIM + 1) / 2`.
     */
    void setFromCellCenteredStressTensor(int tau_cc_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * \brief Compute the side-centered divergence of the stored staggered
     * stress tensor.
     */
    void computeDivergence(int F_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy) const;

    /*!
     * \brief Return the patch-data descriptor index for the diagonal stress
     * variable.
     */
    int getDiagonalStressPatchDataIndex() const;

    /*!
     * \brief Return the patch-data descriptor index for the shear stress
     * variable.
     */
    int getShearStressPatchDataIndex() const;

private:
    /*!
     * Object name used in subordinate variable and context names.
     */
    std::string d_object_name;

    /*!
     * Cell-centered diagonal stresses `tau_xx`, `tau_yy`, `tau_zz`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_tau_diag_var;

#if (NDIM == 2)
    /*!
     * Node-centered shear stress `tau_xy`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_tau_shear_var;
#endif

#if (NDIM == 3)
    /*!
     * Edge-centered shear stresses. Edge axes `0/1/2` store
     * `tau_yz/tau_xz/tau_xy`.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double>> d_tau_shear_var;
#endif

    /*!
     * Patch-data descriptor index for d_tau_diag_var.
     */
    int d_tau_diag_idx = IBTK::invalid_index;

    /*!
     * Patch-data descriptor index for d_tau_shear_var.
     */
    int d_tau_shear_idx = IBTK::invalid_index;
};
} // namespace IBAMR

#endif
