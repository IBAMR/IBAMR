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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_StaggeredStokesIBOperator
#define included_IBAMR_StaggeredStokesIBOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

#include <ibtk/GeneralOperator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <IntVector.h>

#include <string>
#include <vector>

namespace IBAMR
{
class IBImplicitStrategy;
class StaggeredStokesOperator;
} // namespace IBAMR
namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace math
{
template <int DIM, class TYPE>
class HierarchyDataOpsReal;
} // namespace math
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!\brief Shared state used by Stokes-IB velocity-path operator components. */
struct StaggeredStokesIBOperatorContext
{
    SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_implicit_ops = nullptr;
    SAMRAI::tbox::Pointer<StaggeredStokesOperator> stokes_op = nullptr;
    IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op = nullptr;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double>> hier_velocity_data_ops = nullptr;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>> u_synch_scheds;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> u_ghost_fill_scheds;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> f_prolongation_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level = nullptr;
    int u_idx = IBTK::invalid_index;
    int f_idx = IBTK::invalid_index;
    int u_current_idx = IBTK::invalid_index;
    int u_dof_index_idx = IBTK::invalid_index;
    int p_dof_index_idx = IBTK::invalid_index;
    bool use_fixed_le_operators = true;
    TimeSteppingType time_stepping_type = MIDPOINT_RULE;
};

/*!\brief Velocity-path nonlinear Stokes-IB residual operator. */
class StaggeredStokesIBOperator : public IBTK::GeneralOperator
{
public:
    StaggeredStokesIBOperator(const std::string& object_name, bool homogeneous_bc = false);

    ~StaggeredStokesIBOperator() override;

    void setOperatorContext(const StaggeredStokesIBOperatorContext& ctx);

    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z) override;

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    void deallocateOperatorState() override;

    void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    void imposeSolBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u) override;

private:
    StaggeredStokesIBOperator() = delete;
    StaggeredStokesIBOperator(const StaggeredStokesIBOperator& from) = delete;
    StaggeredStokesIBOperator& operator=(const StaggeredStokesIBOperator& that) = delete;

    StaggeredStokesIBOperatorContext d_ctx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBOperator
