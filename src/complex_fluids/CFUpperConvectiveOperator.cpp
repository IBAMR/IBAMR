// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/CFUpperConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "FaceData.h"
#include "FaceIndex.h"
#include "Index.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Utilities.h"

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

#if (NDIM == 2)
#define UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(upper_convective_op2d, UPPER_CONVECTIVE_OP2D)
#define SQRT_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(sqrt_upper_convective_op2d, SQRT_UPPER_CONVECTIVE_OP2D)
#define LOG_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(log_upper_convective_op2d, LOG_UPPER_CONVECTIVE_OP2D)
#endif
#if (NDIM == 3)
#define UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(upper_convective_op3d, UPPER_CONVECTIVE_OP3D)
#define SQRT_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(sqrt_upper_convective_op3d, SQRT_UPPER_CONVECTIVE_OP3D)
#define LOG_UPPER_CONVECTIVE_OP_FC IBAMR_FC_FUNC_(log_upper_convective_op3d, LOG_UPPER_CONVECTIVE_OP3D)
#endif

extern "C"
{
#if (NDIM == 2)
    void UPPER_CONVECTIVE_OP_FC(const double*,
                                const double*,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&);
    void SQRT_UPPER_CONVECTIVE_OP_FC(const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void LOG_UPPER_CONVECTIVE_OP_FC(const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
#if (NDIM == 3)
    void UPPER_CONVECTIVE_OP_FC(const double*,
                                const double*,
                                const double*,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const double*,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&,
                                const int&);
    void SQRT_UPPER_CONVECTIVE_OP_FC(const double*,
                                     const double*,
                                     const double*,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const double*,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&,
                                     const int&);
    void LOG_UPPER_CONVECTIVE_OP_FC(const double*,
                                    const double*,
                                    const double*,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const double*,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&,
                                    const int&);
#endif
}

namespace IBAMR
{
CFUpperConvectiveOperator::CFUpperConvectiveOperator(const std::string& object_name,
                                                     Pointer<CellVariable<NDIM, double> > Q_var,
                                                     Pointer<Database> input_db,
                                                     const std::string& difference_form,
                                                     const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
                                                     const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
    : ConvectiveOperator(object_name, UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE),
      d_Q_var(Q_var),
      d_u_adv_var(new SideVariable<NDIM, double>("Complex U var")),
      d_difference_form(difference_form),
      d_Q_bc_coefs(Q_bc_coefs),
      d_u_bc_coefs(u_bc_coefs)
{
    if (input_db)
    {
        d_interp_type = input_db->getStringWithDefault("interp_type", d_interp_type);
        if (input_db->keyExists("evolution_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolution_type"));
        if (input_db->keyExists("evolve_type"))
            d_evolve_type = string_to_enum<TensorEvolutionType>(input_db->getString("evolve_type"));
    }
    // Register some scratch variables
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_Q_convec_idx = var_db->registerVariableAndContext(
        d_Q_var, var_db->getContext(d_object_name + "::CONVECTIVE"), IntVector<NDIM>(0));
    Pointer<VariableContext> new_cxt = var_db->getContext(d_object_name + "::U_ADV_CXT");
    d_u_scratch_idx = var_db->registerVariableAndContext(d_u_adv_var, new_cxt, IntVector<NDIM>(2));
    Pointer<VariableContext> src_cxt = var_db->getContext(d_object_name + "::SOURCE");
    d_s_idx = var_db->registerVariableAndContext(d_Q_var, src_cxt);
    auto convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
    d_convec_oper = convective_op_manager->allocateOperator(
        d_difference_form, d_object_name + "::convec_oper", d_Q_var, input_db, ADVECTIVE, d_Q_bc_coefs);
} // Constructor

CFUpperConvectiveOperator::~CFUpperConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // Destructor

void
CFUpperConvectiveOperator::applyConvectiveOperator(int Q_idx, int Y_idx)
{
    if (!d_s_fcn)
    {
        TBOX_ERROR("CFUpperConvectiveOperator::applyConvectiveOperator():\n"
                   << "  Source function must be register prior to call to "
                      "applyConvectiveOperator\n");
    }
    if (!d_is_initialized)
    {
        TBOX_ERROR("CFUpperConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
    // Set up velocity information:
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    d_convec_oper->setSolutionTime(d_solution_time);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    // Copy data to side centered velocity field
    //
    // TODO: This is only done because we currently only have operators to
    // fill in physical boundary conditions for cell- and side-centered
    // quantities.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            Pointer<FaceData<NDIM, double> > u_f_data = patch->getPatchData(d_u_idx);
            Pointer<SideData<NDIM, double> > u_s_data = patch->getPatchData(d_u_scratch_idx);
            const Box<NDIM> box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> i(box, axis); i; i++)
                {
                    const SideIndex<NDIM>& si = i();
                    FaceIndex<NDIM> fi(si.toCell(0), axis, 1);
                    (*u_s_data)(si) = (*u_f_data)(fi);
                }
            }
        }
    }
    // Fill boundary conditions for side centered velocity
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        d_u_scratch_idx, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, d_u_bc_coefs, NULL, d_interp_type);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, d_hierarchy);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_u_bc_coefs, nullptr, d_u_scratch_idx, -1, false);
    ghost_fill_op.fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);

    d_convec_oper->applyConvectiveOperator(Q_idx, d_Q_convec_idx);

    d_s_fcn->setPatchDataIndex(Q_idx);
    d_s_fcn->setDataOnPatchHierarchy(d_s_idx, d_Q_var, d_hierarchy, d_solution_time, false, d_coarsest_ln, d_finest_ln);

    for (int level_num = d_coarsest_ln; level_num <= d_finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_scratch_idx);
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* dx = p_geom->getDx();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM> patch_lower = patch_box.lower();
            const IntVector<NDIM> patch_upper = patch_box.upper();
            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double> > Y_data = patch->getPatchData(Y_idx);
            const IntVector<NDIM> Q_data_gcw = Q_data->getGhostCellWidth();

            const IntVector<NDIM> u_data_gcw = u_data->getGhostCellWidth();
            const IntVector<NDIM> Y_data_gcw = Y_data->getGhostCellWidth();
            Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_Q_convec_idx);
            const IntVector<NDIM> C_data_gcw = C_data->getGhostCellWidth();
            Pointer<CellData<NDIM, double> > S_data = patch->getPatchData(d_s_idx);
            const IntVector<NDIM> S_data_gcw = S_data->getGhostCellWidth();

            switch (d_evolve_type)
            {
            case SQUARE_ROOT:
#if (NDIM == 2)
                SQRT_UPPER_CONVECTIVE_OP_FC(dx,
                                            u_data->getPointer(0),
                                            u_data->getPointer(1),
                                            u_data_gcw.max(),
                                            Q_data->getPointer(0),
                                            Q_data_gcw.max(),
                                            S_data->getPointer(0),
                                            S_data_gcw.max(),
                                            C_data->getPointer(0),
                                            C_data_gcw.max(),
                                            Y_data->getPointer(0),
                                            Y_data_gcw.max(),
                                            patch_lower(0),
                                            patch_upper(0),
                                            patch_lower(1),
                                            patch_upper(1));
#endif
#if (NDIM == 3)
                SQRT_UPPER_CONVECTIVE_OP_FC(dx,
                                            u_data->getPointer(0),
                                            u_data->getPointer(1),
                                            u_data->getPointer(2),
                                            u_data_gcw.max(),
                                            Q_data->getPointer(0),
                                            Q_data_gcw.max(),
                                            S_data->getPointer(0),
                                            S_data_gcw.max(),
                                            C_data->getPointer(0),
                                            C_data_gcw.max(),
                                            Y_data->getPointer(0),
                                            Y_data_gcw.max(),
                                            patch_lower(0),
                                            patch_upper(0),
                                            patch_lower(1),
                                            patch_upper(1),
                                            patch_lower(2),
                                            patch_upper(2));
#endif
                break;
            case LOGARITHM:
#if (NDIM == 2)
                LOG_UPPER_CONVECTIVE_OP_FC(dx,
                                           u_data->getPointer(0),
                                           u_data->getPointer(1),
                                           u_data_gcw.max(),
                                           Q_data->getPointer(0),
                                           Q_data_gcw.max(),
                                           S_data->getPointer(0),
                                           S_data_gcw.max(),
                                           C_data->getPointer(0),
                                           C_data_gcw.max(),
                                           Y_data->getPointer(0),
                                           Y_data_gcw.max(),
                                           patch_lower(0),
                                           patch_upper(0),
                                           patch_lower(1),
                                           patch_upper(1));
#endif
#if (NDIM == 3)
                LOG_UPPER_CONVECTIVE_OP_FC(dx,
                                           u_data->getPointer(0),
                                           u_data->getPointer(1),
                                           u_data->getPointer(2),
                                           u_data_gcw.max(),
                                           Q_data->getPointer(0),
                                           Q_data_gcw.max(),
                                           S_data->getPointer(0),
                                           S_data_gcw.max(),
                                           C_data->getPointer(0),
                                           C_data_gcw.max(),
                                           Y_data->getPointer(0),
                                           Y_data_gcw.max(),
                                           patch_lower(0),
                                           patch_upper(0),
                                           patch_lower(1),
                                           patch_upper(1),
                                           patch_lower(2),
                                           patch_upper(2));
#endif
                break;
            case STANDARD:
#if (NDIM == 2)
                UPPER_CONVECTIVE_OP_FC(dx,
                                       u_data->getPointer(0),
                                       u_data->getPointer(1),
                                       u_data_gcw.max(),
                                       Q_data->getPointer(0),
                                       Q_data_gcw.max(),
                                       S_data->getPointer(0),
                                       S_data_gcw.max(),
                                       C_data->getPointer(0),
                                       C_data_gcw.max(),
                                       Y_data->getPointer(0),
                                       Y_data_gcw.max(),
                                       patch_lower(0),
                                       patch_upper(0),
                                       patch_lower(1),
                                       patch_upper(1));
#endif
#if (NDIM == 3)
                UPPER_CONVECTIVE_OP_FC(dx,
                                       u_data->getPointer(0),
                                       u_data->getPointer(1),
                                       u_data->getPointer(2),
                                       u_data_gcw.max(),
                                       Q_data->getPointer(0),
                                       Q_data_gcw.max(),
                                       S_data->getPointer(0),
                                       S_data_gcw.max(),
                                       C_data->getPointer(0),
                                       C_data_gcw.max(),
                                       Y_data->getPointer(0),
                                       Y_data_gcw.max(),
                                       patch_lower(0),
                                       patch_upper(0),
                                       patch_lower(1),
                                       patch_upper(1),
                                       patch_lower(2),
                                       patch_upper(2));
#endif
                break;
            case UNKNOWN_TENSOR_EVOLUTION_TYPE:
                TBOX_ERROR(d_object_name << ":\n"
                                         << "Unknown tensor evolution type.\n");
                break;
            }
        } // end Patch loop
    }     // end Level loop
} // applyConvectiveOperator

void
CFUpperConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                   const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();
    // Get Hierarchy Information
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif
    d_convec_oper->initializeOperatorState(in, out);
    d_convec_oper->setAdvectionVelocity(d_u_idx);
    // Allocate Patch Data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_scratch_idx)) level->allocatePatchData(d_u_scratch_idx);
        if (!level->checkAllocated(d_Q_convec_idx)) level->allocatePatchData(d_Q_convec_idx);
        if (!level->checkAllocated(d_s_idx)) level->allocatePatchData(d_s_idx);
    }
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
CFUpperConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    d_convec_oper->deallocateOperatorState();
    // Deallocate scratch data
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Q_convec_idx)) level->deallocatePatchData(d_Q_convec_idx);
        if (level->checkAllocated(d_u_scratch_idx)) level->deallocatePatchData(d_u_scratch_idx);
        if (level->checkAllocated(d_s_idx)) level->deallocatePatchData(d_s_idx);
    }
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
CFUpperConvectiveOperator::registerSourceFunction(Pointer<CFRelaxationOperator> source_fcn)
{
    d_s_fcn = source_fcn;
}

} // namespace IBAMR
