#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/VOFSurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_vof_utilities.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <utility>

#include <ibamr/namespaces.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC_(mollify_ib_4_2d, MOLLIFY_IB_4_2D)
#define HF_CURVATURE_FC IBAMR_FC_FUNC_(hf_curvature_2d, HF_CURVATURE_2D)
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_2d, SC_NORMAL_2D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC_(cc_curvature_2d, CC_CURVATURE_2D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC_(sc_surface_tension_force_2d, SC_SURFACE_TENSION_FORCE_2D)
#define SC_SURFACE_TENSION_FORCE_VOF_FC IBAMR_FC_FUNC_(sc_surface_tension_force_vof_2d, SC_SURFACE_TENSION_FORCE_VOF_2D)

#endif

#if (NDIM == 3)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC_(mollify_ib_4_3d, MOLLIFY_IB_4_3D)
#define HF_CURVATURE_FC IBAMR_FC_FUNC_(hf_curvature_2d, HF_CURVATURE_2D) // Add 3D version of HF_CURVATURE_FC if needed
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_3d, SC_NORMAL_3D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC_(cc_curvature_3d, CC_CURVATURE_3D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC_(sc_surface_tension_force_3d, SC_SURFACE_TENSION_FORCE_3D)
#endif

extern "C"
{
    void MOLLIFY_IB_4_FC(double* V,
        const int& V_gcw,
        const double* U,
        const int& U_gcw,
        const int& ilower0,
        const int& iupper0,
        const int& ilower1,
        const int& iupper1
#if (NDIM == 3)
        ,
        const int& ilower2,
        const int& iupper2
#endif
);

void CC_CURVATURE_FC(double* K,
    const int& K_gcw,
    const double* N00,
    const double* N01,
#if (NDIM == 3)
    const double* N02,
#endif
    const double* N10,
    const double* N11,
#if (NDIM == 3)
    const double* N12,
    const double* N20,
    const double* N21,
    const double* N22,
#endif
    const int& N_gcw,
    const int& ilower0,
    const int& iupper0,
    const int& ilower1,
    const int& iupper1,
#if (NDIM == 3)
    const int& ilower2,
    const int& iupper2,
#endif
    const double* dx);

void SC_SURFACE_TENSION_FORCE_FC(
#if (NDIM == 2)
double* F0,
double* F1,
const int& F_gcw,
const double* K,
const int& K_gcw,
const double* N00,
const double* N11,
const int& N_gcw,
const int& ilower0,
const int& iupper0,
const int& ilower1,
const int& iupper1
#endif
#if (NDIM == 3)
double* F0,
double* F1,
double* F2,
const int& F_gcw,
const double* K,
const int& K_gcw,
const double* N00,
const double* N11,
const double* N22,
const int& N_gcw,
const int& ilower0,
const int& iupper0,
const int& ilower1,
const int& iupper1,
const int& ilower2,
const int& iupper2
#endif
);

    void HF_CURVATURE_FC(double* K,
                         double* K_valid,
                         const int& K_gcw,
                         const double* alpha,
                         const int& alpha_gcw,
                         const int& ilower0,
                         const int& iupper0,
                         const int& ilower1,
                         const int& iupper1,

#if (NDIM == 3)
                         const int& ilower2,
                         const int& iupper2,
#endif

                         const double* dx,
                         const int& hf_radius, // stencil size
                         const double& alpha_tol,
                         const double& gradient_tol);

    void SC_NORMAL_FC(double* N00,
                      double* N01,
#if (NDIM == 3)
                      double* N02,
#endif
                      double* N10,
                      double* N11,
#if (NDIM == 3)
                      double* N12,
                      double* N20,
                      double* N21,
                      double* N22,
#endif
                      const int& N_gcw,
                      const double* U,
                      const int& U_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1,
#if (NDIM == 3)
                      const int& ilower2,
                      const int& iupper2,
#endif
                      const double* dx);

                      void SC_SURFACE_TENSION_FORCE_VOF_FC(
                        double* F0,
                        double* F1,
                        const int& F_gcw,
                        const double* K,
                        const double* K_valid,
                        const int& K_gcw,
                        const double* C,
                        const int& C_gcw,
                        const double* N00,
                        const double* N11,
                        const int& N_gcw,
                        const int& ilower0,
                        const int& iupper0,
                        const int& ilower1,
                        const int& iupper1);
    
}
// Okay, this is okay for a prototype
namespace IBAMR
{
VOFSurfaceTensionForceFunction::VOFSurfaceTensionForceFunction(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    const AdvDiffHierarchyIntegrator* adv_diff_solver,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> level_set_var,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var)
    : SurfaceTensionForceFunction(object_name, input_db, adv_diff_solver, level_set_var), d_vof_var(vof_var)
{
    if (!d_ls_var.isNull() && !d_vof_var.isNull())
    {
        TBOX_ERROR(
            "VOFSurfaceTensionForceFunction: Both level set and VOF variables are provided. Please provide only one.");
    }
    else if (d_ls_var.isNull() && d_vof_var.isNull())
    {
        TBOX_ERROR(
            "VOFSurfaceTensionForceFunction: Neither level set nor VOF variable is provided. Please provide one.");
    }


    if (input_db)
    {
        if (input_db->keyExists("time_stepping_type"))
        {
            d_ts_type = string_to_enum<TimeSteppingType>(input_db->getString("time_stepping_type"));
        }

        // d_kernel_fcn = input_db->getStringWithDefault("kernel", d_kernel_fcn);
        // d_kernel_fcn = input_db->getStringWithDefault("smoother", d_kernel_fcn);
        // d_kernel_fcn = input_db->getStringWithDefault("kernel_fcn", d_kernel_fcn);
        // d_kernel_fcn = input_db->getStringWithDefault("smoother_fcn", d_kernel_fcn);

        d_sigma = input_db->getDoubleWithDefault("sigma", d_sigma);
        d_sigma = input_db->getDoubleWithDefault("surface_tension_coef", d_sigma);

        d_hf_column_half_width = input_db->getIntegerWithDefault("hf_column_half_width", d_hf_column_half_width);
        d_pure_cell_tol = input_db->getDoubleWithDefault("pure_cell_tol", d_pure_cell_tol);
        d_gradient_tol = input_db->getDoubleWithDefault("gradient_tol", d_gradient_tol);
    }


#if (NDIM != 2)
    TBOX_ERROR("VOFSurfaceTensionForceFunction is only implemented for 2D problems.");
#endif

VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

d_alpha_scratch_var = new CellVariable<NDIM, double>(d_object_name + "::alpha_scratch", 1);
d_alpha_idx = var_db->registerVariableAndContext(d_alpha_scratch_var,
                                          var_db->getContext(d_object_name + "::alpha_scratch"),
                                          IntVector<NDIM>(getAlphaGhostWidth()));

d_alpha_smooth_var = new CellVariable<NDIM, double>(d_object_name + "::alpha_smooth", 1);
d_alpha_smooth_idx = var_db->registerVariableAndContext(d_alpha_smooth_var,
                                        var_db->getContext(d_object_name + "::alpha_smooth"),
                                        IntVector<NDIM>(getAlphaGhostWidth()));
if (!d_ls_var.isNull())
{
    Pointer<CellVariable<NDIM,double>> phi_cc_var = d_ls_var; 
    TBOX_ASSERT(!phi_cc_var.isNull());
    d_phi_idx = var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::Phi"), IntVector<NDIM>(getPhiGhostWidth()));
}


    // intentionally blank
    return;
}

int
VOFSurfaceTensionForceFunction::getPhiGhostWidth() const
{

    return d_hf_column_half_width + 2;// +2 for the gradient computation ( since alpha is computed from phi). Can also make it return getAlphaGhostWidth() + 1;
}

int
VOFSurfaceTensionForceFunction::getAlphaGhostWidth() const
{
    return d_hf_column_half_width + 1;
}

bool
VOFSurfaceTensionForceFunction::isTimeDependent() const
{
    return true;
}

void
VOFSurfaceTensionForceFunction::computeAlphaFromLevelSetOnPatch(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch)
{
    Pointer<CellData<NDIM, double>> alpha = patch->getPatchData(d_alpha_idx);
    Pointer<CellData<NDIM, double>> Phi = patch->getPatchData(d_phi_idx);

#if !defined(NDEBUG)
    TBOX_ASSERT(Phi);
    TBOX_ASSERT(alpha);
#endif

#if (NDIM == 2)

    const Box<NDIM>& alpha_box = alpha->getGhostBox();

    for (Box<NDIM>::Iterator it(alpha_box); it; it++)

    {
        const SAMRAI::pdat::CellIndex<NDIM> ci(it());
        const double phi = (*Phi)(ci);

        std::array<double, NDIM> dPhi{}; 
        SAMRAI::pdat::CellIndex<NDIM> ci_xplus(ci);
    SAMRAI::pdat::CellIndex<NDIM> ci_xminus(ci);
    SAMRAI::pdat::CellIndex<NDIM> ci_yplus(ci);
    SAMRAI::pdat::CellIndex<NDIM> ci_yminus(ci);

    ci_xplus(0) += 1;
ci_xminus(0) -= 1;

ci_yplus(1) += 1;
ci_yminus(1) -= 1;
        // This is strictly for 2D, it is okay for testing but not final implementation
        dPhi[0] = 0.5 * ((*Phi)(ci_xplus) - (*Phi)(ci_xminus)); // xplus, xminus
        dPhi[1] = 0.5 * ((*Phi)(ci_yplus) - (*Phi)(ci_yminus));

        (*alpha)(ci) = VCINSVOFUtilities::vof_alpha(phi, dPhi);
    }

#else
    TBOX_ERROR(d_object_name << "::computeAlphaFromLevelSetOnPatch() is only implemented for 2D problems.");
#endif
}
void
VOFSurfaceTensionForceFunction::setDataOnPatchHierarchy(
    int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
    double data_time,
    bool initial_time,
    int coarsest_ln_in,
    int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif

    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    d_hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy, coarsest_ln, finest_ln);

    // Pointer<CellVariable<NDIM, double>> phi_cc_var = d_ls_var; // Common for both classes

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    // using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    RobinBcCoefStrategy<NDIM>* interface_bc_coef = nullptr;

    Pointer<CellVariable<NDIM, double>> alpha_cc_var = d_alpha_scratch_var;

#if !defined(NDEBUG)
    TBOX_ASSERT(!alpha_cc_var.isNull());
#endif

    const IntVector<NDIM> alpha_ghosts(getAlphaGhostWidth());

    // d_alpha_idx =
    //     var_db->registerVariableAndContext(alpha_cc_var, var_db->getContext(d_object_name + "::alpha"), alpha_ghosts);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_alpha_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_alpha_smooth_idx, data_time);
    }

    if (!d_ls_var.isNull())
    {
        Pointer<CellVariable<NDIM, double>> phi_cc_var = d_ls_var;

#if !defined(NDEBUG)
        TBOX_ASSERT(!phi_cc_var.isNull());
#endif

        const int phi_new_idx = var_db->mapVariableAndContextToIndex(phi_cc_var, d_adv_diff_solver->getNewContext());
        const int phi_current_idx =
            var_db->mapVariableAndContextToIndex(phi_cc_var, d_adv_diff_solver->getCurrentContext());

#if !defined(NDEBUG)
        TBOX_ASSERT(phi_new_idx >= 0);
        TBOX_ASSERT(phi_current_idx >= 0);
#endif

        const IntVector<NDIM> cell_ghosts(getPhiGhostWidth());

        // d_phi_idx =
        //     var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::Phi"), cell_ghosts);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(d_phi_idx, data_time);
        }

        // Copy level set into phi and C.
        if (d_ts_type == MIDPOINT_RULE)
        {
            hier_cc_data_ops.linearSum(d_phi_idx,
                                       0.5,
                                       phi_new_idx,
                                       0.5,
                                       phi_current_idx,
                                       /*interior_only*/ true);
        }
        else if (d_ts_type == BACKWARD_EULER)
        {
            hier_cc_data_ops.copyData(d_phi_idx,
                                      phi_new_idx,
                                      /*interior_only*/ true);
        }
        else if (d_ts_type == FORWARD_EULER)
        {
            hier_cc_data_ops.copyData(d_phi_idx,
                                      phi_current_idx,
                                      /*interior_only*/ true);
        }
        else
        {
            TBOX_ERROR("SurfaceTensionForceFunction::setDataOnPatchHierarchy : "
                       << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                          "MIDPOINT_RULE"
                       << std::endl);
        }

        // Fill ghost cells
        RobinBcCoefStrategy<NDIM>* phi_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(phi_cc_var)).front();
        interface_bc_coef = phi_bc_coef;
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent phi_transaction(
            d_phi_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, phi_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> phi_fill_op = new HierarchyGhostCellInterpolation();
        phi_fill_op->initializeOperatorState(phi_transaction, hierarchy, coarsest_ln, finest_ln);

        phi_fill_op->fillData(data_time);



        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                TBOX_ASSERT(patch->checkAllocated(d_alpha_idx));
                TBOX_ASSERT(patch->checkAllocated(d_phi_idx));
                computeAlphaFromLevelSetOnPatch(patch);
            }
        }

        // using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        // InterpolationTransactionComponent alpha_smooth_transaction(
        //     d_alpha_smooth_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, phi_bc_coef);
        // Pointer<HierarchyGhostCellInterpolation> alpha_smooth_fill_op = new HierarchyGhostCellInterpolation();
        // alpha_smooth_fill_op->initializeOperatorState(alpha_smooth_transaction, hierarchy, coarsest_ln, finest_ln);

        // alpha_smooth_fill_op->fillData(data_time);
        // mollifyData(d_alpha_smooth_idx, coarsest_ln, finest_ln, data_time, hierarchy, alpha_smooth_fill_op);
    }

    else // if d_ls_var is null, d_vof_var is not null
    {
        Pointer<CellVariable<NDIM, double>> vof_cc_var = d_vof_var;
        if (vof_cc_var.isNull())
        {
            TBOX_ERROR("VOFSurfaceTensionForceFunction: VOF variable is null. Please provide a valid VOF variable.");
        }

        const int vof_new_idx = var_db->mapVariableAndContextToIndex(vof_cc_var, d_adv_diff_solver->getNewContext());
        const int vof_current_idx =
            var_db->mapVariableAndContextToIndex(vof_cc_var, d_adv_diff_solver->getCurrentContext());

#if !defined(NDEBUG)
        TBOX_ASSERT(vof_new_idx >= 0);
        TBOX_ASSERT(vof_current_idx >= 0);
#endif

        // Copy level set into phi and C.
        if (d_ts_type == MIDPOINT_RULE)
        {
            hier_cc_data_ops.linearSum(d_alpha_idx,
                                       0.5,
                                       vof_new_idx,
                                       0.5,
                                       vof_current_idx,
                                       /*interior_only*/ true);
        }
        else if (d_ts_type == BACKWARD_EULER)
        {
            hier_cc_data_ops.copyData(d_alpha_idx,
                                      vof_new_idx,
                                      /*interior_only*/ true);
        }
        else if (d_ts_type == FORWARD_EULER)
        {
            hier_cc_data_ops.copyData(d_alpha_idx,
                                      vof_current_idx,
                                      /*interior_only*/ true);
        }
        else
        {
            TBOX_ERROR("SurfaceTensionForceFunction::setDataOnPatchHierarchy : "
                       << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                          "MIDPOINT_RULE"
                       << std::endl);
        }

        // Fill ghost cells
        RobinBcCoefStrategy<NDIM>* vof_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(vof_cc_var)).front();
        interface_bc_coef = vof_bc_coef;
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent alpha_transaction(
            d_alpha_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, vof_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> alpha_fill_op = new HierarchyGhostCellInterpolation();
        alpha_fill_op->initializeOperatorState(alpha_transaction, hierarchy, coarsest_ln, finest_ln);

        alpha_fill_op->fillData(data_time);

      

    }

    


    TBOX_ASSERT(interface_bc_coef);

    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent alpha_smooth_transaction(
        d_alpha_smooth_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, interface_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> alpha_smooth_fill_op = new HierarchyGhostCellInterpolation();
    alpha_smooth_fill_op->initializeOperatorState(alpha_smooth_transaction, hierarchy, coarsest_ln, finest_ln);

    // alpha_smooth_fill_op->fillData(data_time);

    hier_cc_data_ops.copyData(d_alpha_smooth_idx, d_alpha_idx, /*interior_only*/ false);
    mollifyData(d_alpha_smooth_idx, coarsest_ln, finest_ln, data_time, hierarchy, alpha_smooth_fill_op);

    // Mollify C.

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // Limit the surface tension force if necessary. This is mainly used in phase change problems
    // to activate the surface tension force only at the liquid-gas interface.
    if (d_mask_surface_tension_force)
    {
        const double apply_time = data_time;
        const double current_time = data_time;
        const double new_time = data_time;
        d_mask_surface_tension_force(data_idx,
                                     d_hier_math_ops,
                                     -1 /*cycle_num*/,
                                     apply_time,
                                     current_time,
                                     new_time,
                                     d_mask_surface_tension_force_ctx);
    }

    // Deallocate and remove scratch/smooth phi.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_alpha_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_alpha_smooth_idx);

        if (d_phi_idx != IBTK::invalid_index)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(d_phi_idx);
        }
    }
    // var_db->removePatchDataIndex(d_alpha_idx);
    // d_alpha_idx = IBTK::invalid_index;

    // if (d_phi_idx != IBTK::invalid_index)
    // {
    //     var_db->removePatchDataIndex(d_phi_idx);
    //     d_phi_idx = IBTK::invalid_index;
    // }
}

// double VOFSurfaceTensionForceFunction::copmuteFallbackCurvature(
//     const SAMRAI::pdat::CellData<NDIM, double>& alpha,
//     const SAMRAI::pdat::CellIndex<NDIM>& ci,
//     const double* const dx) const

// Need to think about how to implement this.

void
VOFSurfaceTensionForceFunction::setDataOnPatch(const int data_idx,
                                               Pointer<Variable<NDIM>> /*var*/,
                                               Pointer<Patch<NDIM>> patch,
                                               const double data_time,
                                               const bool initial_time,
                                               Pointer<PatchLevel<NDIM>> level)
{
    Pointer<PatchData<NDIM>> f_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_data);
#endif
    Pointer<CellData<NDIM, double>> f_cc_data = f_data;
    Pointer<SideData<NDIM, double>> f_sc_data = f_data;
#if !defined(NDEBUG)
    TBOX_ASSERT(f_cc_data || f_sc_data);
#endif
    if (f_cc_data) f_cc_data->fillAll(0.0);
    if (f_sc_data) f_sc_data->fillAll(0.0);

    if (initial_time) return;

    if (f_cc_data) setDataOnPatchCell(f_cc_data, patch, data_time, initial_time, level);
    if (f_sc_data)
    {
        // runtime polymorphism to call the child class (vof or level set) function to set the data on the patch
        // interior
        setDataOnPatchSide(f_sc_data, patch, data_time, initial_time, level);

        PatchSideDataOpsReal<NDIM, double> patch_sc_data_ops;
        if (d_compute_surface_tension_coef)
        {
            // Compute variable surface tension coefficient sigma as F = sigma*F.
            const double apply_time = data_time;
            const double current_time = data_time;
            const double new_time = data_time;
            d_compute_surface_tension_coef(data_idx,
                                           patch,
                                           -1 /*cycle_num*/,
                                           apply_time,
                                           current_time,
                                           new_time,
                                           d_compute_surface_tension_coef_ctx);
        }
        else
        {
            patch_sc_data_ops.scale(f_sc_data, d_sigma, f_sc_data, patch->getBox());
        }
    }
    return;
} // setDataOnPatch

void
VOFSurfaceTensionForceFunction::setDataOnPatchSide(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_data,
                                                   SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                   const double,
                                                   const bool,
                                                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>)
{
TBOX_ASSERT(d_alpha_idx >= 0);
TBOX_ASSERT(patch->checkAllocated(d_alpha_idx));
TBOX_ASSERT(patch->checkAllocated(d_alpha_smooth_idx));

    // conditional directive
    #if (NDIM == 2)
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();



    Pointer<CellData<NDIM, double>> alpha = patch->getPatchData(d_alpha_idx);
    Pointer<CellData<NDIM, double>> alpha_smooth = patch->getPatchData(d_alpha_smooth_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(alpha);
#endif

    CellData<NDIM, double> K_data(patch_box, 1, IntVector<NDIM>(2));  // cell centered curvature
    CellData<NDIM, double> K_valid(patch_box, 1, IntVector<NDIM>(2)); // cell centered curvature validity

    K_data.fillAll(0.0);  // initialize curvature to zero
    K_valid.fillAll(0.0); // initialize curvature to zero

    HF_CURVATURE_FC(K_data.getPointer(),
                    K_valid.getPointer(),
                    K_data.getGhostCellWidth().max(), // Add here #if (NDIM == 2)
                    alpha->getPointer(),
                    alpha->getGhostCellWidth().max(),
                    patch_box.lower(0),
                    patch_box.upper(0),
                    patch_box.lower(1),
                    patch_box.upper(1),
                    dx,
                    d_hf_column_half_width,
                    d_pure_cell_tol,
                    d_gradient_tol);



    // SC_SURFACE_TENSION_FORCE_FC(F_data->getPointer(0), // Add here #if (NDIM == 2)
    //                             F_data->getPointer(1),
    //                             F_data->getGhostCellWidth().max(),
    //                             K_data.getPointer(),
    //                             K_data.getGhostCellWidth().max(),
    //                             grad_alpha.getPointer(0, 0),
    //                             grad_alpha.getPointer(1, 1),
    //                             grad_alpha.getGhostCellWidth().max(),
    //                             patch_box.lower(0),
    //                             patch_box.upper(0),
    //                             patch_box.lower(1),
    //                             patch_box.upper(1));

    SideData<NDIM, double> grad_alpha_smooth(patch_box, NDIM, IntVector<NDIM>(2)); // side centered surface tension force

SC_NORMAL_FC(grad_alpha_smooth.getPointer(0, 0), // Add here #if (NDIM == 2)
grad_alpha_smooth.getPointer(0, 1),
grad_alpha_smooth.getPointer(1, 0),
             grad_alpha_smooth.getPointer(1, 1),
             grad_alpha_smooth.getGhostCellWidth().max(),
             alpha_smooth->getPointer(),
             alpha_smooth->getGhostCellWidth().max(),
             patch_box.lower(0),
             patch_box.upper(0),
             patch_box.lower(1),
             patch_box.upper(1),
             dx);
             

    CellData<NDIM, double> K_cv(patch_box, 1, IntVector<NDIM>(1));  // cell centered curvature

    K_cv.fillAll(0.0);

    CC_CURVATURE_FC(K_cv.getPointer(),
                    K_cv.getGhostCellWidth().max(),
                    grad_alpha_smooth.getPointer(0, 0),
                    grad_alpha_smooth.getPointer(0, 1),
#if (NDIM == 3)
grad_alpha_smooth.getPointer(0, 2),
#endif
grad_alpha_smooth.getPointer(1, 0),
grad_alpha_smooth.getPointer(1, 1),
#if (NDIM == 3)
grad_alpha_smooth.getPointer(1, 2),
grad_alpha_smooth.getPointer(2, 0),
grad_alpha_smooth.getPointer(2, 1),
grad_alpha_smooth.getPointer(2, 2),
#endif
grad_alpha_smooth.getGhostCellWidth().max(),
                    patch_box.lower(0),
                    patch_box.upper(0),
                    patch_box.lower(1),
                    patch_box.upper(1),
#if (NDIM == 3)
                    patch_box.lower(2),
                    patch_box.upper(2),
#endif
                    dx);

Box<NDIM> K_box = patch_box;
K_box.grow(IntVector<NDIM>(1));

for (Box<NDIM>::Iterator it(K_box); it; it++)
{
    const CellIndex<NDIM> ci(it());

    const double a = (*alpha)(ci);

    const bool interface_cell =
        a > d_pure_cell_tol &&
        a < 1.0 - d_pure_cell_tol;

    if (interface_cell &&
        K_valid(ci) < 0.5)
    {
        K_data(ci) = K_cv(ci);

        /*
         * From this point onward K_valid means that K_data contains
         * usable curvature, either SHF or CV.
         */
        K_valid(ci) = 1.0;
    }
}

SideData<NDIM, double> grad_alpha(patch_box, NDIM, IntVector<NDIM>(2)); // side centered surface tension force

SC_NORMAL_FC(grad_alpha.getPointer(0, 0), // Add here #if (NDIM == 2)
             grad_alpha.getPointer(0, 1),
             grad_alpha.getPointer(1, 0),
             grad_alpha.getPointer(1, 1),
             grad_alpha.getGhostCellWidth().max(),
             alpha->getPointer(),
             alpha->getGhostCellWidth().max(),
             patch_box.lower(0),
             patch_box.upper(0),
             patch_box.lower(1),
             patch_box.upper(1),
             dx);

//         SC_SURFACE_TENSION_FORCE_FC(F_data->getPointer(0),
//         F_data->getPointer(1),
// #if (NDIM == 3)
//         F_data->getPointer(2),
// #endif
//         F_data->getGhostCellWidth().max(),
//         K.getPointer(),
//         K.getGhostCellWidth().max(),
//         N.getPointer(0, 0),
//         N.getPointer(1, 1),
// #if (NDIM == 3)
//         N.getPointer(2, 2),
// #endif
//         N.getGhostCellWidth().max(),
//         patch_box.lower(0),
//         patch_box.upper(0),
//         patch_box.lower(1),
//         patch_box.upper(1)
// #if (NDIM == 3)
//             ,
//         patch_box.lower(2),
//         patch_box.upper(2)
// #endif
//     );


//     CellData<NDIM, double> K(patch_box, 1, IntVector<NDIM>(1));

//     CC_CURVATURE_FC(K.getPointer(),
//                     K.getGhostCellWidth().max(),
//                     N.getPointer(0, 0),
//                     N.getPointer(0, 1),
// #if (NDIM == 3)
//                     N.getPointer(0, 2),
// #endif
//                     N.getPointer(1, 0),
//                     N.getPointer(1, 1),
// #if (NDIM == 3)
//                     N.getPointer(1, 2),
//                     N.getPointer(2, 0),
//                     N.getPointer(2, 1),
//                     N.getPointer(2, 2),
// #endif
//                     N.getGhostCellWidth().max(),
//                     patch_box.lower(0),
//                     patch_box.upper(0),
//                     patch_box.lower(1),
//                     patch_box.upper(1),
// #if (NDIM == 3)
//                     patch_box.lower(2),
//                     patch_box.upper(2),
// #endif
//                     dx);

//     Pointer<CellData<NDIM, double>> C = patch->getPatchData(d_C_idx);
//     // N = Grad C
//     SC_NORMAL_FC(N.getPointer(0, 0),
//                  N.getPointer(0, 1),
// #if (NDIM == 3)
//                  N.getPointer(0, 2),
// #endif
//                  N.getPointer(1, 0),
//                  N.getPointer(1, 1),
// #if (NDIM == 3)
//                  N.getPointer(1, 2),
//                  N.getPointer(2, 0),
//                  N.getPointer(2, 1),
//                  N.getPointer(2, 2),
// #endif
//                  N.getGhostCellWidth().max(),
//                  C->getPointer(),
//                  C->getGhostCellWidth().max(),
//                  patch_box.lower(0),
//                  patch_box.upper(0),
//                  patch_box.lower(1),
//                  patch_box.upper(1),
// #if (NDIM == 3)
//                  patch_box.lower(2),
//                  patch_box.upper(2),
// #endif
//                  dx);

SC_SURFACE_TENSION_FORCE_VOF_FC(
    F_data->getPointer(0),
    F_data->getPointer(1),
    F_data->getGhostCellWidth().max(),
    K_data.getPointer(),
    K_valid.getPointer(),
    K_data.getGhostCellWidth().max(),
    alpha->getPointer(),
    alpha->getGhostCellWidth().max(),
    grad_alpha.getPointer(0, 0),
    grad_alpha.getPointer(1, 1),
    grad_alpha.getGhostCellWidth().max(),
    patch_box.lower(0),
    patch_box.upper(0),
    patch_box.lower(1),
    patch_box.upper(1));

    #else
    TBOX_ERROR("Only implemented for 2D");
    #endif
    return;
}



void
VOFSurfaceTensionForceFunction::mollifyData(int smooth_C_idx,
                                                 int coarsest_ln,
                                                 int finest_ln,
                                                 double data_time,
                                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                 Pointer<HierarchyGhostCellInterpolation> fill_op)
{
   

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> smooth_C_data = patch->getPatchData(smooth_C_idx);
            CellData<NDIM, double> C_data(patch_box, /*depth*/ 1, smooth_C_data->getGhostCellWidth());

            C_data.copy(*smooth_C_data);

            double* const V = smooth_C_data->getPointer(0);
            const double* const U = C_data.getPointer(0);
            const int V_gcw = (smooth_C_data->getGhostCellWidth()).max();
            const int U_gcw = (C_data.getGhostCellWidth()).max();


                MOLLIFY_IB_4_FC(V,
                                V_gcw,
                                U,
                                U_gcw,
                                patch_box.lower(0),
                                patch_box.upper(0),
                                patch_box.lower(1),
                                patch_box.upper(1)
#if (NDIM == 3)
                                    ,
                                patch_box.lower(2),
                                patch_box.upper(2)
#endif
                );
          
        }
    }

    if (fill_op) fill_op->fillData(data_time);

    return;
} // mollifyData


void
VOFSurfaceTensionForceFunction::registerSurfaceTensionForceMasking(MaskSurfaceTensionForcePtr callback, void* ctx)
{
    d_mask_surface_tension_force = callback;
    d_mask_surface_tension_force_ctx = ctx;
}

void
VOFSurfaceTensionForceFunction::setDataOnPatchCell(Pointer<CellData<NDIM, double>> /*F_data*/,
                                                   Pointer<Patch<NDIM>> /*patch*/,
                                                   const double /*data_time*/,
                                                   const bool /*initial_time*/,
                                                   Pointer<PatchLevel<NDIM>> /*level*/)
{
    TBOX_ERROR(d_object_name << "::setDataOnPatchCell(): "
                             << "cell-centered surface tension force "
                             << "is not implemented.\n");
}

void
VOFSurfaceTensionForceFunction::registerSurfaceTensionCoefficientFunction(ComputeSurfaceTensionCoefficientPtr callback,
                                                                          void* ctx)
{
    d_compute_surface_tension_coef = callback;
    d_compute_surface_tension_coef_ctx = ctx;
}

// SideData<NDIM, double> N(patch_box, NDIM, IntVector<NDIM>(2)); // side centered normal vector

// SC_NORMAL_FC(N.getPointer(0, 0),
//              N.getPointer(0, 1),
//              N.getPointer(1, 0),
//              N.getPointer(1, 1),
//              N.getGhostCellWidth().max(),
//              alpha.getPointer(),
//              alpha.getGhostCellWidth().max(),
//              patch_box.lower(0),
//              patch_box.upper(0),
//              patch_box.lower(1),
//              patch_box.upper(1),
//              dx);

// Box<NDIM> K_box = patch_box;
// K_box.grow(IntVector<NDIM>(1)); // grow by 1 for the curvature computation

// for (Box<NDIM>::Iterator it(K_box); it; it++)
// {
//     const CellIndex<NDIM> ci(it());

//     const CellIndex<NDIM> ci_plus_x = ci + CellIndex<NDIM>(1, 0);
//     const CellIndex<NDIM> ci_minus_x = ci - CellIndex<NDIM>(1, 0);
//     const CellIndex<NDIM> ci_plus_y = ci + CellIndex<NDIM>(0, 1);
//     const CellIndex<NDIM> ci_minus_y = ci - CellIndex<NDIM>(0, 1);

//     // compute curvature using height function method
//     // K(ci) = computeHeightFunctionCurvature(alpha, ci, dx);

//     const double alpha_x = (alpha(ci_plus_x) - alpha(ci_minus_x)) / (2.0 * dx[0]);
//     const double alpha_y = (alpha(ci_plus_y) - alpha(ci_minus_y)) / (2.0 * dx[1]);

//     const double grad_alpha_mag = std::sqrt(alpha_x * alpha_x + alpha_y * alpha_y);

//     if (grad_alpha_mag<= d_gradient_tol)
//     {
//         K(ci) = 0.0;
//         constinue;
//         // No interface in this neighborhood, so set curvature to zero.
//     }
//     else
//     {
//         //Need to find out if normal is primary in y or in x.
//         // If |apha_x| > |alpha_y|, then normal is primary in x direction, else normal is primary in y direction.

//         const bool use_vertical_column = (std::abs(alpha_x) > std::abs(alpha_y));
//         double h_minus = 0.0;
//         double h_center = 0.0;
//         double h_plus = 0.0;

//         bool valid_column = true;
//         // height function curvature computation
//         // K(ci) = computeHeightFunctionCurvature(alpha, ci, dx);
//     }
// }

// TBOX_ERROR("VOFSurfaceTensionForceFunction::setDataOnPatchSide() is not implemented yet.");

} // namespace IBAMR
