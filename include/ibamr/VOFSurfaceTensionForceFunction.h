// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_VOFSurfaceTensionForceFunction
#define included_IBAMR_VOFSurfaceTensionForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>
#include <CellVariable.h>

#include <CartesianGridGeometry.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <string>

namespace IBAMR
{
class VOFSurfaceTensionForceFunction : public SurfaceTensionForceFunction
{
public:
    VOFSurfaceTensionForceFunction(const std::string& object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                   const AdvDiffHierarchyIntegrator* adv_diff_solver,
                                   const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> level_set_var,
                                   const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var =
                                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>>(nullptr));

    ~VOFSurfaceTensionForceFunction() override = default;

    /*
     * Common surface-tension interface.
     */
    void setSurfaceTensionCoef(double sigma) override
    {
        d_sigma = sigma;
    }

    double getSurfaceTensionCoef() const override
    {
        return d_sigma;
    }

    bool isTimeDependent() const override;

    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = IBTK::invalid_level_number,
                                 int finest_ln = IBTK::invalid_level_number) override;

    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override;

    void registerSurfaceTensionForceMasking(MaskSurfaceTensionForcePtr callback, void* ctx) override;

    void registerSurfaceTensionCoefficientFunction(ComputeSurfaceTensionCoefficientPtr callback, void* ctx) override;

private:
    int getPhiGhostWidth() const;
    int getAlphaGhostWidth() const;

    void computeAlphaFromLevelSetOnPatch(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch);

    void mollifyData(int C_idx,
        int coarsest_ln,
        int finest_ln,
        double data_time,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> fill_op);

    void setDataOnPatchSide(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                            const double data_time,
                            const bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    void setDataOnPatchCell(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                            double data_time,
                            bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);


    TimeSteppingType d_ts_type = MIDPOINT_RULE;

    // // maybe shift the above two functions to protected?
    VOFSurfaceTensionForceFunction() = delete;
    VOFSurfaceTensionForceFunction(const VOFSurfaceTensionForceFunction&) = delete;

    VOFSurfaceTensionForceFunction& operator=(const VOFSurfaceTensionForceFunction&) = delete;

    // static double computeVOFFromLevelSet(double phi, const std::array<double, NDIM>& dphi); // This to be replaced by
    // VOFFromLevelSetInitializer class in vc_ins_vof_utilities.h

    // double computeFallbackCurvature (const SAMRAI::pdat::CellData<NDIM, double>& alpha,
    //                                  const SAMRAI::pdat::CellIndex<NDIM>& ci,
    //                                  const double* dx) const;

    int d_phi_idx = IBTK::invalid_index;
    int d_alpha_idx = IBTK::invalid_index;
    int d_alpha_smooth_idx = IBTK::invalid_index;

    
    double d_sigma = 1.0;

    int d_hf_column_half_width = 3;
    double d_pure_cell_tol = 1.0e-6;
    double d_gradient_tol = 1.0e-12;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_alpha_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_alpha_smooth_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_vof_var;

    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    MaskSurfaceTensionForcePtr d_mask_surface_tension_force = nullptr;

    void* d_mask_surface_tension_force_ctx = nullptr;

    ComputeSurfaceTensionCoefficientPtr d_compute_surface_tension_coef = nullptr;

    void* d_compute_surface_tension_coef_ctx = nullptr;
};

} // namespace IBAMR
#endif
