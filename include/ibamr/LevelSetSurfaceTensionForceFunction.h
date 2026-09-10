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

#ifndef included_IBAMR_LevelSetSurfaceTensionForceFunction
#define included_IBAMR_LevelSetSurfaceTensionForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/SurfaceTensionForceFunction.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 */
class LevelSetSurfaceTensionForceFunction : public SurfaceTensionForceFunction
{
public:
    LevelSetSurfaceTensionForceFunction(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                        const AdvDiffHierarchyIntegrator* adv_diff_solver,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> level_set_var);

    ~LevelSetSurfaceTensionForceFunction() override = default;

    /*
     * Level-set-specific configuration.
     */
    void setSmoother(const std::string& kernel_fcn)
    {
        d_kernel_fcn = kernel_fcn;
    }

    std::string getSmoother() const
    {
        return d_kernel_fcn;
    }

    void setNumberOfInterfaceCells(double m)
    {
        d_num_interface_cells = m;
    }

    double getNumberOfInterfaceCells() const
    {
        return d_num_interface_cells;
    }

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

    TimeSteppingType d_ts_type = MIDPOINT_RULE;

    int d_phi_idx = IBTK::invalid_index;
    int d_C_idx = IBTK::invalid_index;

    std::string d_kernel_fcn = "none";

    double d_sigma = 1.0;
    double d_num_interface_cells = 2.0;

    int getStencilSize(const std::string& kernel_fcn) const;

    int getMinimumGhostWidth(const std::string& kernel_fcn) const;

private:
    void setDataOnPatchSide(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                            double data_time,
                            bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    void setDataOnPatchCell(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> F_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                            double data_time,
                            bool initial_time,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    void convertToHeaviside(int C_idx,
                            int coarsest_ln,
                            int finest_ln,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    void mollifyData(int C_idx,
                     int coarsest_ln,
                     int finest_ln,
                     double data_time,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                     SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> fill_op);

    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    MaskSurfaceTensionForcePtr d_mask_surface_tension_force = nullptr;
    void* d_mask_surface_tension_force_ctx = nullptr;

    ComputeSurfaceTensionCoefficientPtr d_compute_surface_tension_coef = nullptr;
    void* d_compute_surface_tension_coef_ctx = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_SurfaceTensionForceFunction
