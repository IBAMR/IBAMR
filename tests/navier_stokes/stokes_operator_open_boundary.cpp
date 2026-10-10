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

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Apply the staggered Stokes operator to non-smooth data in two ways and compare the results.  The first is the
// library's application, which uses the reflected normal velocity ghost values, the pressure boundary value -g, and the
// added viscous term that imposes TRACTION conditions where the normal velocity is not prescribed.  The second is
// assembled here with divergence-free normal velocity ghost values and the pressure ghost values p_G = 2*p_b - p_I with
// p_b = 2*mu*du_n/dx_n - g.
void
check_open_boundary_ghost_values(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                 Pointer<INSHierarchyIntegrator> ins_integrator,
                                 const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    const double fill_time = 0.5;
    const double mu = ins_integrator->getStokesSpecifications()->getMu();

    // Configure the integrator's boundary condition objects.
    const vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(ins_integrator->getStokesSpecifications());
        U_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
        U_bc_coef->setSolutionTime(fill_time);
    }
    RobinBcCoefStrategy<NDIM>* P_bc_coef = ins_integrator->getPressureBoundaryConditions();
    auto P_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(P_bc_coef);
    P_ins_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
    P_ins_bc_coef->setSolutionTime(fill_time);

    // Data: u and p are the operator's argument, u2 is a copy of u with divergence-free ghost values, and y1 and y2 are
    // the two results.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("open_boundary_ghost_values");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    auto register_data = [&](Pointer<Variable<NDIM>> var, const int ghosts)
    {
        const int idx = var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(ghosts));
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(idx, fill_time);
        }
        return idx;
    };
    const int u_idx = register_data(u_var, 1), p_idx = register_data(p_var, 1);
    const int u2_idx = register_data(new SideVariable<NDIM, double>("u2"), 1);
    Pointer<SideVariable<NDIM, double>> y1u_var = new SideVariable<NDIM, double>("y1u");
    Pointer<CellVariable<NDIM, double>> y1p_var = new CellVariable<NDIM, double>("y1p");
    const int y1u_idx = register_data(y1u_var, 0), y1p_idx = register_data(y1p_var, 0);
    const int y2u_idx = register_data(new SideVariable<NDIM, double>("y2u"), 0);
    const int y2p_idx = register_data(new CellVariable<NDIM, double>("y2p"), 0);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            u_data->fillAll(0.0);
            p_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                {
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                        std::sin(1.7 * it()(0) + 2.9 * it()(1) + 0.8 * axis);
                }
            }
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                (*p_data)(it()) = std::cos(2.3 * it()(0) + 0.9 * it()(1));
            }
        }
    }

    // Weights, vectors, and the operator.
    HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
    SAMRAIVectorReal<NDIM, double> x("x", patch_hierarchy, 0, finest_ln), y1("y1", patch_hierarchy, 0, finest_ln);
    x.addComponent(u_var, u_idx, wgt_sc_idx);
    x.addComponent(p_var, p_idx, wgt_cc_idx);
    y1.addComponent(y1u_var, y1u_idx, wgt_sc_idx);
    y1.addComponent(y1p_var, y1p_idx, wgt_cc_idx);
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, fill_time, patch_hierarchy);
    PoissonSpecifications spec("spec");
    spec.setCConstant(1.0);
    spec.setDConstant(-mu);
    StaggeredStokesOperator op("StaggeredStokesOperator", /*homogeneous_bc*/ false);
    op.setVelocityPoissonSpecifications(spec);
    op.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    op.setPhysicalBoundaryHelper(bc_helper);
    op.setSolutionTime(fill_time);
    op.setTimeInterval(fill_time, fill_time);
    op.initializeOperatorState(x, y1);

    // (i) The library's application; it fills the ghost values of x.
    op.apply(x, y1);

    // (ii) Divergence-free normal velocity ghost values and the corresponding pressure ghost values (which replace
    // those of x) where the normal velocity is not prescribed, and the operator assembled from its terms.
    HierarchySideDataOpsReal<NDIM, double> sc_ops(patch_hierarchy, 0, finest_ln);
    HierarchyCellDataOpsReal<NDIM, double> cc_ops(patch_hierarchy, 0, finest_ln);
    sc_ops.copyData(u2_idx, u_idx, /*interior_only*/ false);
    bc_helper->enforceDivergenceFreeConditionAtBoundary(u2_idx);
    double max_ghost_difference = 0.0;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> u2_data = patch->getPatchData(u2_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            const tbox::Array<BoundaryBox<NDIM>>& bdry_boxes = pgeom->getCodimensionBoundaries(1);
            for (int k = 0; k < bdry_boxes.size(); ++k)
            {
                const unsigned int axis = bdry_boxes[k].getLocationIndex() / 2;
                const bool is_lower = bdry_boxes[k].getLocationIndex() % 2 == 0;
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_boxes[k]);
                Pointer<ArrayData<NDIM, double>> a_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> b_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> g_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                u_bc_coefs[axis]->setBcCoefs(
                    a_data, b_data, g_data, Pointer<Variable<NDIM>>(), *patch, bdry_boxes[k], fill_time);
                // Limit the boundary cells to the tangential extent of the patch.
                Box<NDIM> cell_box = patch->getBox();
                cell_box.lower(axis) -= 1;
                cell_box.upper(axis) += 1;
                for (Box<NDIM>::Iterator it(bc_coef_box * cell_box); it; it++)
                {
                    // This is the interior cell abutting a lower boundary and the ghost cell abutting an upper
                    // boundary, so its lower face is the boundary face.
                    if (!IBTK::rel_equal_eps((*b_data)(it(), 0), 1.0))
                    {
                        continue;
                    }
                    SideIndex<NDIM> s_in(it(), axis, SideIndex<NDIM>::Lower), s_out(s_in);
                    s_in(axis) += is_lower ? 1 : -1;
                    s_out(axis) += is_lower ? -1 : 1;
                    max_ghost_difference =
                        std::max(max_ghost_difference, std::abs((*u2_data)(s_out) - (*u_data)(s_out)));
                    const double du_dx =
                        (is_lower ? 1.0 : -1.0) * ((*u2_data)(s_in) - (*u2_data)(s_out)) / (2.0 * pgeom->getDx()[axis]);
                    hier::Index<NDIM> i_in(it()), i_ghost(it());
                    (is_lower ? i_ghost : i_in)(axis) -= 1;
                    (*p_data)(i_ghost) = 2.0 * (2.0 * mu * du_dx - (*g_data)(it(), 0)) - (*p_data)(i_in);
                }
            }
        }
    }
    hier_math_ops.grad(y2u_idx, u_var, /*cf_bdry_synch*/ false, 1.0, p_idx, p_var, nullptr, fill_time);
    hier_math_ops.laplace(y2u_idx, u_var, spec, u2_idx, u_var, nullptr, fill_time, 1.0, y2u_idx, u_var);
    hier_math_ops.div(y2p_idx, p_var, -1.0, u2_idx, u_var, nullptr, fill_time, /*cf_bdry_synch*/ true);
    bc_helper->copyDataAtDirichletBoundaries(y2u_idx, u2_idx);

    // Compare.
    const double velocity_norm = sc_ops.maxNorm(y1u_idx, wgt_sc_idx);
    const double result_norm = std::max(velocity_norm, cc_ops.maxNorm(y1p_idx, wgt_cc_idx));
    sc_ops.subtract(y2u_idx, y1u_idx, y2u_idx);
    cc_ops.subtract(y2p_idx, y1p_idx, y2p_idx);
    const double difference = std::max(sc_ops.maxNorm(y2u_idx, wgt_sc_idx), cc_ops.maxNorm(y2p_idx, wgt_cc_idx));
    if (max_ghost_difference == 0.0 || !std::isfinite(max_ghost_difference) || !std::isfinite(velocity_norm) ||
        !std::isfinite(difference))
    {
        TBOX_ERROR("the normal velocity ghost values do not differ or a result is not finite\n");
    }
    if (difference > 1.0e-10 * result_norm)
    {
        TBOX_ERROR("the two results differ by " << difference << " and the largest magnitude of the result is "
                                                << result_norm << "\n");
    }
    plog << std::setprecision(10)
         << "Largest difference between the reflected and divergence-free normal velocity ghost values: "
         << max_ghost_difference << "\nLargest magnitude of the velocity component of the result: " << velocity_norm
         << "\n";
} // check_open_boundary_ghost_values
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                         app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                                         grid_geometry);
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        check_open_boundary_ghost_values(patch_hierarchy, ins_integrator, u_bc_coefs);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
