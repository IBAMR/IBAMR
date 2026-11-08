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

#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredVelocityBcCoef.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <muParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Check the velocity boundary conditions of the variable-coefficient staggered solver at TRACTION boundaries: report
// the largest error in the tangential ghost values of a smooth velocity field outside the x boundaries.
void
check_traction_velocity_bc(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                           Pointer<INSHierarchyIntegrator> ins_integrator,
                           Pointer<Database> exact_velocity_db)
{
    std::array<double, NDIM> X;
    std::vector<mu::Parser> parsers(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        parsers[d].SetExpr(exact_velocity_db->getString("function_" + std::to_string(d)));
        for (unsigned int k = 0; k < NDIM; ++k)
        {
            parsers[d].DefineVar("X_" + std::to_string(k), &X[k]);
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(3);
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("traction_velocity_bc"), ghost_width);
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    CartSideRobinPhysBdryOp bdry_op(u_idx, U_bc_coefs, /*homogeneous_bc*/ false);
    const double fill_time = ins_integrator->getIntegratorTime();

    double ghost_error = 0.0;
    int ghost_count = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        if (level->getNumberOfPatches() != 1)
        {
            TBOX_ERROR("this test requires one patch per level\n"
                       << "  set largest_patch_size in the input so that each level has one patch\n");
        }
        level->allocatePatchData(u_idx, fill_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            const auto exact_velocity = [&](const unsigned int axis, const hier::Index<NDIM>& i)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (i(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                }
                return parsers[axis].Eval();
            };
            // Call fcn(axis, i, is_interior) for each side index i of the patch data.
            const auto for_each_side = [&](const auto& fcn)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                    {
                        fcn(axis, b(), side_box.contains(b()));
                    }
                }
            };

            for_each_side([&](const unsigned int axis, const hier::Index<NDIM>& i, const bool is_interior)
                          { u_data->getArrayData(axis)(i, 0) = is_interior ? exact_velocity(axis, i) : 0.0; });
            bdry_op.setPatchDataIndex(u_idx);
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, ghost_width);
            for (unsigned int axis = 1; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    bool beyond_x_boundary_only = i(0) < side_box.lower(0) || i(0) > side_box.upper(0);
                    for (unsigned int d = 1; d < NDIM; ++d)
                    {
                        beyond_x_boundary_only =
                            beyond_x_boundary_only && side_box.lower(d) <= i(d) && i(d) <= side_box.upper(d);
                    }
                    if (beyond_x_boundary_only)
                    {
                        const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact_velocity(axis, i));
                        if (std::isnan(error))
                        {
                            TBOX_ERROR("the tangential ghost value at " << i << " of axis " << axis << " is NaN\n");
                        }
                        ghost_error = std::max(ghost_error, error);
                        ++ghost_count;
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
    }

    if (ghost_count == 0 || !(ghost_error <= exact_velocity_db->getDouble("error_bound")))
    {
        TBOX_ERROR("the largest error in "
                   << ghost_count << " tangential ghost values is " << ghost_error << "\n"
                   << "  the test requires at least one value and an error of at most error_bound\n");
    }
    pout << "tangential ghost values checked outside the x boundaries = " << ghost_count << '\n'
         << std::setprecision(6) << std::scientific << "largest error in those ghost values = " << ghost_error << '\n';
    return;
} // check_traction_velocity_bc
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");

        Pointer<INSVCStaggeredHierarchyIntegrator> ins_integrator =
            new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
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

        // The integrator requires the density or the viscosity to be a field.
        Pointer<CellVariable<NDIM, double>> rho_var = new CellVariable<NDIM, double>("rho");
        Pointer<CartGridFunction> rho_fcn = new muParserCartGridFunction(
            "rho_fcn", app_initializer->getComponentDatabase("DensityFunction"), grid_geometry);
        ins_integrator->registerMassDensityVariable(rho_var);
        ins_integrator->registerMassDensityInitialConditions(rho_fcn);

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Configure the integrator's velocity boundary condition objects, as it does at the start of a time step.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            auto U_bc_coef =
                dynamic_cast<INSVCStaggeredVelocityBcCoef*>(ins_integrator->getVelocityBoundaryConditions()[d]);
            U_bc_coef->setStokesSpecifications(ins_integrator->getStokesSpecifications());
            U_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
            U_bc_coef->setSolutionTime(ins_integrator->getIntegratorTime());
        }
        check_traction_velocity_bc(
            patch_hierarchy, ins_integrator, app_initializer->getComponentDatabase("ExactVelocity"));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
