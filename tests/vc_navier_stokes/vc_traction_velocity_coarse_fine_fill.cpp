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
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <cmath>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
// Fill the ghost values of a smooth finite velocity field with the coarse-fine interpolation and the velocity boundary
// conditions of the integrator, as the integrator fills its velocity, and count the NaN ghost values inside the domain.
// The interpolation applies the boundary conditions to temporary coarse patches, which have no viscosity data.
void
fill_with_coarse_fine_interpolation(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                    Pointer<INSVCStaggeredHierarchyIntegrator> ins_integrator,
                                    const int u_idx,
                                    const bool homogeneous_bc)
{
    const double fill_time = ins_integrator->getIntegratorTime();
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<SideData<NDIM, double>> u_data = level->getPatch(p())->getPatchData(u_idx);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    double phase = axis;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        phase += (0.3 + 0.2 * d) * b()(d);
                    }
                    u_data->getArrayData(axis)(b(), 0) = std::sin(phase);
                }
            }
        }
    }

    // The coarse-fine interpolation fills at most as many layers as the refinement ratio.
    Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
    const IntVector<NDIM> ratio = finest_level->getRatioToCoarserLevel();
    HierarchyGhostCellInterpolation fill_op;
    fill_op.initializeOperatorState(HierarchyGhostCellInterpolation::InterpolationTransactionComponent(
                                        u_idx,
                                        "NONE",
                                        /*use_cf_bdry_interpolation*/ true,
                                        "CUBIC_COARSEN",
                                        "LINEAR",
                                        /*consistent_type_2_bdry*/ false,
                                        ins_integrator->getVelocityBoundaryConditions(),
                                        new SideNoCornersFillPattern(ratio.min(), /*overwrite_interior*/ false)),
                                    patch_hierarchy);
    fill_op.setHomogeneousBc(homogeneous_bc);
    fill_op.fillData(fill_time);

    int nan_count = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        Box<NDIM> domain_box = patch_hierarchy->getGridGeometry()->getPhysicalDomain()[0];
        domain_box.refine(level->getRatio());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> patch_side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                const Box<NDIM> domain_side_box = SideGeometry<NDIM>::toSideBox(domain_box, axis);
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    if (!patch_side_box.contains(b()) && domain_side_box.contains(b()) &&
                        std::isnan(u_data->getArrayData(axis)(b(), 0)))
                    {
                        ++nan_count;
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
    }
    pout << "NaN ghost values inside the domain after the coarse-fine fill = " << nan_count << '\n';
} // fill_with_coarse_fine_interpolation
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // A case that ends with the expected error records its message in the output.
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);
    // SAMRAI warnings contain source file paths, which differ between builds.
    Logger::getInstance()->setWarning(false);

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

        // The viscosity is a field, which requires a variable, initial values, and boundary conditions.
        Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<CartGridFunction> mu_fcn = new muParserCartGridFunction(
            "mu_fcn", app_initializer->getComponentDatabase("ViscosityFunction"), grid_geometry);
        ins_integrator->registerViscosityVariable(mu_var);
        ins_integrator->registerViscosityInitialConditions(mu_fcn);
        muParserRobinBcCoefs mu_bc_coef(
            "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBoundaryConditions"), grid_geometry);
        ins_integrator->registerViscosityBoundaryConditions(&mu_bc_coef);

        // The integrator requires the density or the viscosity to be a field.
        Pointer<CellVariable<NDIM, double>> rho_var = new CellVariable<NDIM, double>("rho");
        Pointer<CartGridFunction> rho_fcn = new muParserCartGridFunction(
            "rho_fcn", app_initializer->getComponentDatabase("DensityFunction"), grid_geometry);
        ins_integrator->registerMassDensityVariable(rho_var);
        ins_integrator->registerMassDensityInitialConditions(rho_fcn);

        // The ghost cells of the velocity are wider than those of the viscosity. The boundary boxes are sized for the
        // largest registered ghost width.
        const int ghost_width = app_initializer->getInputDatabase()->getInteger("ghost_width");
        ins_integrator->initializeHierarchyIntegrator(patch_hierarchy, gridding_algorithm);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->registerVariableAndContext(ins_integrator->getVelocityVariable(),
                                                             var_db->getContext("traction_velocity_bc"),
                                                             IntVector<NDIM>(ghost_width));
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
        fill_with_coarse_fine_interpolation(
            patch_hierarchy,
            ins_integrator,
            u_idx,
            app_initializer->getInputDatabase()->getBool("coarse_fine_fill_homogeneous_bc"));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
