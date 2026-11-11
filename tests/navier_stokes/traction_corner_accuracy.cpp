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

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <muParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
void
check_traction_corner_accuracy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               Pointer<INSHierarchyIntegrator> ins_integrator,
                               Pointer<Database> exact_velocity_db)
{
    // Fill the ghost values of a velocity field that satisfies the boundary conditions and check the error in the
    // tangential ghost values outside the x boundaries, at corners and elsewhere.
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
        var_db->registerVariableAndContext(u_var, var_db->getContext("traction_corner_accuracy"), ghost_width);
    CartSideRobinPhysBdryOp bdry_op(u_idx, ins_integrator->getVelocityBoundaryConditions(), /*homogeneous_bc*/ false);
    const double fill_time = ins_integrator->getIntegratorTime();

    double lower_end_error = 0.0;
    double upper_end_error = 0.0;
    double other_error = 0.0;
    int lower_end_count = 0;
    int upper_end_count = 0;
    int other_count = 0;
    int outermost_count = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
            domain.refine(pgeom->getRatio());
            const Box<NDIM> domain_box = domain.getBoundingBox();
            const auto exact_velocity = [&](const unsigned int axis, const hier::Index<NDIM>& i)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (i(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                }
                return parsers[axis].Eval();
            };

            // Set the values in the physical domain, including ghost cells, as a ghost cell fill would.
            u_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    hier::Index<NDIM> i_lower = b();
                    i_lower(axis) -= 1;
                    if (domain.contains(b()) || domain.contains(i_lower))
                    {
                        u_data->getArrayData(axis)(b(), 0) = exact_velocity(axis, b());
                    }
                }
            }
            bdry_op.setPatchDataIndex(u_idx);
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, ghost_width);

            // Tangential ghost values outside the x boundaries, including the outermost faces of the patch data.
            for (unsigned int axis = 1; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis);
                for (Box<NDIM>::Iterator b(ghost_side_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    bool beyond_x_boundary = i(0) < domain_box.lower(0) || i(0) > domain_box.upper(0);
                    bool outermost = false;
                    for (unsigned int d = 1; d < NDIM; ++d)
                    {
                        const bool within_patch = side_box.lower(d) <= i(d) && i(d) <= side_box.upper(d);
                        if (d == axis && !within_patch)
                        {
                            outermost = i(d) == ghost_side_box.lower(d) || i(d) == ghost_side_box.upper(d);
                            beyond_x_boundary = beyond_x_boundary && outermost;
                        }
                        else
                        {
                            beyond_x_boundary = beyond_x_boundary && within_patch;
                        }
                    }
                    if (!beyond_x_boundary)
                    {
                        continue;
                    }
                    // The ghost value is at a corner if the layer i(axis) - 1 or i(axis) of the normal velocity lies
                    // outside the physical domain.
                    hier::Index<NDIM> i_boundary = i;
                    i_boundary(0) = i(0) < patch_box.lower(0) ? patch_box.lower(0) : patch_box.upper(0);
                    hier::Index<NDIM> i_boundary_lower = i_boundary;
                    i_boundary_lower(axis) -= 1;
                    const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact_velocity(axis, i));
                    const bool lower_layer_in_domain = domain.contains(i_boundary_lower);
                    const bool upper_layer_in_domain = domain.contains(i_boundary);
                    if (outermost)
                    {
                        // The fill sets the outermost ghost faces to NaN, because their difference involves a normal
                        // velocity face outside the patch data.
                        if (lower_layer_in_domain && upper_layer_in_domain)
                        {
                            if (!std::isnan(u_data->getArrayData(axis)(i, 0)))
                            {
                                TBOX_ERROR("check_traction_corner_accuracy():\n"
                                           << "  the outermost ghost face " << i << " of axis " << axis
                                           << " is not NaN\n");
                            }
                            ++outermost_count;
                        }
                    }
                    else if (!lower_layer_in_domain)
                    {
                        lower_end_error = std::max(lower_end_error, error);
                        ++lower_end_count;
                    }
                    else if (!upper_layer_in_domain)
                    {
                        upper_end_error = std::max(upper_end_error, error);
                        ++upper_end_count;
                    }
                    else
                    {
                        other_error = std::max(other_error, error);
                        ++other_count;
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
    }
    // Check the largest error against the bound in the input and print the number of faces checked and the largest
    // error in each region.
    lower_end_error = IBTK_MPI::maxReduction(lower_end_error);
    upper_end_error = IBTK_MPI::maxReduction(upper_end_error);
    other_error = IBTK_MPI::maxReduction(other_error);
    const double error_bound = exact_velocity_db->getDouble("error_bound");
    const double max_error = std::max(std::max(lower_end_error, upper_end_error), other_error);
    if (!(max_error <= error_bound))
    {
        TBOX_ERROR("check_traction_corner_accuracy():\n"
                   << "  the largest tangential ghost error " << max_error << " exceeds the bound " << error_bound
                   << "\n");
    }
    const auto print = [](const std::string& region, const int count, const double error)
    {
        plog << "ghost faces " << region << ": " << IBTK_MPI::sumReduction(count) << " checked, largest error " << error
             << '\n';
    };
    print("at the lower corners of the x boundaries", lower_end_count, lower_end_error);
    print("at the upper corners of the x boundaries", upper_end_count, upper_end_error);
    print("elsewhere on the x boundaries", other_count, other_error);
    plog << "outermost ghost faces at the ends of the ghost box: " << IBTK_MPI::sumReduction(outermost_count)
         << " are NaN, as the fill does not compute them\n";
    return;
} // check_traction_corner_accuracy
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

        std::vector<std::unique_ptr<muParserRobinBcCoefs>> u_bc_coef_owners;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coef_owners.push_back(std::make_unique<muParserRobinBcCoefs>(
                "u_bc_coefs_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry));
            u_bc_coefs.push_back(u_bc_coef_owners.back().get());
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // The integrator sets up its velocity boundary condition objects at the start of a time step.
        const double current_time = ins_integrator->getIntegratorTime();
        ins_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);
        check_traction_corner_accuracy(
            patch_hierarchy, ins_integrator, app_initializer->getComponentDatabase("ExactVelocity"));
    }
} // main
