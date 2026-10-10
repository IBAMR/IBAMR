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
#include <iomanip>
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
    // Fill the ghost values of a velocity field that satisfies the boundary conditions and report the
    // error in the tangential ghost values outside the x boundaries: at the lower corners, at the upper corners,
    // elsewhere within the extent of the patch along the boundary, and at the outermost ghost faces beyond that
    // extent.
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
    double outermost_error = 0.0;
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
            const auto exact = [&](const unsigned int axis, const hier::Index<NDIM>& i)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (i(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                }
                return parsers[axis].Eval();
            };

            // Set the values at the faces of cells in the physical domain,
            // including in the ghost cells, as a ghost cell fill from the other
            // patches would.
            u_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    hier::Index<NDIM> i_lower = b();
                    i_lower(axis) -= 1;
                    if (domain.contains(b()) || domain.contains(i_lower))
                    {
                        u_data->getArrayData(axis)(b(), 0) = exact(axis, b());
                    }
                }
            }
            bdry_op.setPatchDataIndex(u_idx);
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, ghost_width);

            // Tangential ghost values outside the x boundaries, within the
            // extent of the patch in the other directions, and the outermost
            // faces of the patch data along the component's own axis.
            for (unsigned int axis = 1; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis);
                for (Box<NDIM>::Iterator b(ghost_side_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    bool beyond_x_boundary = i(0) < patch_box.lower(0) || i(0) > patch_box.upper(0);
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
                    // The ghost value is set using the normal velocity at the
                    // boundary faces in the rows i(axis) - 1 and i(axis). It is
                    // at a lower or upper corner of the boundary if the lower or
                    // upper row is outside the physical domain.
                    hier::Index<NDIM> i_boundary = i;
                    i_boundary(0) = i(0) < patch_box.lower(0) ? patch_box.lower(0) : patch_box.upper(0);
                    hier::Index<NDIM> i_boundary_lower = i_boundary;
                    i_boundary_lower(axis) -= 1;
                    const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact(axis, i));
                    const bool lower_row_in_domain = domain.contains(i_boundary_lower);
                    const bool upper_row_in_domain = domain.contains(i_boundary);
                    if (outermost)
                    {
                        if (lower_row_in_domain && upper_row_in_domain)
                        {
                            outermost_error = std::max(outermost_error, error);
                            ++outermost_count;
                        }
                    }
                    else if (!lower_row_in_domain)
                    {
                        lower_end_error = std::max(lower_end_error, error);
                        ++lower_end_count;
                    }
                    else if (!upper_row_in_domain)
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
    // The errors that the input expects to vanish are at the level of round-off, which the expected output must not
    // depend on, so they are checked against the bound in the input and the output reports the number of faces
    // checked along with the errors.
    const double error_bound = exact_velocity_db->getDouble("error_bound");
    const double max_error = IBTK_MPI::maxReduction(
        std::max(std::max(lower_end_error, upper_end_error), std::max(other_error, outermost_error)));
    if (!(max_error <= error_bound))
    {
        TBOX_ERROR("traction_corner_accuracy: the largest tangential ghost error " << max_error << " exceeds the bound "
                                                                                   << error_bound << "\n");
    }
    pout << std::setprecision(6) << std::scientific
         << "ghost faces checked at the lower corners of the x boundaries = " << IBTK_MPI::sumReduction(lower_end_count)
         << '\n'
         << "ghost faces checked at the upper corners of the x boundaries = " << IBTK_MPI::sumReduction(upper_end_count)
         << '\n'
         << "ghost faces checked elsewhere on the x boundaries            = " << IBTK_MPI::sumReduction(other_count)
         << '\n'
         << "ghost faces checked at the outermost ghost faces             = " << IBTK_MPI::sumReduction(outermost_count)
         << '\n'
         << "max tangential ghost error at the lower corners of the x boundaries = "
         << IBTK_MPI::maxReduction(lower_end_error) << '\n'
         << "max tangential ghost error at the upper corners of the x boundaries = "
         << IBTK_MPI::maxReduction(upper_end_error) << '\n'
         << "max tangential ghost error elsewhere on the x boundaries            = "
         << IBTK_MPI::maxReduction(other_error) << '\n'
         << "max tangential ghost error at the outermost ghost faces             = "
         << IBTK_MPI::maxReduction(outermost_error) << '\n';
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

        // The velocity boundary condition objects of the integrator are set up at the start of a time step.
        const double current_time = ins_integrator->getIntegratorTime();
        ins_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);
        check_traction_corner_accuracy(
            patch_hierarchy, ins_integrator, app_initializer->getComponentDatabase("ExactVelocity"));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
