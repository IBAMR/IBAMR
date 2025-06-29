// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <tbox/Array.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure
#include <ibamr/SOAcousticStreamingBcCoefs.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

// A test program to check that the second-order acoustic Stokes drift boundary condition is working correctly.

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "so_acoustic_bcs.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create side-centered velocity variables
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<SideVariable<NDIM, double> > U1_real_var = new SideVariable<NDIM, double>("U1_real", /*depth*/ 1);
        Pointer<SideVariable<NDIM, double> > U1_imag_var = new SideVariable<NDIM, double>("U1_imag", /*depth*/ 1);
        const int gcw = 2;
        const int U1_real_idx = var_db->registerVariableAndContext(U1_real_var, context, gcw);
        const int U1_imag_idx = var_db->registerVariableAndContext(U1_imag_var, context, gcw);

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Configure boundary condition objects for the second-order system.
        std::array<IBAMR::SOAcousticStreamingBcCoefs, NDIM> so_bc_coefs;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            auto& bc_coef = so_bc_coefs[d];
            bc_coef.setSOVelocityComponent(d);
            bc_coef.setFOVelocityPressurePatchDataIndices(
                U1_real_idx, U1_imag_idx, IBTK::invalid_index, IBTK::invalid_index);
            bc_coef.setDensityPatchDataIndex(IBTK::invalid_index);
            bc_coef.setSoundSpeed(1.0);
            bc_coef.setAcousticAngularFrequency(1.0);
            bc_coef.useStokesDriftVelocityForm(true);
        }

        // Loop over patch box and fill the data
        // U1Real = "X0^3 + X1^3"
        // V1Real = "X1^2 + X0^2"
        // U1Imag = "(X0^2 + X1^2)"
        // V1Imag = "(X1^3 + X0^3)"
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U1_real_idx);
            level->allocatePatchData(U1_imag_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > U1_real_data = patch->getPatchData(U1_real_idx);
                Pointer<SideData<NDIM, double> > U1_imag_data = patch->getPatchData(U1_imag_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const patch_x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();

                std::array<double, NDIM> posn;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> ic(patch_box, axis); ic; ic++)
                    {
                        const SideIndex<NDIM>& i = ic();
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d == axis)
                            {
                                posn[d] = patch_x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                            }
                            else
                            {
                                posn[d] = patch_x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                            }
                        }
                        const double& X0 = posn[0];
                        const double& X1 = posn[1];

                        if (axis == 0)
                        {
                            (*U1_real_data)(i) = std::pow(X0, 3) + std::pow(X1, 3);
                            (*U1_imag_data)(i) = std::pow(X0, 2) + std::pow(X1, 2);
                        }

                        if (axis == 1)
                        {
                            (*U1_real_data)(i) = std::pow(X0, 2) + std::pow(X1, 2);
                            (*U1_imag_data)(i) = std::pow(X0, 3) + std::pow(X1, 3);
                        }
                    }
                }
            }
        }

        std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2> U1_bc_coefs;
        U1_bc_coefs[0].resize(NDIM, nullptr);
        U1_bc_coefs[1].resize(NDIM, nullptr);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "U1_real_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "FOVelocityRealBcCoefs_" + std::to_string(d);

            U1_bc_coefs[0][d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "U1_imag_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "FOVelocityImagBcCoefs_" + std::to_string(d);

            U1_bc_coefs[1][d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }

        // Fill the ghost cells of U1_real and U1_imag
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> comp_transactions(2);
        comp_transactions[0] = InterpolationTransactionComponent(U1_real_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 true,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "QUADRATIC",
                                                                 false,
                                                                 U1_bc_coefs[0],
                                                                 Pointer<VariableFillPattern<NDIM> >(nullptr));
        comp_transactions[1] = InterpolationTransactionComponent(U1_imag_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 true,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "QUADRATIC",
                                                                 false,
                                                                 U1_bc_coefs[1],
                                                                 Pointer<VariableFillPattern<NDIM> >(nullptr));
        HierarchyGhostCellInterpolation comp_fill_op;
        comp_fill_op.initializeOperatorState(comp_transactions, patch_hierarchy);
        comp_fill_op.fillData(/*time*/ 0.0);

        // Loop over ghost box and check that there are no NaNs/unintialized large values
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > U1_real_data = patch->getPatchData(U1_real_idx);
                Pointer<SideData<NDIM, double> > U1_imag_data = patch->getPatchData(U1_imag_idx);
                const Box<NDIM>& ghost_box = U1_real_data->getGhostBox();

                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> ic(ghost_box, axis); ic; ic++)
                    {
                        const SideIndex<NDIM>& i = ic();
                        TBOX_ASSERT(std::abs((*U1_real_data)(i)) < 1e300);
                        TBOX_ASSERT(std::abs((*U1_imag_data)(i)) < 1e300);
                    }
                }
            }
        }

        double weighted_sum = 0.0;
        double normal_area = 0.0;
        double tangential_area = 0.0;
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const tbox::Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                if (n_physical_codim1_boxes == 0) continue;

                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const patch_x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();

                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const unsigned int location_index = bdry_box.getLocationIndex();
                    const unsigned int bdry_normal_axis = location_index / 2;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, /*gcw_to_fill*/ 1);
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                    Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                    auto& bc_coef = so_bc_coefs[bdry_normal_axis];
                    bc_coef.setBcCoefs(acoef_data,
                                       bcoef_data,
                                       gcoef_data,
                                       Pointer<SideVariable<NDIM, double> >(),
                                       *patch,
                                       trimmed_bdry_box,
                                       /*fill_time*/ 0.0);

                    std::array<double, NDIM> posn;
                    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const hier::Index<NDIM>& i = it();
                        const SideIndex<NDIM> is(i, bdry_normal_axis, SideIndex<NDIM>::Lower);

                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d != bdry_normal_axis)
                            {
                                posn[d] = patch_x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                            }
                            else
                            {
                                posn[d] = patch_x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                            }
                        }
                        const double& X0 = posn[0];
                        const double& X1 = posn[1];
                        double U2 = std::numeric_limits<double>::signaling_NaN();
                        if (bdry_normal_axis == 0)
                        {
                            U2 = (-3 * std::pow(X0, 2) * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                  2 * X1 * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                  2 * X0 * (std::pow(X0, 3) + std::pow(X1, 3)) -
                                  3 * std::pow(X1, 2) * (std::pow(X0, 3) + std::pow(X1, 3))) /
                                 2;
                        }

                        if (bdry_normal_axis == 1)
                        {
                            U2 = (-2 * X0 * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                  3 * std::pow(X1, 2) * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                  3 * std::pow(X0, 2) * (std::pow(X0, 3) + std::pow(X1, 3)) -
                                  2 * X1 * (std::pow(X0, 3) + std::pow(X1, 3))) /
                                 2;
                        }

                        const double error = std::abs((*gcoef_data)(i, 0) - U2);
                        double area = 1.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d != bdry_normal_axis) area *= dx[d];
                        }
                        weighted_sum += error * area;
                        normal_area += area;
                    }
                }
            }
        }

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const tbox::Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                if (n_physical_codim1_boxes == 0) continue;

                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const patch_x_lower = pgeom->getXLower();
                const double* const patch_x_upper = pgeom->getXUpper();
                const Box<NDIM>& patch_box = patch->getBox();
                const hier::Index<NDIM>& patch_lower = patch_box.lower();

                const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
                tbox::Array<tbox::Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    touches_regular_bdry[axis].resizeArray(2);
                    touches_periodic_bdry[axis].resizeArray(2);
                    for (int upperlower = 0; upperlower < 2; ++upperlower)
                    {
                        touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                        touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
                    }
                }

                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const unsigned int location_index = bdry_box.getLocationIndex();
                    const unsigned int bdry_normal_axis = location_index / 2;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, /*gcw_to_fill*/ 1);
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        if (axis != bdry_normal_axis)
                        {
                            Box<NDIM> bc_coef_box = compute_tangential_extension(
                                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                            Pointer<ArrayData<NDIM, double> > error_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                            // Temporarily reset the patch geometry object associated with
                            // the patch so that boundary conditions are set at the correct
                            // spatial locations.
                            std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                shifted_patch_x_lower[d] = patch_x_lower[d];
                                shifted_patch_x_upper[d] = patch_x_upper[d];
                            }
                            shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                            shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                            patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                                     touches_regular_bdry,
                                                                                     touches_periodic_bdry,
                                                                                     dx,
                                                                                     shifted_patch_x_lower.data(),
                                                                                     shifted_patch_x_upper.data()));

                            auto& bc_coef = so_bc_coefs[axis];
                            bc_coef.setBcCoefs(acoef_data,
                                               bcoef_data,
                                               gcoef_data,
                                               Pointer<SideVariable<NDIM, double> >(),
                                               *patch,
                                               trimmed_bdry_box,
                                               /*fill_time*/ 0.0);

                            std::array<double, NDIM> posn;
                            for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                            {
                                const hier::Index<NDIM>& i = b();
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d != bdry_normal_axis)
                                    {
                                        posn[d] = shifted_patch_x_lower[d] +
                                                  dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                                    }
                                    else
                                    {
                                        posn[d] = shifted_patch_x_lower[d] +
                                                  dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                                    }
                                }
                                const double& X0 = posn[0];
                                const double& X1 = posn[1];
                                double U2 = std::numeric_limits<double>::signaling_NaN();
                                if (axis == 0)
                                {
                                    U2 = (-3 * std::pow(X0, 2) * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                          2 * X1 * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                          2 * X0 * (std::pow(X0, 3) + std::pow(X1, 3)) -
                                          3 * std::pow(X1, 2) * (std::pow(X0, 3) + std::pow(X1, 3))) /
                                         2;
                                }

                                if (axis == 1)
                                {
                                    U2 = (-2 * X0 * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                          3 * std::pow(X1, 2) * (std::pow(X0, 2) + std::pow(X1, 2)) +
                                          3 * std::pow(X0, 2) * (std::pow(X0, 3) + std::pow(X1, 3)) -
                                          2 * X1 * (std::pow(X0, 3) + std::pow(X1, 3))) /
                                         2;
                                }

                                const double error = std::abs((*gcoef_data)(i, 0) - U2);
                                double area = 1.0;
                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    if (d != bdry_normal_axis) area *= dx[d];
                                }
                                weighted_sum += error * area;
                                tangential_area += area;
                            }
                        }
                    }
                }
                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);
            }
        }

        weighted_sum = IBTK_MPI::sumReduction(weighted_sum);
        normal_area = IBTK_MPI::sumReduction(normal_area);
        tangential_area = IBTK_MPI::sumReduction(tangential_area);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "total error is = " << weighted_sum << "\n"
                << "normal surface area is = " << normal_area << "\n"
                << "tangential surface area is = " << tangential_area << "\n";
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
