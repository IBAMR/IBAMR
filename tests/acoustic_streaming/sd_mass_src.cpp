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
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ACOUSTIC_SD_MASS_COUPLING_FC IBAMR_FC_FUNC(acoustic_sd_mass_coupling_2d, ACOUSTIC_SD_MASS_COUPLING_2D)
#endif
#if (NDIM == 3)
#define ACOUSTIC_SD_MASS_COUPLING_FC IBAMR_FC_FUNC(acoustic_sd_mass_coupling_3d, ACOUSTIC_SD_MASS_COUPLING_3D)
#endif

extern "C"
{
    void ACOUSTIC_SD_MASS_COUPLING_FC(const double* U0_real,
                                      const double* U1_real,
#if (NDIM == 3)
                                      const double* U2_real,
#endif
                                      const int& U_real_gcw,
                                      const double* U0_imag,
                                      const double* U1_imag,
#if (NDIM == 3)
                                      const double* U2_imag,
#endif
                                      const int& U_imag_gcw,
                                      const double* rho0,
                                      const double* rho1,
#if (NDIM == 3)
                                      const double* rho2,
#endif
                                      const int& rho_gcw,
                                      const double& omega,
                                      double* m,
                                      const int& m_gcw,
                                      const int& ilower0,
                                      const int& iupper0,
                                      const int& ilower1,
                                      const int& iupper1,
#if (NDIM == 3)
                                      const int& ilower2,
                                      const int& iupper2,
#endif
                                      const double* dx);
}

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
        Pointer<SideVariable<NDIM, double> > rho_var = new SideVariable<NDIM, double>("rho0", /*depth*/ 1);
        Pointer<CellVariable<NDIM, double> > mass_var = new CellVariable<NDIM, double>("mass", /*depth*/ 1);
        Pointer<CellVariable<NDIM, double> > mass_err_var = new CellVariable<NDIM, double>("mass_err", /*depth*/ 1);
        const int gcw = 2;
        const int U1_real_idx = var_db->registerVariableAndContext(U1_real_var, context, gcw);
        const int U1_imag_idx = var_db->registerVariableAndContext(U1_imag_var, context, gcw);
        const int rho_idx = var_db->registerVariableAndContext(rho_var, context, gcw);
        const int mass_idx = var_db->registerVariableAndContext(mass_var, context, /*ghost_width*/ 0);
        const int mass_err_idx = var_db->registerVariableAndContext(mass_err_var, context, /*ghost_width*/ 0);

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

        // Loop over patch box and fill the data
        // U1Real = "X0^3 + X1^3"
        // V1Real = "X1^2 + X0^2"
        // U1Imag = "(X0^2 + X1^2)"
        // V1Imag = "(X1^3 + X0^3)"
        // rho0 = "rho0 + rho1(X0^2*X1)"
        const double rho0 = input_db->getDouble("RHO0");
        const double rho1 = input_db->getDouble("RHO1");
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(U1_real_idx);
            level->allocatePatchData(U1_imag_idx);
            level->allocatePatchData(rho_idx);
            level->allocatePatchData(mass_idx);
            level->allocatePatchData(mass_err_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > U1_real_data = patch->getPatchData(U1_real_idx);
                Pointer<SideData<NDIM, double> > U1_imag_data = patch->getPatchData(U1_imag_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
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

                        (*rho_data)(i) = rho0 + rho1 * (X0 * X0 * X1);
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
        std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs(NDIM, nullptr);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const std::string bc_coefs_name = "rho_bc_coefs_" + std::to_string(d);
            const std::string bc_coefs_db_name = "DensityBcCoefs_" + std::to_string(d);

            rho_bc_coefs[d] = new muParserRobinBcCoefs(
                bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }

        // Fill the ghost cells of U1_real and U1_imag
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> comp_transactions(3);
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
        comp_transactions[2] = InterpolationTransactionComponent(
            rho_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "CONSTANT", false, rho_bc_coefs);
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
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
                const Box<NDIM>& ghost_box = U1_real_data->getGhostBox();

                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> ic(ghost_box, axis); ic; ic++)
                    {
                        const SideIndex<NDIM>& i = ic();
                        TBOX_ASSERT(std::abs((*U1_real_data)(i)) < 1e300);
                        TBOX_ASSERT(std::abs((*U1_imag_data)(i)) < 1e300);
                        TBOX_ASSERT(std::abs((*rho_data)(i)) < 1e300);
                    }
                }
            }
        }

        // Compute the mass source term using Stokes drift velocity div (rho U_SD)
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();

                Pointer<SideData<NDIM, double> > U1_real_data = patch->getPatchData(U1_real_idx);
                Pointer<SideData<NDIM, double> > U1_imag_data = patch->getPatchData(U1_imag_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
                Pointer<CellData<NDIM, double> > m_data = patch->getPatchData(mass_idx);

                const int U_real_gcw = (U1_real_data->getGhostCellWidth()).max();
                const int U_imag_gcw = (U1_imag_data->getGhostCellWidth()).max();
                const int m_gcw = (m_data->getGhostCellWidth()).max();
                const int rho_gcw = (rho_data->getGhostCellWidth()).max();

                const double* U0_real = U1_real_data->getPointer(0, 0);
                const double* U1_real = U1_real_data->getPointer(1, 0);
#if (NDIM == 3)
                const double* U2_real = U1_real_data->getPointer(2, 0);
#endif

                const double* U0_imag = U1_imag_data->getPointer(0, 0);
                const double* U1_imag = U1_imag_data->getPointer(1, 0);
#if (NDIM == 3)
                const double* U2_imag = U1_imag_data->getPointer(2, 0);
#endif

                const double* rho0 = rho_data->getPointer(0, 0);
                const double* rho1 = rho_data->getPointer(1, 0);
#if (NDIM == 3)
                const double* rho2 = rho_data->getPointer(2, 0);
#endif

                double* m = m_data->getPointer(0);

                ACOUSTIC_SD_MASS_COUPLING_FC(U0_real,
                                             U1_real,
#if (NDIM == 3)
                                             U2_real,
#endif
                                             U_real_gcw,
                                             U0_imag,
                                             U1_imag,
#if (NDIM == 3)
                                             U2_imag,
#endif
                                             U_imag_gcw,
                                             rho0,
                                             rho1,
#if (NDIM == 3)
                                             rho2,
#endif
                                             rho_gcw,
                                             /*acoustic_freq*/ 1.0,
                                             m,
                                             m_gcw,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2),
                                             patch_box.upper(2),
#endif
                                             dx);
            }
        }

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy, coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        muParserCartGridFunction mass_src_exact(
            "mass_src_exact", app_initializer->getComponentDatabase("MassSrc_fcn"), grid_geometry);
        mass_src_exact.setDataOnPatchHierarchy(mass_err_idx, mass_err_var, patch_hierarchy, /*time*/ 0.0);

        std::vector<double> mass_err(3);
        hier_cc_data_ops.subtract(mass_err_idx, mass_err_idx, mass_idx);
        mass_err[0] = hier_cc_data_ops.L1Norm(mass_err_idx, wgt_cc_idx);
        mass_err[1] = hier_cc_data_ops.L2Norm(mass_err_idx, wgt_cc_idx);
        mass_err[2] = hier_cc_data_ops.maxNorm(mass_err_idx, wgt_cc_idx);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "Error in mass src at time " << 0.0 << ":\n"
                << "  L1-norm:  " << std::setprecision(10) << mass_err[0] << "\n"
                << "  L2-norm:  " << mass_err[1] << "\n"
                << "  max-norm: " << mass_err[2] << "\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
