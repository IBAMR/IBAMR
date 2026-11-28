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

#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonBoxRelaxationFACOperator.h>
#include <ibtk/CCPoissonLevelRelaxationFACOperator.h>
#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PoissonFACPreconditionerStrategy.h>
#include <ibtk/SCPoissonPointRelaxationFACOperator.h>
#include <ibtk/SideNoCornersFillPattern.h>

#include <tbox/MemoryDatabase.h>
#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchyDataOpsReal.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <ProcessorMapping.h>
#include <SAMRAIVectorReal.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../coarse_fine_ghost_values.h"

#include <ibtk/app_namespaces.h>

using namespace coarse_fine_ghost_values;

// Compare the coarse-fine ghost values that the smoothers of the FAC preconditioner strategies set with the values that
// the ghost cell fill of the residual sets, after a prolongation and a post-smoothing sweep and after a pre-smoothing
// sweep from zero followed by a post-smoothing sweep.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "fac_coarse_fine_ghost_values.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        const bool is_cell = input_db->getString("CENTERING") == "CELL";
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        const int depth = input_db->getIntegerWithDefault("DEPTH", 1);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("CTX");
        Pointer<Variable<NDIM>> u_var, f_var, w_var;
        if (is_cell)
        {
            u_var = new CellVariable<NDIM, double>("u", depth);
            f_var = new CellVariable<NDIM, double>("f", depth);
            w_var = new CellVariable<NDIM, double>("w", depth);
        }
        else
        {
            u_var = new SideVariable<NDIM, double>("u", depth);
            f_var = new SideVariable<NDIM, double>("f", depth);
            w_var = new SideVariable<NDIM, double>("w", depth);
        }
        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(ghost_width));
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(ghost_width));
        const int w_idx = var_db->registerVariableAndContext(w_var, ctx, IntVector<NDIM>(ghost_width));

        // The fine level has the given boxes, spread over the processes.
        gridding_algorithm->makeCoarsestLevel(hierarchy, 0.0);
        const tbox::Array<tbox::DatabaseBox> db_boxes = input_db->getDatabaseBoxArray("FINE_BOXES");
        const tbox::Array<int> ratio_values = input_db->getIntegerArray("REF_RATIO");
        IntVector<NDIM> ratio;
        for (int d = 0; d < NDIM; ++d)
        {
            ratio(d) = ratio_values[ratio_values.size() == 1 ? 0 : d];
        }
        BoxArray<NDIM> fine_boxes(db_boxes.size());
        ProcessorMapping mapping(db_boxes.size());
        for (int k = 0; k < db_boxes.size(); ++k)
        {
            fine_boxes[k] = Box<NDIM>(db_boxes[k]);
            mapping.setProcessorAssignment(k, k % IBTK_MPI::getNodes());
        }
        hierarchy->makeNewPatchLevel(1, ratio, fine_boxes, mapping);
        for (int ln = 0; ln <= 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (const int idx : { u_idx, f_idx, w_idx })
            {
                level->allocatePatchData(idx);
            }
        }

        // The error satisfies homogeneous Dirichlet conditions on the physical boundary.
        std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs;
        for (int d = 0; d < (is_cell ? depth : NDIM); ++d)
        {
            bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc_coefs", nullptr));
            for (int location = 0; location < 2 * NDIM; ++location)
            {
                bc_storage.back()->setBoundaryValue(location, 0.0);
            }
            bc_coefs.push_back(bc_storage.back().get());
        }
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCZero();
        poisson_spec.setDConstant(1.0);

        SAMRAIVectorReal<NDIM, double> u_vec("u_vec", hierarchy, 0, 1), f_vec("f_vec", hierarchy, 0, 1);
        u_vec.addComponent(u_var, u_idx);
        f_vec.addComponent(f_var, f_idx);
        Pointer<HierarchyDataOpsReal<NDIM, double>> ops =
            HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(u_var, hierarchy, true);

        // The operators, with their names.
        Pointer<Database> op_db = new MemoryDatabase("op_db");
        op_db->putString("coarse_solver_type", "PATCH_GAUSS_SEIDEL");
        Pointer<Database> level_op_db = new MemoryDatabase("level_op_db");
        std::vector<std::pair<std::string, Pointer<PoissonFACPreconditionerStrategy>>> operators;
        if (is_cell)
        {
            operators.emplace_back("point relaxation",
                                   new CCPoissonPointRelaxationFACOperator("cc_point", op_db, "cc_point_"));
            operators.emplace_back("box relaxation", new CCPoissonBoxRelaxationFACOperator("cc_box", op_db, "cc_box_"));
            operators.emplace_back("level relaxation",
                                   new CCPoissonLevelRelaxationFACOperator("cc_level", level_op_db, "cc_level_"));
        }
        else
        {
            operators.emplace_back("point relaxation",
                                   new SCPoissonPointRelaxationFACOperator("sc_point", op_db, "sc_point_"));
        }

        // The ghost cell fill of the residual of the operators, without coarsening.
        Pointer<VariableFillPattern<NDIM>> fill_pattern;
        if (is_cell)
        {
            fill_pattern = new CellNoCornersFillPattern(1, /*overwrite_interior*/ false);
        }
        else
        {
            fill_pattern = new SideNoCornersFillPattern(1, /*overwrite_interior*/ false);
        }
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        const InterpolationTransactionComponent w_comp(
            w_idx, "NONE", /*use_cf_bdry_interpolation*/ true, "NONE", "LINEAR", false, bc_coefs, fill_pattern);

        for (auto& entry : operators)
        {
            const std::string& name = entry.first;
            Pointer<PoissonFACPreconditionerStrategy> op = entry.second;
            op->setPoissonSpecifications(poisson_spec);
            op->setPhysicalBcCoefs(bc_coefs);
            op->setSolutionTime(0.0);
            op->setTimeInterval(0.0, 1.0);
            op->initializeOperatorState(u_vec, f_vec);
            HierarchyGhostCellInterpolation w_fill;
            w_fill.initializeOperatorState(w_comp, hierarchy, 0, 1);
            w_fill.setHomogeneousBc(true);

            // A prolongation and a post-smoothing sweep.
            set_coarse_fine_ghost_test_data(hierarchy, u_idx, 0.0);
            set_coarse_fine_ghost_test_data(hierarchy, f_idx, 1.0);
            op->prolongErrorAndCorrect(u_vec, u_vec, 1);
            ops->copyData(w_idx, u_idx, /*interior_only*/ false);
            op->smoothError(u_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ false, /*performing_post_sweeps*/ true);
            w_fill.fillData(0.0);
            plog << name << ": post-smoothing compares "
                 << compare_coarse_fine_ghost_values(hierarchy, 1, u_idx, w_idx, name + " post-smoothing")
                 << " coarse-fine ghost values\n";

            // Pre-smoothing from zero, then a post-smoothing sweep; the coarse correction is zero.
            u_vec.setToScalar(0.0, /*interior_only*/ false);
            set_coarse_fine_ghost_test_data(hierarchy, f_idx, 1.0);
            op->smoothError(u_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ true, /*performing_post_sweeps*/ false);
            ops->copyData(w_idx, u_idx, /*interior_only*/ false);
            op->smoothError(u_vec, f_vec, 1, 1, /*performing_pre_sweeps*/ false, /*performing_post_sweeps*/ true);
            w_fill.fillData(0.0);
            plog << name << ": zero-start smoothing compares "
                 << compare_coarse_fine_ghost_values(hierarchy, 1, u_idx, w_idx, name + " zero-start smoothing")
                 << " coarse-fine ghost values\n";

            w_fill.deallocateOperatorState();
            op->deallocateOperatorState();
        }
    }
    return 0;
}
