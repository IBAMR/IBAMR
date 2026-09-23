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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellVariable.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

// For each of the additive shell, multiplicative shell, and ASM preconditioners, reinitialize one
// level solver on two levels whose patch layouts differ, and compare it with a solver that is
// initialized only on the second level. The subdomains must be regenerated for the new layout.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(appender);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    if (hierarchy->getFinestLevelNumber() != 1)
    {
        TBOX_ERROR("The reinitialization test requires two levels.\n");
    }
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("reinitialization_test");
    Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h = new CellVariable<NDIM, double>("h");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("p_dof");
    // The reference solution needs its own data.
    Pointer<SideVariable<NDIM, double>> u_reference = new SideVariable<NDIM, double>("u_reference");
    Pointer<CellVariable<NDIM, double>> p_reference = new CellVariable<NDIM, double>("p_reference");
    const int uri = variables->registerVariableAndContext(u_reference, context, IntVector<NDIM>(1));
    const int pri = variables->registerVariableAndContext(p_reference, context, IntVector<NDIM>(1));
    const int ui = variables->registerVariableAndContext(u, context, IntVector<NDIM>(1));
    const int pi = variables->registerVariableAndContext(p, context, IntVector<NDIM>(1));
    const int fi = variables->registerVariableAndContext(f, context, IntVector<NDIM>(1));
    const int hi = variables->registerVariableAndContext(h, context, IntVector<NDIM>(1));
    const int udi = variables->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int pdi = variables->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    for (int ln = 0; ln <= 1; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (int index : { ui, pi, fi, hi, udi, pdi, uri, pri })
        {
            level->allocatePatchData(index);
        }
    }

    SAMRAIVectorReal<NDIM, double> x0("x0", hierarchy, 0, 0), b0("b0", hierarchy, 0, 0);
    SAMRAIVectorReal<NDIM, double> x1("x1", hierarchy, 1, 1), b1("b1", hierarchy, 1, 1);
    SAMRAIVectorReal<NDIM, double> x_reference("x_reference", hierarchy, 1, 1);
    for (SAMRAIVectorReal<NDIM, double>* x : { &x0, &x1 })
    {
        x->addComponent(u, ui);
        x->addComponent(p, pi);
    }
    x_reference.addComponent(u_reference, uri);
    x_reference.addComponent(p_reference, pri);
    for (SAMRAIVectorReal<NDIM, double>* b : { &b0, &b1 })
    {
        b->addComponent(f, fi);
        b->addComponent(h, hi);
    }
    b0.setToScalar(0.0);
    b1.setToScalar(0.0);
    for (int ln = 0; ln <= 1; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        const double wavenumber = 2.0 * std::acos(-1.0) / (ln == 0 ? 8.0 : 16.0);
        for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
            Pointer<SideData<NDIM, double>> force = patch->getPatchData(fi);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                for (Box<NDIM>::Iterator it(box); it; it++)
                {
                    (*force)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                        std::sin(wavenumber * it()(axis)) + 0.25 * (axis + 1.0);
                }
            }
        }
    }

    PoissonSpecifications coefficients("coefficients");
    coefficients.setCConstant(1.0);
    coefficients.setDConstant(-1.0);

    // Each preconditioner builds its subdomains from the level of the solver's current state.
    struct Preconditioner
    {
        const char* name;
        const char* pc_type;
        const char* shell_pc_type;
    };
    const std::vector<Preconditioner> preconditioners = { { "additive", "shell", "additive" },
                                                          { "multiplicative", "shell", "multiplicative" },
                                                          { "asm", "asm", "" } };
    plog << std::setprecision(12);
    for (const Preconditioner& preconditioner : preconditioners)
    {
        Pointer<MemoryDatabase> db = new MemoryDatabase("solver");
        db->putString("ksp_type", "preonly");
        db->putString("pc_type", preconditioner.pc_type);
        if (std::string(preconditioner.shell_pc_type).size() > 0)
        {
            db->putString("shell_pc_type", preconditioner.shell_pc_type);
        }
        db->putBool("initial_guess_nonzero", false);
        db->putInteger("max_iterations", 1);
        int box_size[NDIM];
        std::fill_n(box_size, NDIM, 4);
        db->putIntegerArray("subdomain_box_size", box_size, NDIM);
        const auto make_solver = [&](const std::string& name)
        {
            auto solver = std::make_unique<StaggeredStokesPETScLevelSolver>(name, db, "shell_");
            solver->setVelocityPoissonSpecifications(coefficients);
            solver->setComponentsHaveNullSpace(false, true);
            return solver;
        };

        std::unique_ptr<StaggeredStokesPETScLevelSolver> reused = make_solver("reused_solver");
        x0.setToScalar(0.0);
        reused->initializeSolverState(x0, b0);
        if (!reused->solveSystem(x0, b0))
        {
            TBOX_ERROR("Failed check: !reused->solveSystem(x0, b0).\n");
        }
        reused->deallocateSolverState();

        x1.setToScalar(0.0);
        reused->initializeSolverState(x1, b1);
        if (!reused->solveSystem(x1, b1))
        {
            TBOX_ERROR("Failed check: !reused->solveSystem(x1, b1).\n");
        }
        reused->deallocateSolverState();

        std::unique_ptr<StaggeredStokesPETScLevelSolver> fresh = make_solver("fresh_solver");
        x_reference.setToScalar(0.0);
        fresh->initializeSolverState(x_reference, b1);
        if (!fresh->solveSystem(x_reference, b1))
        {
            TBOX_ERROR("Failed check: !fresh->solveSystem(x_reference, b1).\n");
        }
        fresh->deallocateSolverState();

        const double solution_norm = x1.maxNorm();
        x1.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&x1, false),
                    Pointer<SAMRAIVectorReal<NDIM, double>>(&x_reference, false));
        const double error = x1.maxNorm();
        if (!std::isfinite(error) || error > 1.0e-10 || !(solution_norm > 0.0))
        {
            TBOX_ERROR("Failed check: !std::isfinite(error) || error > 1.0e-10 || !(solution_norm > 0.0).\n");
        }
        plog << preconditioner.name << " solution_norm = " << solution_norm << "\n"
             << preconditioner.name << " error = " << error << '\n';
    }
    return 0;
}
