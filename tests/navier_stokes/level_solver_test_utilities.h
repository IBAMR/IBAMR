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

#ifndef included_IBAMR_tests_level_solver_test_utilities
#define included_IBAMR_tests_level_solver_test_utilities

// Fixtures shared by the tests of PETScLevelSolver and its subclasses: a one-level periodic hierarchy with coupled
// velocity-pressure DOF indices, a probe that exposes the protected shell storage of a level solver, and a solve
// check on the installed KSP.

#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>

#include <tbox/MemoryDatabase.h>

#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <PatchHierarchy.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace level_solver_test
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

// Expose existing protected state for lifecycle assertions, without changing
// the production API or replacing any solver operation.
template <class Solver>
class LevelSolverProbe : public Solver
{
public:
    LevelSolverProbe(const std::string& name, Pointer<Database> db) : Solver(name, db, "")
    {
    }
    Mat matrixBeforeKSP() const
    {
        return d_matrix_before_ksp;
    }
    PetscInt referencesBeforeKSP() const
    {
        return d_references_before_ksp;
    }
    bool shellStorageEmpty() const
    {
        return this->d_sub_x.empty() && this->d_sub_y.empty() && this->d_sub_ksp.empty() &&
               this->d_restriction.empty() && this->d_prolongation.empty() && this->d_sub_mat == nullptr;
    }
    std::vector<Vec> retainShellVectors()
    {
        std::vector<Vec> result = this->d_sub_x;
        result.insert(result.end(), this->d_sub_y.begin(), this->d_sub_y.end());
        for (Vec v : result)
        {
            PetscErrorCode ierr = PetscObjectReference(reinterpret_cast<PetscObject>(v));
            IBTK_CHKERRQ(ierr);
        }
        return result;
    }

protected:
    void initializeSolverStateSpecialized(const HierarchyVector& x, const HierarchyVector& b) override
    {
        Solver::initializeSolverStateSpecialized(x, b);
        d_matrix_before_ksp = this->d_petsc_mat;
        PetscErrorCode ierr =
            PetscObjectGetReference(reinterpret_cast<PetscObject>(d_matrix_before_ksp), &d_references_before_ksp);
        IBTK_CHKERRQ(ierr);
    }

private:
    Mat d_matrix_before_ksp = nullptr;
    PetscInt d_references_before_ksp = 0;
};

struct LevelFixture
{
    Pointer<PatchHierarchy<NDIM>> hierarchy;
    Pointer<PatchLevel<NDIM>> level;
    Pointer<HierarchyVector> x, b;
    std::vector<int> indices;
    std::vector<PetscInt> velocity_ids;
    std::vector<int> dof_counts, velocity_counts;
    int full_size = 0;

    explicit LevelFixture(Pointer<Database> geometry_db, int ln = 0, bool full = true, bool distributed = false)
    {
        if (IBTK_MPI::getNodes() != (distributed ? 2 : 1))
        {
            TBOX_ERROR("Unexpected rank count for level fixture\n");
        }
        Pointer<CartesianGridGeometry<NDIM>> geometry = new CartesianGridGeometry<NDIM>("level_geometry", geometry_db);
        hierarchy = new PatchHierarchy<NDIM>("level_hierarchy", geometry);
        const int ranks = IBTK_MPI::getNodes();
        BoxArray<NDIM> boxes(ranks);
        ProcessorMapping mapping(ranks);
        for (int rank = 0; rank < ranks; ++rank)
        {
            SAMRAI::hier::Index<NDIM> lower(0), upper(15);
            lower(0) = rank * 16 / ranks;
            upper(0) = (rank + 1) * 16 / ranks - 1;
            boxes[rank] = Box<NDIM>(lower, upper);
            mapping.setProcessorAssignment(rank, rank);
        }
        hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), boxes, mapping);
        if (ln == 1)
        {
            boxes[0] = full ? Box<NDIM>(SAMRAI::hier::Index<NDIM>(0), SAMRAI::hier::Index<NDIM>(31)) :
                              Box<NDIM>(SAMRAI::hier::Index<NDIM>(8), SAMRAI::hier::Index<NDIM>(23));
            hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), boxes, mapping);
        }
        level = hierarchy->getPatchLevel(ln);
        VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = db->getContext("level_fixture");
        Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("level_u");
        Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("level_p");
        if (db->checkVariableExists("level_u"))
        {
            u = db->getVariable("level_u");
        }
        if (db->checkVariableExists("level_p"))
        {
            p = db->getVariable("level_p");
        }
        const int ui = db->registerVariableAndContext(u, context, IntVector<NDIM>(1));
        const int pi = db->registerVariableAndContext(p, context, IntVector<NDIM>(1));
        x = new HierarchyVector("level_x", hierarchy, ln, ln);
        x->addComponent(u, ui, -1, new HierarchySideDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->addComponent(p, pi, -1, new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->allocateVectorData();
        x->setToScalar(0.0);
        b = x->cloneVector("level_b");
        b->allocateVectorData();
        b->setToScalar(0.0);
        Pointer<SideVariable<NDIM, int>> ud = new SideVariable<NDIM, int>("level_ud");
        Pointer<CellVariable<NDIM, int>> pd = new CellVariable<NDIM, int>("level_pd");
        if (db->checkVariableExists("level_ud"))
        {
            ud = db->getVariable("level_ud");
        }
        if (db->checkVariableExists("level_pd"))
        {
            pd = db->getVariable("level_pd");
        }
        const int udi = db->registerVariableAndContext(ud, context, IntVector<NDIM>(1));
        const int pdi = db->registerVariableAndContext(pd, context, IntVector<NDIM>(1));
        indices = { udi, pdi };
        for (int idx : indices)
        {
            level->allocatePatchData(idx);
        }
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(dof_counts, udi, pdi, level);
        const int rank = IBTK_MPI::getRank();
        int first_dof = 0;
        for (int r = 0; r < ranks; ++r)
        {
            full_size += dof_counts[r];
            if (r < rank)
            {
                first_dof += dof_counts[r];
            }
        }
        std::set<int> velocity;
        Pointer<SideData<NDIM, int>> data = level->getPatch(rank)->getPatchData(udi);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator i(SideGeometry<NDIM>::toSideBox(level->getPatch(rank)->getBox(), axis)); i; i++)
            {
                const int id = (*data)(SideIndex<NDIM>(i(), axis, SideIndex<NDIM>::Lower));
                if (id >= first_dof && id < first_dof + dof_counts[rank])
                {
                    velocity.insert(id);
                }
            }
        }
        velocity_ids.assign(velocity.begin(), velocity.end());
        const int local_velocity = velocity_ids.size();
        velocity_counts.resize(ranks);
        int ierr = MPI_Allgather(&local_velocity, 1, MPI_INT, velocity_counts.data(), 1, MPI_INT, PETSC_COMM_WORLD);
        IBTK_CHKERRQ(ierr);
        std::vector<int> offsets(ranks, 0);
        for (int r = 1; r < ranks; ++r)
        {
            offsets[r] = offsets[r - 1] + velocity_counts[r - 1];
        }
        std::vector<PetscInt> all_velocity(offsets.back() + velocity_counts.back());
        ierr = MPI_Allgatherv(velocity_ids.data(),
                              local_velocity,
                              MPIU_INT,
                              all_velocity.data(),
                              velocity_counts.data(),
                              offsets.data(),
                              MPIU_INT,
                              PETSC_COMM_WORLD);
        IBTK_CHKERRQ(ierr);
        velocity_ids = std::move(all_velocity);
    }
    ~LevelFixture()
    {
        free_vector_components(*b);
        free_vector_components(*x);
        for (int idx : indices)
        {
            level->deallocatePatchData(idx);
            VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(idx);
        }
    }
};

inline Pointer<Database>
level_solver_database(const std::string& pc = "none", int overlap = 0)
{
    Pointer<Database> db = new MemoryDatabase("level_solver");
    db->putString("ksp_type", "gmres");
    db->putString("options_prefix", "level_");
    db->putString("pc_type", pc);
    db->putString("shell_pc_type", "additive");
    db->putBool("initial_guess_nonzero", false);
    db->putDouble("rel_residual_tol", 1.0e-12);
    db->putInteger("max_iterations", 100);
    const int size[NDIM] = { 4, 8 }, width[NDIM] = { overlap, overlap };
    db->putIntegerArray("subdomain_box_size", size, NDIM);
    db->putIntegerArray("subdomain_overlap_size", width, NDIM);
    return db;
}

// The largest solution error of check_level_solve().
inline double max_level_solve_error = 0.0;

// Exercise the installed solver KSP and its actual configured matrix/PC.
inline bool
check_level_solve(PETScLevelSolver& solver)
{
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, solution;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &solution);
    IBTK_CHKERRQ(ierr);
    PetscInt first, last;
    ierr = VecGetOwnershipRange(exact, &first, &last);
    IBTK_CHKERRQ(ierr);
    PetscScalar* values;
    ierr = VecGetArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = first; i < last; ++i)
    {
        values[i - first] = std::sin(0.13 * i) + 0.5;
    }
    ierr = VecRestoreArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSolve(solver.getPETScKSP(), rhs, solution);
    IBTK_CHKERRQ(ierr);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(solver.getPETScKSP(), &reason);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(solution, -1.0, exact);
    IBTK_CHKERRQ(ierr);
    PetscReal error;
    ierr = VecNorm(solution, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &solution })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    // std::max keeps its first argument when the second is NaN, so a NaN error must make the maximum infinite.
    max_level_solve_error =
        std::isfinite(error) ? std::max<double>(max_level_solve_error, error) : std::numeric_limits<double>::infinity();
    return reason > 0 && std::isfinite(error) && error < 1.0e-9;
}
} // namespace level_solver_test

#endif
