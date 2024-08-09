// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AcousticStreamingPETScMatUtilities.h"
#include "ibamr/AcousticStreamingPETScVecUtilities.h"
#include "ibamr/FOAcousticStreamingPETScLevelSolver.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"

#include "BoundaryBox.h"
#include "CellData.h"
#include "CellVariable.h"
#include "CoarseFineBoundary.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscvec.h"
#include <petsclog.h>

#include <algorithm>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

static const int REAL = 0;
static const int IMAG = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FOAcousticStreamingPETScLevelSolver::FOAcousticStreamingPETScLevelSolver(const std::string& object_name,
                                                                         Pointer<Database> input_db,
                                                                         const std::string& default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    if (input_db)
    {
        if (input_db->keyExists("acoustic_frequency")) d_omega = 2.0 * M_PI * input_db->getDouble("acoustic_frequency");
        if (input_db->keyExists("frequency")) d_omega = 2.0 * M_PI * input_db->getDouble("frequency");
        if (input_db->keyExists("acoustic_angular_frequency"))
            d_omega = input_db->getDouble("acoustic_angular_frequency");
        if (input_db->keyExists("angular_frequency")) d_omega = input_db->getDouble("angular_frequency");
        if (input_db->keyExists("sound_speed")) d_sound_speed = input_db->getDouble("sound_speed");
    }

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");

    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index", /*depth*/ 2);
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), SIDEG);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);

    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index", /*depth*/ 2);
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), CELLG);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    return;
} // FOAcousticStreamingPETScLevelSolver

FOAcousticStreamingPETScLevelSolver::~FOAcousticStreamingPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~FOAcousticStreamingPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FOAcousticStreamingPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                      const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    AcousticStreamingPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);

    AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp(d_petsc_mat,
                                                                                 d_omega,
                                                                                 d_sound_speed,
                                                                                 d_rho_idx,
                                                                                 d_mu_idx,
                                                                                 d_lambda_idx,
                                                                                 d_U_bc_coefs,
                                                                                 d_new_time,
                                                                                 d_num_dofs_per_proc,
                                                                                 d_u_dof_index_idx,
                                                                                 d_p_dof_index_idx,
                                                                                 d_level,
                                                                                 d_mu_interp_type);

    d_petsc_pc = d_petsc_mat;

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = AcousticStreamingPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    d_ghost_fill_sched = AcousticStreamingPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
FOAcousticStreamingPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    return;
} // deallocateSolverStateSpecialized

void
FOAcousticStreamingPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    AcousticStreamingPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
FOAcousticStreamingPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    AcousticStreamingPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

void
FOAcousticStreamingPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                                  Vec& petsc_b,
                                                  SAMRAIVectorReal<NDIM, double>& x,
                                                  SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int f_idx = b.getComponentDescriptorIndex(0);
    const int h_idx = b.getComponentDescriptorIndex(1);
    const auto f_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(f_idx);
    const auto h_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(h_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
        Pointer<SideData<NDIM, double> > f_adj_data = patch->getPatchData(f_adj_idx);
        Pointer<CellData<NDIM, double> > h_adj_data = patch->getPatchData(h_adj_idx);
        f_adj_data->copy(*f_data);
        h_adj_data->copy(*h_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        if (at_physical_bdry)
        {
            for (int comp = 0; comp < 2; ++comp)
            {
                const int other_comp = (comp == REAL ? IMAG : REAL);

                PoissonUtilities::adjustVCSCViscousDilatationalOpRHSAtPhysicalBoundary(*f_adj_data,
                                                                                       comp,
                                                                                       patch,
                                                                                       d_mu_idx,
                                                                                       d_lambda_idx,
                                                                                       d_U_bc_coefs[other_comp],
                                                                                       d_solution_time,
                                                                                       d_homogeneous_bc,
                                                                                       d_mu_interp_type);

                enforceNormalVelocityBoundaryConditions(
                    f_adj_idx, comp, patch, d_U_bc_coefs[comp], d_solution_time, d_homogeneous_bc);
            }
        }
        const Array<BoundaryBox<NDIM> >& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM> >() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            for (int comp = 0; comp < 2; ++comp)
            {
                PoissonUtilities::adjustVCSCViscousDilatationalOpRHSAtCoarseFineBoundary(
                    *f_adj_data, *u_data, comp, patch, d_mu_idx, d_lambda_idx, type_1_cf_bdry, d_mu_interp_type);
            }
        }
    }

    AcousticStreamingPETScVecUtilities::copyToPatchLevelVec(
        petsc_b, f_adj_idx, d_u_dof_index_idx, h_adj_idx, d_p_dof_index_idx, d_level);

    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FOAcousticStreamingPETScLevelSolver::enforceNormalVelocityBoundaryConditions(
    const int u_data_idx,
    const int data_depth,
    Pointer<Patch<NDIM> > patch,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const bool homogeneous_bc)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif

    Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_data_idx);

    // Data structures required to set physical boundary conditions.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;

        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        u_bc_coefs[bdry_normal_axis]->setBcCoefs(
            acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM> >(), *patch, trimmed_bdry_box, fill_time);
        auto const extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
        if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            const double& alpha = (*acoef_data)(i, 0);
            const double gamma = homogeneous_bc && !extended_bc_coef ? 0.0 : (*gcoef_data)(i, 0);
#if !defined(NDEBUG)
            const double& beta = (*bcoef_data)(i, 0);
            TBOX_ASSERT(IBTK::rel_equal_eps(alpha + beta, 1.0));
            TBOX_ASSERT(IBTK::rel_equal_eps(alpha, 1.0) || IBTK::rel_equal_eps(beta, 1.0));
#endif
            if (IBTK::rel_equal_eps(alpha, 1.0))
                (*u_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower), data_depth) = gamma;
        }
    }
    return;
} // enforceNormalVelocityBoundaryConditions

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
