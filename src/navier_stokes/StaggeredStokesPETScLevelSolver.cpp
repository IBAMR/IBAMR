// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StokesBcCoefStrategy.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/GeneralSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/ShiftedPatchGeometry.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petsclog.h>
#include <petscvec.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchDescriptor.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "./ins_staggered_traction_stencil.h"
#include "./ins_velocity_bc.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NOGHOST = 0;

// Subtract from the right-hand side of each velocity side on a coarse-fine interface the term of the pressure in the
// cell across the interface, which has no degree of freedom.
void
adjust_rhs_for_pressure_at_cf_boundary(SideData<NDIM, double>& rhs_data,
                                       const CellData<NDIM, double>& p_data,
                                       const CellData<NDIM, int>& p_dof_index_data,
                                       Pointer<Patch<NDIM>> patch,
                                       const Array<BoundaryBox<NDIM>>& type_1_cf_bdry)
{
    TBOX_ASSERT(p_data.getGhostCellWidth().min() >= 1);
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (int n = 0; n < type_1_cf_bdry.size(); ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = type_1_cf_bdry[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(trimmed_bdry_box, patch_box, IntVector<NDIM>(1));
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int axis = location_index / 2;
        const bool is_upper = location_index % 2 == 1;

        // The pressure gradient is +1/dx from the upper cell of a side and -1/dx from the lower cell (see
        // constructPatchLevelMACStokesOp()).
        const double coef = (is_upper ? +1.0 : -1.0) / dx[axis];
        for (Box<NDIM>::Iterator b(bc_fill_box); b; b++)
        {
            const CellIndex<NDIM> i_c_bdry(b());
            if (p_dof_index_data(i_c_bdry) >= 0)
            {
                continue;
            }
            hier::Index<NDIM> i_s = b();
            if (!is_upper)
            {
                i_s(axis) += 1;
            }
            rhs_data(SideIndex<NDIM>(i_s, axis, SideIndex<NDIM>::Lower)) -= coef * p_data(i_c_bdry);
        }
    }
    return;
} // adjust_rhs_for_pressure_at_cf_boundary

// Move the term f_g*g of the ghost pressure p_G = f_i*p_I + f_g*g, with interior pressure p_I and boundary data g, in
// the normal-velocity equations at open boundaries to the right-hand side.
void
adjust_rhs_for_ghost_pressure(SideData<NDIM, double>& rhs_data,
                              Patch<NDIM>& patch,
                              const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                              RobinBcCoefStrategy<NDIM>* P_bc_coef,
                              const double data_time,
                              const bool homogeneous_bc)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        // The normal velocity is not prescribed where the velocity boundary condition is a traction condition.
        Pointer<ArrayData<NDIM, double>> u_acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> u_bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> u_gcoef_data;
        U_bc_coefs[bdry_normal_axis]->setBcCoefs(
            u_acoef_data, u_bcoef_data, u_gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);

        // For homogeneous conditions, setupBcCoefObjects() makes an ExtendedRobinBcCoefStrategy return homogeneous
        // data. This function zeroes g of any other object, as the other users of Robin data do.
        Pointer<ArrayData<NDIM, double>> p_acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> p_bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> p_gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        auto extended_p_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(P_bc_coef);
        P_bc_coef->setBcCoefs(p_acoef_data, p_bcoef_data, p_gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);
        if (homogeneous_bc && !extended_p_bc_coef)
        {
            p_gcoef_data->fillAll(0.0);
        }

        const double h = dx[bdry_normal_axis];
        for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            if (velocity_bc::get_type((*u_acoef_data)(i, 0),
                                      (*u_bcoef_data)(i, 0),
                                      "StaggeredStokesPETScLevelSolver::setupKSPVecs()") != velocity_bc::Type::TRACTION)
            {
                continue;
            }
            const double f_g = 2.0 * h / ((*p_acoef_data)(i, 0) * h + 2.0 * (*p_bcoef_data)(i, 0));
            const double ghost_pressure_datum = f_g * (*p_gcoef_data)(i, 0);
            rhs_data(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) +=
                (is_lower ? 1.0 : -1.0) * ghost_pressure_datum / h;
        }
    }
    return;
} // adjust_rhs_for_ghost_pressure

// Move the gamma part of the normal-stress term D/h^2*(u_I - u_div) (traction_stencil::get_normal_stress_stencil()) at
// a TRACTION boundary to the right-hand side. A tangential ghost value is f_i*u + f_g*gamma, in which f_i and f_g are
// weights computed from the Robin coefficients and the grid spacing, u is the value on the face inside the domain, and
// gamma is the inhomogeneous coefficient. The matrix contains the u terms.
void
adjust_rhs_for_normal_stress(SideData<NDIM, double>& rhs_data,
                             Patch<NDIM>& patch,
                             const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                             const double D,
                             const double data_time,
                             const bool homogeneous_bc)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        const auto normal_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(U_bc_coefs[bdry_normal_axis]);
        if (!normal_bc_coef || normal_bc_coef->getTractionBcType() != TRACTION)
        {
            continue;
        }
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
        const Box<NDIM> normal_bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        // The boundary conditions of the tangential components at the locations of their ghost values.
        std::array<Pointer<ArrayData<NDIM, double>>, NDIM> acoef_data, bcoef_data, gcoef_data;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (axis == bdry_normal_axis)
            {
                continue;
            }
            Box<NDIM> bc_coef_box = normal_bc_coef_box;
            bc_coef_box.upper()(axis) += 1;
            acoef_data[axis] = new ArrayData<NDIM, double>(bc_coef_box, 1);
            bcoef_data[axis] = new ArrayData<NDIM, double>(bc_coef_box, 1);
            gcoef_data[axis] = new ArrayData<NDIM, double>(bc_coef_box, 1);
            std::array<double, NDIM> shift;
            shift.fill(0.0);
            shift[axis] = -0.5 * dx[axis];
            const ShiftedPatchGeometry shifted_geometry(patch, shift);
            const ShiftedPatchGeometry::ScopedGeometryReplacement geometry_replacement(patch, shifted_geometry);
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(U_bc_coefs[axis]);
            U_bc_coefs[axis]->setBcCoefs(
                acoef_data[axis], bcoef_data[axis], gcoef_data[axis], nullptr, patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef)
            {
                gcoef_data[axis]->fillAll(0.0);
            }
        }

        // The normal velocity is not prescribed where the velocity boundary condition is a traction condition.
        Pointer<ArrayData<NDIM, double>> normal_acoef_data = new ArrayData<NDIM, double>(normal_bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> normal_bcoef_data = new ArrayData<NDIM, double>(normal_bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> normal_gcoef_data;
        normal_bc_coef->setBcCoefs(
            normal_acoef_data, normal_bcoef_data, normal_gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);
        const double h = dx[bdry_normal_axis];
        for (Box<NDIM>::Iterator b(normal_bc_coef_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            if (velocity_bc::get_type((*normal_acoef_data)(i, 0),
                                      (*normal_bcoef_data)(i, 0),
                                      "StaggeredStokesPETScLevelSolver::setupKSPVecs()") != velocity_bc::Type::TRACTION)
            {
                continue;
            }
            double normal_stress_datum = 0.0;
            for (const auto& entry : traction_stencil::get_normal_stress_stencil(i, bdry_normal_axis, is_lower, dx))
            {
                const unsigned int axis = entry.first.getAxis();
                if (axis == bdry_normal_axis)
                {
                    continue;
                }
                hier::Index<NDIM> i_bc = entry.first;
                if (is_lower)
                {
                    i_bc(bdry_normal_axis) += 1;
                }
                const double a = (*acoef_data[axis])(i_bc, 0);
                const double b_coef = (*bcoef_data[axis])(i_bc, 0);
                const double f_g = 2.0 * h / (a * h + 2.0 * b_coef);
                normal_stress_datum += entry.second * f_g * (*gcoef_data[axis])(i_bc, 0);
            }
            rhs_data(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) -= D / (h * h) * normal_stress_datum;
        }
    }
    return;
} // adjust_rhs_for_normal_stress
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), SIDEG);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);
    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), CELLG);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    // Construct the nullspace variable/index.
    d_u_nullspace_var = new SideVariable<NDIM, double>(object_name + "::u_nullspace_var");
    if (var_db->checkVariableExists(d_u_nullspace_var->getName()))
    {
        d_u_nullspace_var = var_db->getVariable(d_u_nullspace_var->getName());
        d_u_nullspace_idx = var_db->mapVariableAndContextToIndex(d_u_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_u_nullspace_idx);
    }
    d_u_nullspace_idx = var_db->registerVariableAndContext(d_u_nullspace_var, d_context, NOGHOST);
    d_p_nullspace_var = new CellVariable<NDIM, double>(object_name + "::p_nullspace_var");
    if (var_db->checkVariableExists(d_p_nullspace_var->getName()))
    {
        d_p_nullspace_var = var_db->getVariable(d_p_nullspace_var->getName());
        d_p_nullspace_idx = var_db->mapVariableAndContextToIndex(d_p_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_p_nullspace_idx);
    }
    d_p_nullspace_idx = var_db->registerVariableAndContext(d_p_nullspace_var, d_context, NOGHOST);

    // Construct the variable/index of the velocity that the velocity boundary conditions use.
    d_u_bc_target_var = new SideVariable<NDIM, double>(object_name + "::u_bc_target_var");
    if (var_db->checkVariableExists(d_u_bc_target_var->getName()))
    {
        d_u_bc_target_var = var_db->getVariable(d_u_bc_target_var->getName());
        d_u_bc_target_idx = var_db->mapVariableAndContextToIndex(d_u_bc_target_var, d_context);
        var_db->removePatchDataIndex(d_u_bc_target_idx);
    }
    d_u_bc_target_idx = var_db->registerVariableAndContext(d_u_bc_target_var, d_context, SIDEG);

    return;
} // StaggeredStokesPETScLevelSolver

StaggeredStokesPETScLevelSolver::~StaggeredStokesPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~StaggeredStokesPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesPETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                                       std::vector<std::set<int>>& nonoverlap_is)
{
    // Construct subdomains for ASM and MSM preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(overlap_is,
                                                                       nonoverlap_is,
                                                                       d_box_size,
                                                                       d_overlap_size,
                                                                       d_num_dofs_per_proc,
                                                                       d_u_dof_index_idx,
                                                                       d_p_dof_index_idx,
                                                                       d_level,
                                                                       d_cf_boundary);

    return;
} // generateASMSubdomains

void
StaggeredStokesPETScLevelSolver::generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                                              std::vector<std::set<int>>& field_is)
{
    // Set IS'es for field split preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    return;
} // generateFieldSplitSubdomains

void
StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAIVectorReal<NDIM, double>& b)
{
    // Physical boundary DOFs need the boundary helper. Fully periodic hierarchies
    // do not require a boundary helper.
    if (!d_bc_helper)
    {
        for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                TBOX_ERROR(d_object_name << "::initializeSolverStateSpecialized():\n"
                                         << "  physical boundaries require a StaggeredStokesPhysicalBoundaryHelper "
                                         << "provided via setPhysicalBoundaryHelper().\n");
            }
        }
    }

    // The velocity and pressure data of both vectors must have a ghost layer. The data may be unallocated.
    const Pointer<PatchDescriptor<NDIM>> patch_descriptor = VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor();
    for (const SAMRAIVectorReal<NDIM, double>* const vec : { &x, &b })
    {
        for (int comp = 0; comp < 2; ++comp)
        {
            const int ghost_width = patch_descriptor->getPatchDataFactory(vec->getComponentDescriptorIndex(comp))
                                        ->getGhostCellWidth()
                                        .min();
            if (ghost_width < 1)
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                         << "  the " << (comp == 0 ? "velocity" : "pressure") << " data of the "
                                         << (vec == &x ? "solution" : "right-hand side")
                                         << " vector has ghost cell width " << ghost_width << "\n"
                                         << "  allocate the velocity and pressure data with a ghost cell width of at "
                                            "least one\n");
            }
        }
    }

    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    if (!d_level->checkAllocated(d_u_bc_target_idx))
    {
        d_level->allocatePatchData(d_u_bc_target_idx);
    }
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<SideData<NDIM, double>> u_bc_target_data = d_level->getPatch(p())->getPatchData(d_u_bc_target_idx);
        u_bc_target_data->fillAll(0.0);
    }
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_petsc_mat,
                                                                     d_U_problem_coefs,
                                                                     d_U_bc_coefs,
                                                                     d_new_time,
                                                                     d_num_dofs_per_proc,
                                                                     d_u_dof_index_idx,
                                                                     d_p_dof_index_idx,
                                                                     d_level,
                                                                     d_P_bc_coef);
    d_petsc_pc = d_petsc_mat;

    // Set pressure nullspace if the level covers the entire domain.
    if (d_has_pressure_nullspace)
    {
        bool level_covers_entire_domain = d_level_num == 0;
        if (d_level_num > 0)
        {
            int local_cf_bdry_box_size = 0;
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                              /* boundary type */ 1);
                local_cf_bdry_box_size += type_1_cf_bdry.size();
            }
            level_covers_entire_domain = IBTK_MPI::sumReduction(local_cf_bdry_box_size) == 0;
        }

        if (level_covers_entire_domain)
        {
            // Allocate pressure nullspace data.
            if (!d_level->checkAllocated(d_u_nullspace_idx)) d_level->allocatePatchData(d_u_nullspace_idx);
            if (!d_level->checkAllocated(d_p_nullspace_idx)) d_level->allocatePatchData(d_p_nullspace_idx);

            Pointer<SAMRAIVectorReal<NDIM, double>> nullspace_vec = new SAMRAIVectorReal<NDIM, double>(
                d_object_name + "nullspace_vec", d_hierarchy, d_level_num, d_level_num);
            nullspace_vec->addComponent(d_u_nullspace_var, d_u_nullspace_idx);
            nullspace_vec->addComponent(d_p_nullspace_var, d_p_nullspace_idx);
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_patch_data = nullspace_vec->getComponentPatchData(0, *patch);
                u_patch_data->fill(0.0);
                Pointer<CellData<NDIM, double>> p_patch_data = nullspace_vec->getComponentPatchData(1, *patch);
                p_patch_data->fill(1.0);
            }

            LinearSolver::setNullSpace(
                /*const vec*/ false, std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>(1, nullspace_vec));
        }
    }

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    if (d_level->checkAllocated(d_u_bc_target_idx))
    {
        d_level->deallocatePatchData(d_u_bc_target_idx);
    }
    return;
} // deallocateSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
StaggeredStokesPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched);
    return;
} // copyFromPETScVec

void
StaggeredStokesPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                              Vec& petsc_b,
                                              SAMRAIVectorReal<NDIM, double>& x,
                                              SAMRAIVectorReal<NDIM, double>& b)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);

    if (d_initial_guess_nonzero)
    {
        copyToPETScVec(petsc_x, x);
    }
    const bool level_zero = (d_level_num == 0);
    const int f_idx = b.getComponentDescriptorIndex(0);
    const int h_idx = b.getComponentDescriptorIndex(1);
    const auto f_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(f_idx);
    const auto h_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(h_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(d_p_dof_index_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
        Pointer<SideData<NDIM, double>> f_adj_data = patch->getPatchData(f_adj_idx);
        Pointer<CellData<NDIM, double>> h_adj_data = patch->getPatchData(h_adj_idx);
        f_adj_data->copy(*f_data);
        h_adj_data->copy(*h_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = level_zero ?
                                                             Array<BoundaryBox<NDIM>>() :
                                                             d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                          /* boundary type */ 1);
        const Array<BoundaryBox<NDIM>>& type_2_cf_bdry = level_zero ?
                                                             Array<BoundaryBox<NDIM>>() :
                                                             d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                          /* boundary type */ 2);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;

        // Because the matrix contains the dependence of TRACTION conditions on the normal velocity, the boundary
        // condition objects use zero velocities, except in coarse-fine ghost cells, where they use the solution vector.
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, d_u_bc_target_idx, p_idx, d_homogeneous_bc);
        if (at_physical_bdry)
        {
            Pointer<SideData<NDIM, double>> u_bc_target_data = patch->getPatchData(d_u_bc_target_idx);
            if (at_cf_bdry)
            {
                Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(d_u_dof_index_idx);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(u_bc_target_data->getGhostBox()); b; b++)
                    {
                        const SideIndex<NDIM> is(b(), axis, SideIndex<NDIM>::Lower);
                        if ((*u_dof_index_data)(is) < 0)
                        {
                            (*u_bc_target_data)(is) = (*u_data)(is);
                        }
                    }
                }
            }
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *f_adj_data, patch, d_U_problem_coefs, d_U_bc_coefs, d_solution_time, d_homogeneous_bc);
            adjust_rhs_for_normal_stress(
                *f_adj_data, *patch, d_U_bc_coefs, d_U_problem_coefs.getDConstant(), d_solution_time, d_homogeneous_bc);
            if (at_cf_bdry)
            {
                u_bc_target_data->fillAll(0.0);
            }
            adjust_rhs_for_ghost_pressure(
                *f_adj_data, *patch, d_U_bc_coefs, d_P_bc_coef, d_solution_time, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                f_adj_idx, h_adj_idx, d_U_bc_coefs, d_solution_time, d_homogeneous_bc, d_level_num, d_level_num);
        }
        if (at_cf_bdry)
        {
            // enforceNormalVelocityBoundaryConditions() resets the targets of the boundary condition objects, and the
            // adjustment evaluates them.
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_U_bc_coefs, d_P_bc_coef, u_idx, p_idx, d_homogeneous_bc);
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(*f_adj_data,
                                                            *u_data,
                                                            patch,
                                                            d_U_problem_coefs,
                                                            type_1_cf_bdry,
                                                            type_2_cf_bdry,
                                                            d_U_bc_coefs,
                                                            d_solution_time,
                                                            d_homogeneous_bc);
            adjust_rhs_for_pressure_at_cf_boundary(*f_adj_data, *p_data, *p_dof_index_data, patch, type_1_cf_bdry);
        }
    }
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);

    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_b, f_adj_idx, d_u_dof_index_idx, h_adj_idx, d_p_dof_index_idx, d_level);

    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
