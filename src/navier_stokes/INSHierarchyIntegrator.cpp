// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/ibtk_utilities.h"

#include "FaceVariable.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of INSHierarchyIntegrator restart file data.
static const int INS_HIERARCHY_INTEGRATOR_VERSION = 4;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSHierarchyIntegrator::~INSHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U_star_bc_coefs[d];
        d_U_star_bc_coefs[d] = nullptr;
    }
    return;
} // ~INSHierarchyIntegrator

TimeSteppingType
INSHierarchyIntegrator::getViscousTimeSteppingType() const
{
    return d_viscous_time_stepping_type;
} // getViscousTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getConvectiveTimeSteppingType() const
{
    return d_convective_time_stepping_type;
} // getConvectiveTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getInitialConvectiveTimeSteppingType() const
{
    return d_init_convective_time_stepping_type;
} // getInitialConvectiveTimeSteppingType

void
INSHierarchyIntegrator::registerAdvDiffHierarchyIntegrator(Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(adv_diff_hier_integrator);
    // Bad things happen if the same integrator is registered twice.
    auto pointer_compare = [adv_diff_hier_integrator](Pointer<AdvDiffHierarchyIntegrator> integrator) -> bool
    { return adv_diff_hier_integrator.getPointer() == integrator.getPointer(); };
    TBOX_ASSERT(std::find_if(d_adv_diff_hier_integrators.begin(), d_adv_diff_hier_integrators.end(), pointer_compare) ==
                d_adv_diff_hier_integrators.end());
#endif
    d_adv_diff_hier_integrators.push_back(adv_diff_hier_integrator);
    registerChildHierarchyIntegrator(adv_diff_hier_integrator);
    adv_diff_hier_integrator->registerAdvectionVelocity(d_U_adv_diff_var);
    adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, !d_Q_fcn);
    return;
} // registerAdvDiffHierarchyIntegrator

void
INSHierarchyIntegrator::setStokesSpecifications(StokesSpecifications problem_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_problem_coefs = problem_coefs;
    return;
} // setStokesSpecifications

const StokesSpecifications*
INSHierarchyIntegrator::getStokesSpecifications() const
{
    return &d_problem_coefs;
} // getStokesSpecifications

void
INSHierarchyIntegrator::registerPhysicalBoundaryConditions(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // registerPhysicalBoundaryConditions

const std::vector<RobinBcCoefStrategy<NDIM>*>&
INSHierarchyIntegrator::getVelocityBoundaryConditions() const
{
    return d_U_bc_coefs;
} // getVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSHierarchyIntegrator::getPressureBoundaryConditions() const
{
    return d_P_bc_coef;
} // getPressureBoundaryConditions

void
INSHierarchyIntegrator::registerVelocityInitialConditions(Pointer<CartGridFunction> U_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_U_init = U_init;
    return;
} // registerVelocityInitialConditions

void
INSHierarchyIntegrator::registerPressureInitialConditions(Pointer<CartGridFunction> P_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_P_init = P_init;
    return;
} // registerPressureInitialConditions

void
INSHierarchyIntegrator::registerBodyForceFunction(Pointer<CartGridFunction> F_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_F_fcn)
    {
        Pointer<CartGridFunctionSet> p_F_fcn = d_F_fcn;
        if (!p_F_fcn)
        {
            pout << d_object_name << "::registerBodyForceFunction(): WARNING:\n"
                 << "  body force function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the body force term value.\n";
            p_F_fcn = new CartGridFunctionSet(d_object_name + "::body_force_function_set");
            p_F_fcn->addFunction(d_F_fcn);
        }
        p_F_fcn->addFunction(F_fcn);
        d_F_fcn = p_F_fcn;
    }
    else
    {
        d_F_fcn = F_fcn;
    }
    return;
} // registerBodyForceFunction

void
INSHierarchyIntegrator::registerFluidSourceFunction(Pointer<CartGridFunction> Q_fcn)
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2(
        "INSHierarchyIntegrator", "registerFluidSourceFunction", "registerVelocityDivergenceFunction");
    registerVelocityDivergenceFunction(Q_fcn);
    return;
} // registerFluidSourceFunction

void
INSHierarchyIntegrator::registerVelocityDivergenceFunction(Pointer<CartGridFunction> Q_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_Q_fcn)
    {
        Pointer<CartGridFunctionSet> p_Q_fcn = d_Q_fcn;
        if (!p_Q_fcn)
        {
            pout << d_object_name << "::registerVelocityDivergenceFunction(): WARNING:\n"
                 << "  velocity divergence function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were "
                    "registered with the solver\n"
                 << "  when evaluating the fluid source term value.\n";
            p_Q_fcn = new CartGridFunctionSet(d_object_name + "::fluid_source_function_set");
            p_Q_fcn->addFunction(d_Q_fcn);
        }
        p_Q_fcn->addFunction(Q_fcn);
    }
    else
    {
        d_Q_fcn = Q_fcn;
    }
    return;
} // registerVelocityDivergenceFunction

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getVelocityVariable() const
{
    return d_U_var;
} // getVelocityVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
} // getPressureVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getBodyForceVariable() const
{
    return d_F_var;
} // getBodyForceVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getFluidSourceVariable() const
{
    IBTK_DEPRECATED_MEMBER_FUNCTION2(
        "INSHierarchyIntegrator", "getFluidSourceVariable", "getVelocityDivergenceVariable");
    return d_Q_var;
} // getFluidSourceVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getVelocityDivergenceVariable() const
{
    return d_Q_var;
} // getVelocityDivergenceVariable

Pointer<FaceVariable<NDIM, double> >
INSHierarchyIntegrator::getAdvectionVelocityVariable() const
{
    return d_U_adv_diff_var;
} // getAdvectionVelocityVariable

std::vector<RobinBcCoefStrategy<NDIM>*>
INSHierarchyIntegrator::getIntermediateVelocityBoundaryConditions() const
{
    return d_U_star_bc_coefs;
} // getIntermediateVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSHierarchyIntegrator::getProjectionBoundaryConditions() const
{
    return d_Phi_bc_coef.get();
} // getProjectionBoundaryConditions

void
INSHierarchyIntegrator::setCreepingFlow(bool creeping_flow)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_creeping_flow = creeping_flow;
    d_convective_op.setNull();
    d_convective_difference_form = UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE;
    return;
} // setCreepingFlow

bool
INSHierarchyIntegrator::getCreepingFlow() const
{
    return d_creeping_flow;
} // getCreepingFlow

void
INSHierarchyIntegrator::setConvectiveOperatorType(const std::string& op_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_convective_op_type = op_type;
    return;
} // setConvectiveOperatorType

const std::string&
INSHierarchyIntegrator::getConvectiveOperatorType() const
{
    return d_convective_op_type;
} // getConvectiveOperatorType

void
INSHierarchyIntegrator::setConvectiveDifferencingType(ConvectiveDifferencingType difference_form)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_convective_difference_form = difference_form;
    return;
} // setConvectiveDifferencingType

ConvectiveDifferencingType
INSHierarchyIntegrator::getConvectiveDifferencingType() const
{
    return d_convective_difference_form;
} // getConvectiveDifferencingType

void
INSHierarchyIntegrator::setConvectiveOperator(Pointer<ConvectiveOperator> convective_op)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
#endif
    d_convective_op = convective_op;
    d_creeping_flow = !d_convective_op;
    return;
} // setConvectiveOperator

void
INSHierarchyIntegrator::setConvectiveOperatorNeedsInit()
{
    d_convective_op_needs_init = true;
    return;
}

void
INSHierarchyIntegrator::setVelocitySubdomainSolver(Pointer<PoissonSolver> velocity_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_velocity_solver);
#endif
    d_velocity_solver = velocity_solver;
    return;
} // setVelocitySubdomainSolver

void
INSHierarchyIntegrator::setVelocitySubdomainSolverNeedsInit()
{
    d_velocity_solver_needs_init = true;
    return;
}

void
INSHierarchyIntegrator::setPressureSubdomainSolver(Pointer<PoissonSolver> pressure_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_pressure_solver);
#endif
    d_pressure_solver = pressure_solver;
    return;
} // setPressureSubdomainSolver

void
INSHierarchyIntegrator::setPressureSubdomainSolverNeedsInit()
{
    d_pressure_solver_needs_init = true;
    return;
}

int
INSHierarchyIntegrator::getNumberOfCycles() const
{
    int num_cycles = d_num_cycles;
    if (!d_creeping_flow && IBTK::abs_equal_eps(d_integrator_time, d_start_time) &&
        is_multistep_time_stepping_type(d_convective_time_stepping_type) &&
        d_init_convective_time_stepping_type != FORWARD_EULER)
    {
        num_cycles = std::max(2, num_cycles);
    }
    return num_cycles;
} // getNumberOfCycles

double
INSHierarchyIntegrator::getCurrentCFLNumber() const
{
    return d_cfl_current;
} // getCurrentCFLNumber

void
INSHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                      const double new_time,
                                                      const bool skip_synchronize_new_state_data,
                                                      const int num_cycles)
{
    // The child class has the data indices so we have to look them up manually at this point
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    updateCurrentCFLNumber(var_db->mapVariableAndContextToIndex(getVelocityVariable(), getNewContext()),
                           new_time - current_time);

    if (!d_parent_integrator && d_enable_logging)
    {
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << d_cfl_current << "\n";
    }

    HierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator::INSHierarchyIntegrator(std::string object_name,
                                               Pointer<Database> input_db,
                                               Pointer<Variable<NDIM> > U_var,
                                               Pointer<Variable<NDIM> > P_var,
                                               Pointer<Variable<NDIM> > F_var,
                                               Pointer<Variable<NDIM> > Q_var,
                                               bool register_for_restart)
    : INSHierarchyIntegrator(std::move(object_name),
                             input_db,
                             U_var,
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             P_var,
                             "CONSERVATIVE_COARSEN",
                             "LINEAR_REFINE",
                             F_var,
                             "CONSERVATIVE_COARSEN",
                             "CONSERVATIVE_LINEAR_REFINE",
                             Q_var,
                             "CONSERVATIVE_COARSEN",
                             "CONSTANT_REFINE",
                             register_for_restart)
{
}

INSHierarchyIntegrator::INSHierarchyIntegrator(std::string object_name,
                                               Pointer<Database> input_db,
                                               Pointer<Variable<NDIM> > U_var,
                                               std::string U_default_coarsen_type,
                                               std::string U_default_refine_type,
                                               Pointer<Variable<NDIM> > P_var,
                                               std::string P_default_coarsen_type,
                                               std::string P_default_refine_type,
                                               Pointer<Variable<NDIM> > F_var,
                                               std::string F_default_coarsen_type,
                                               std::string F_default_refine_type,
                                               Pointer<Variable<NDIM> > Q_var,
                                               std::string Q_default_coarsen_type,
                                               std::string Q_default_refine_type,
                                               bool register_for_restart)
    : HierarchyIntegrator(std::move(object_name), input_db, register_for_restart),
      d_U_var(U_var),
      d_U_coarsen_type(std::move(U_default_coarsen_type)),
      d_U_refine_type(std::move(U_default_refine_type)),
      d_P_var(P_var),
      d_P_coarsen_type(std::move(P_default_coarsen_type)),
      d_P_refine_type(std::move(P_default_refine_type)),
      d_F_var(F_var),
      d_F_coarsen_type(std::move(F_default_coarsen_type)),
      d_F_refine_type(std::move(F_default_refine_type)),
      d_Q_var(Q_var),
      d_Q_coarsen_type(std::move(Q_default_coarsen_type)),
      d_Q_refine_type(std::move(Q_default_refine_type)),
      d_default_bc_coefs(d_object_name + "::default_bc_coefs", Pointer<Database>(nullptr)),
      d_bc_coefs(NDIM, nullptr)
{
    // Set some default values.
    d_convective_op_input_db = new MemoryDatabase(d_object_name + "::convective_op_input_db");

    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet (solid-wall) boundary conditions for the velocity.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coefs.setBoundaryValue(2 * d, 0.0);
        d_default_bc_coefs.setBoundaryValue(2 * d + 1, 0.0);
    }
    registerPhysicalBoundaryConditions(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, &d_default_bc_coefs));

    // Setup physical boundary conditions objects.
    d_U_star_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_star_bc_coefs[d] = new INSIntermediateVelocityBcCoef(d, d_bc_coefs);
    }
    d_Phi_bc_coef = std::make_unique<INSProjectionBcCoef>(d_bc_coefs);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Initialize an advection velocity variable.  NOTE: Patch data are
    // allocated for this variable only when an advection-diffusion solver is
    // registered with the INSHierarchyIntegrator.
    d_U_adv_diff_var = new FaceVariable<NDIM, double>(d_object_name + "::U_adv_diff");
    return;
} // INSHierarchyIntegrator

void
INSHierarchyIntegrator::updateCurrentCFLNumber(const int data_idx, const double dt)
{
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM, double> patch_cc_ops;
    PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx, dx + NDIM));
            Pointer<CellData<NDIM, double> > u_cc_new_data = patch->getPatchData(data_idx);
            Pointer<SideData<NDIM, double> > u_sc_new_data = patch->getPatchData(data_idx);
#ifndef NDEBUG
            TBOX_ASSERT(u_cc_new_data || u_sc_new_data);
#endif
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }

    d_cfl_current = IBTK_MPI::maxReduction(cfl_max);
} // updateCurrentCFLNumber

double
INSHierarchyIntegrator::getMaximumVorticityMagnitude(const int Omega_idx)
{
    TBOX_ASSERT(d_hier_math_ops);
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    double max_vorticity_norm = 0.0;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > Omega_data_ptr = patch->getPatchData(Omega_idx);
            Pointer<CellData<NDIM, double> > cc_wgt_data_ptr = patch->getPatchData(wgt_cc_idx);
            TBOX_ASSERT(Omega_data_ptr);
            TBOX_ASSERT(cc_wgt_data_ptr);
            const CellData<NDIM, double>& Omega_data = *Omega_data_ptr;
            const CellData<NDIM, double>& cc_wgt_data = *cc_wgt_data_ptr;
            for (CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const hier::Index<NDIM>& i = ic();
                if (cc_wgt_data(i) > 0.0)
                {
                    if (NDIM == 2)
                    {
                        max_vorticity_norm = std::max(max_vorticity_norm, std::abs(Omega_data(i)));
                    }
                    else
                    {
                        double norm_Omega_sq = 0.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const double o = Omega_data(i, d);
                            norm_Omega_sq += o * o;
                        }
                        max_vorticity_norm = std::max(max_vorticity_norm, std::sqrt(norm_Omega_sq));
                    }
                }
            }
        }
    }

    return IBTK_MPI::maxReduction(max_vorticity_norm);
}

void
INSHierarchyIntegrator::tagCellsByVorticityMagnitude(const int level_number, const int Omega_idx, const int tag_idx)
{
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
    const double Omega_max = getMaximumVorticityMagnitude(Omega_idx);

    // Tag cells based on the magnitude of the vorticity.
    //
    // Note that if either the relative or absolute threshold is zero for a
    // particular level, no tagging is performed on that level.
    double Omega_rel_thresh = 0.0;
    if (d_Omega_rel_thresh.size() > 0)
    {
        Omega_rel_thresh = d_Omega_rel_thresh[std::max(std::min(level_number, d_Omega_rel_thresh.size() - 1), 0)];
    }
    double Omega_abs_thresh = 0.0;
    if (d_Omega_abs_thresh.size() > 0)
    {
        Omega_abs_thresh = d_Omega_abs_thresh[std::max(std::min(level_number, d_Omega_abs_thresh.size() - 1), 0)];
    }
    if (Omega_rel_thresh > 0.0 || Omega_abs_thresh > 0.0)
    {
        double thresh = std::numeric_limits<double>::max();
        if (Omega_rel_thresh > 0.0) thresh = std::min(thresh, Omega_rel_thresh * Omega_max);
        if (Omega_abs_thresh > 0.0) thresh = std::min(thresh, Omega_abs_thresh);
        thresh += std::sqrt(std::numeric_limits<double>::epsilon());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > Omega_data_ptr = patch->getPatchData(Omega_idx);
            Pointer<CellData<NDIM, int> > tag_data_ptr = patch->getPatchData(tag_idx);
            TBOX_ASSERT(Omega_data_ptr);
            TBOX_ASSERT(tag_data_ptr);
            const CellData<NDIM, double>& Omega_data = *Omega_data_ptr;
            CellData<NDIM, int>& tag_data = *tag_data_ptr;
            for (CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const hier::Index<NDIM>& i = ic();
                double norm_Omega = 0.0;
                if (NDIM == 2)
                {
                    norm_Omega = std::abs(Omega_data(i));
                }
                else
                {
                    double norm_Omega_sq = 0.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const double o = Omega_data(i, d);
                        norm_Omega_sq += o * o;
                    }
                    norm_Omega = std::sqrt(norm_Omega_sq);
                }
                if (norm_Omega > thresh)
                {
                    tag_data(i) = 1;
                }
            }
        }
    }
    return;
}

double
INSHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = HierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_cfl_max * getStableTimestep(level));
    }
    return dt;
} // getMaximumTimeStepSizeSpecialized

double
INSHierarchyIntegrator::getStableTimestep(Pointer<PatchLevel<NDIM> > level) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt, getStableTimestep(patch));
    }
    stable_dt = IBTK_MPI::minReduction(stable_dt);
    return stable_dt;
} // getStableTimestep

void
INSHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("INS_HIERARCHY_INTEGRATOR_VERSION", INS_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_time_stepping_type", enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type));
    db->putString("d_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_convective_time_stepping_type));
    db->putString("d_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_init_convective_time_stepping_type));
    db->putDouble("d_rho", d_problem_coefs.getRho());
    db->putDouble("d_mu", d_problem_coefs.getMu());
    db->putDouble("d_lambda", d_problem_coefs.getLambda());
    db->putDouble("d_cfl_current", d_cfl_current);
    db->putDouble("d_cfl_max", d_cfl_max);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    if (d_Omega_rel_thresh.size() > 0) db->putDoubleArray("d_Omega_rel_thresh", d_Omega_rel_thresh);
    if (d_Omega_abs_thresh.size() > 0) db->putDoubleArray("d_Omega_abs_thresh", d_Omega_abs_thresh);
    db->putBool("d_normalize_pressure", d_normalize_pressure);
    db->putBool("d_normalize_velocity", d_normalize_velocity);
    db->putString("d_convective_op_type", d_convective_op_type);
    db->putString("d_convective_difference_form",
                  enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form));
    db->putBool("d_creeping_flow", d_creeping_flow);
    db->putDouble("d_regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);
    db->putDouble("d_U_scale", d_U_scale);
    db->putDouble("d_P_scale", d_P_scale);
    db->putDouble("d_F_scale", d_F_scale);
    db->putDouble("d_Q_scale", d_Q_scale);
    db->putDouble("d_Omega_scale", d_Omega_scale);
    db->putDouble("d_Div_U_scale", d_Div_U_scale);
    db->putDouble("d_EE_scale", d_EE_scale);
    db->putBool("d_output_U", d_output_U);
    db->putBool("d_output_P", d_output_P);
    db->putBool("d_output_F", d_output_F);
    db->putBool("d_output_Q", d_output_Q);
    db->putBool("d_output_Omega", d_output_Omega);
    db->putBool("d_output_Div_U", d_output_Div_U);
    db->putBool("d_output_EE", d_output_EE);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSHierarchyIntegrator::getFromInput(Pointer<Database> db, const bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->keyExists("viscous_time_stepping_type"))
            d_viscous_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("viscous_time_stepping_type"));
        else if (db->keyExists("viscous_timestepping_type"))
            d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("viscous_timestepping_type"));
        if (db->keyExists("convective_time_stepping_type"))
            d_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_time_stepping_type"));
        else if (db->keyExists("convective_timestepping_type"))
            d_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_timestepping_type"));
        if (db->keyExists("init_convective_time_stepping_type"))
            d_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_time_stepping_type"));
        else if (db->keyExists("init_convective_timestepping_type"))
            d_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_timestepping_type"));
        if (db->keyExists("rho"))
        {
            d_problem_coefs.setRho(db->getDouble("rho"));
        }
        else
        {
            TBOX_WARNING("INSHierarchyIntegrator::getFromInput()\n"
                         << "  no constant density specified;\n"
                         << "  setting to quiet_NaN");
            d_problem_coefs.setRho(std::numeric_limits<double>::quiet_NaN());
        }

        if (db->keyExists("mu"))
        {
            d_problem_coefs.setMu(db->getDouble("mu"));
        }
        else
        {
            TBOX_WARNING("INSHierarchyIntegrator::getFromInput()\n"
                         << "  no constant viscosity specified;\n"
                         << "  setting to quiet_NaN");
            d_problem_coefs.setMu(std::numeric_limits<double>::quiet_NaN());
        }

        if (db->keyExists("lambda"))
        {
            d_problem_coefs.setLambda(db->getDouble("lambda"));
        }
        else
        {
            d_problem_coefs.setLambda(0.0);
        }
    }
    if (db->keyExists("num_cycles")) d_num_cycles = db->getInteger("num_cycles");
    if (db->keyExists("cfl"))
        d_cfl_max = db->getDouble("cfl");
    else if (db->keyExists("cfl_max"))
        d_cfl_max = db->getDouble("cfl_max");
    else if (db->keyExists("CFL"))
        d_cfl_max = db->getDouble("CFL");
    else if (db->keyExists("CFL_max"))
        d_cfl_max = db->getDouble("CFL_max");
    if (db->keyExists("using_vorticity_tagging")) d_using_vorticity_tagging = db->getBool("using_vorticity_tagging");
    if (db->keyExists("Omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("Omega_rel_thresh");
    else if (db->keyExists("omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("omega_rel_thresh");
    else if (db->keyExists("vorticity_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
    if (db->keyExists("Omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("Omega_abs_thresh");
    else if (db->keyExists("omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("omega_abs_thresh");
    else if (db->keyExists("vorticity_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
    if (db->keyExists("normalize_pressure")) d_normalize_pressure = db->getBool("normalize_pressure");
    if (db->keyExists("normalize_velocity")) d_normalize_velocity = db->getBool("normalize_velocity");
    if (db->keyExists("convective_op_type"))
        d_convective_op_type = db->getString("convective_op_type");
    else if (db->keyExists("convective_operator_type"))
        d_convective_op_type = db->getString("convective_operator_type");
    else if (db->keyExists("default_convective_op_type"))
        d_convective_op_type = db->getString("default_convective_op_type");
    else if (db->keyExists("default_convective_operator_type"))
        d_convective_op_type = db->getString("default_convective_operator_type");
    if (db->keyExists("convective_difference_form"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_form"));
    else if (db->keyExists("convective_difference_type"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_type"));
    else if (db->keyExists("default_convective_difference_form"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_form"));
    else if (db->keyExists("default_convective_difference_type"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_type"));
    if (db->keyExists("convective_op_db"))
        d_convective_op_input_db = db->getDatabase("convective_op_db");
    else if (db->keyExists("default_convective_op_db"))
        d_convective_op_input_db = db->getDatabase("default_convective_op_db");
    if (db->keyExists("creeping_flow")) d_creeping_flow = db->getBool("creeping_flow");
    if (db->keyExists("regrid_max_div_growth_factor"))
        d_regrid_max_div_growth_factor = db->getDouble("regrid_max_div_growth_factor");
    if (db->keyExists("U_scale")) d_U_scale = db->getDouble("U_scale");
    if (db->keyExists("P_scale")) d_P_scale = db->getDouble("P_scale");
    if (db->keyExists("F_scale")) d_F_scale = db->getDouble("F_scale");
    if (db->keyExists("Q_scale")) d_Q_scale = db->getDouble("Q_scale");
    if (db->keyExists("Omega_scale")) d_Omega_scale = db->getDouble("Omega_scale");
    if (db->keyExists("Div_U_scale")) d_Div_U_scale = db->getDouble("Div_U_scale");
    if (db->keyExists("EE_scale")) d_EE_scale = db->getDouble("EE_scale");
    if (db->keyExists("output_U")) d_output_U = db->getBool("output_U");
    if (db->keyExists("output_P")) d_output_P = db->getBool("output_P");
    if (db->keyExists("output_F")) d_output_F = db->getBool("output_F");
    if (db->keyExists("output_Q")) d_output_Q = db->getBool("output_Q");
    if (db->keyExists("output_Omega")) d_output_Omega = db->getBool("output_Omega");
    if (db->keyExists("output_Div_U")) d_output_Div_U = db->getBool("output_Div_U");
    if (db->keyExists("output_EE")) d_output_EE = db->getBool("output_EE");
    if (db->keyExists("traction_bc_type"))
        d_traction_bc_type = string_to_enum<TractionBcType>(db->getString("traction_bc_type"));
    if (db->keyExists("use_div_sink_drag_term")) d_use_div_sink_drag_term = db->getBool("use_div_sink_drag_term");

    if (db->keyExists("velocity_solver_type"))
    {
        d_velocity_solver_type = db->getString("velocity_solver_type");
        if (db->keyExists("velocity_solver_db")) d_velocity_solver_db = db->getDatabase("velocity_solver_db");
    }
    if (!d_velocity_solver_db) d_velocity_solver_db = new MemoryDatabase("velocity_solver_db");

    if (db->keyExists("velocity_precond_type"))
    {
        d_velocity_precond_type = db->getString("velocity_precond_type");
        if (db->keyExists("velocity_precond_db")) d_velocity_precond_db = db->getDatabase("velocity_precond_db");
    }
    if (!d_velocity_precond_db) d_velocity_precond_db = new MemoryDatabase("velocity_precond_db");

    if (db->keyExists("velocity_sub_precond_type"))
    {
        d_velocity_sub_precond_type = db->getString("velocity_sub_precond_type");
        if (db->keyExists("velocity_sub_precond_db"))
            d_velocity_sub_precond_db = db->getDatabase("velocity_sub_precond_db");
    }
    if (!d_velocity_sub_precond_db) d_velocity_sub_precond_db = new MemoryDatabase("velocity_sub_precond_db");

    if (db->keyExists("pressure_solver_type"))
    {
        d_pressure_solver_type = db->getString("pressure_solver_type");
        if (db->keyExists("pressure_solver_db")) d_pressure_solver_db = db->getDatabase("pressure_solver_db");
    }
    if (!d_pressure_solver_db) d_pressure_solver_db = new MemoryDatabase("pressure_solver_db");

    if (db->keyExists("pressure_precond_type"))
    {
        d_pressure_precond_type = db->getString("pressure_precond_type");
        if (db->keyExists("pressure_precond_db")) d_pressure_precond_db = db->getDatabase("pressure_precond_db");
    }
    if (!d_pressure_precond_db) d_pressure_precond_db = new MemoryDatabase("pressure_precond_db");

    if (db->keyExists("pressure_sub_precond_type"))
    {
        d_pressure_sub_precond_type = db->getString("pressure_sub_precond_type");
        if (db->keyExists("pressure_sub_precond_db"))
            d_pressure_sub_precond_db = db->getDatabase("pressure_sub_precond_db");
    }
    if (!d_pressure_sub_precond_db) d_pressure_sub_precond_db = new MemoryDatabase("pressure_sub_precond_db");

    if (db->keyExists("regrid_projection_solver_type"))
    {
        d_regrid_projection_solver_type = db->getString("regrid_projection_solver_type");
        if (db->keyExists("regrid_projection_solver_db"))
            d_regrid_projection_solver_db = db->getDatabase("regrid_projection_solver_db");
    }
    if (!d_regrid_projection_solver_db)
        d_regrid_projection_solver_db = new MemoryDatabase("regrid_projection_solver_db");

    if (db->keyExists("regrid_projection_precond_type"))
    {
        d_regrid_projection_precond_type = db->getString("regrid_projection_precond_type");
        if (db->keyExists("regrid_projection_precond_db"))
            d_regrid_projection_precond_db = db->getDatabase("regrid_projection_precond_db");
    }
    if (!d_regrid_projection_precond_db)
        d_regrid_projection_precond_db = new MemoryDatabase("regrid_projection_precond_db");

    if (db->keyExists("regrid_projection_sub_precond_type"))
    {
        d_regrid_projection_sub_precond_type = db->getString("regrid_projection_sub_precond_type");
        if (db->keyExists("regrid_projection_sub_precond_db"))
            d_regrid_projection_sub_precond_db = db->getDatabase("regrid_projection_sub_precond_db");
    }
    if (!d_regrid_projection_sub_precond_db)
        d_regrid_projection_sub_precond_db = new MemoryDatabase("regrid_projection_sub_precond_db");
    if (db->keyExists("U_coarsen_type")) d_U_coarsen_type = db->getString("U_coarsen_type");
    if (db->keyExists("P_coarsen_type")) d_P_coarsen_type = db->getString("P_coarsen_type");
    if (db->keyExists("F_coarsen_type")) d_F_coarsen_type = db->getString("F_coarsen_type");
    if (db->keyExists("Q_coarsen_type")) d_Q_coarsen_type = db->getString("Q_coarsen_type");
    return;
} // getFromInput

void
INSHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("INS_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_viscous_time_stepping_type"));
    d_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_convective_time_stepping_type"));
    d_init_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_init_convective_time_stepping_type"));
    d_problem_coefs.setRho(db->getDouble("d_rho"));
    d_problem_coefs.setMu(db->getDouble("d_mu"));
    d_problem_coefs.setLambda(db->getDouble("d_lambda"));
    d_num_cycles = db->getInteger("d_num_cycles");
    d_cfl_current = db->getDouble("d_cfl_current");
    d_cfl_max = db->getDouble("d_cfl_max");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    if (db->keyExists("d_Omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    else
        d_Omega_rel_thresh.resizeArray(0);
    if (db->keyExists("d_Omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    else
        d_Omega_abs_thresh.resizeArray(0);
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_normalize_velocity = db->getBool("d_normalize_velocity");
    d_convective_op_type = db->getString("d_convective_op_type");
    d_convective_difference_form =
        string_to_enum<ConvectiveDifferencingType>(db->getString("d_convective_difference_form"));
    d_creeping_flow = db->getBool("d_creeping_flow");
    d_regrid_max_div_growth_factor = db->getDouble("d_regrid_max_div_growth_factor");
    d_U_scale = db->getDouble("d_U_scale");
    d_P_scale = db->getDouble("d_P_scale");
    d_F_scale = db->getDouble("d_F_scale");
    d_Q_scale = db->getDouble("d_Q_scale");
    d_Omega_scale = db->getDouble("d_Omega_scale");
    d_Div_U_scale = db->getDouble("d_Div_U_scale");
    d_EE_scale = db->getDouble("d_EE_scale");
    d_output_U = db->getBool("d_output_U");
    d_output_P = db->getBool("d_output_P");
    d_output_F = db->getBool("d_output_F");
    d_output_Q = db->getBool("d_output_Q");
    d_output_Omega = db->getBool("d_output_Omega");
    d_output_Div_U = db->getBool("d_output_Div_U");
    d_output_EE = db->getBool("d_output_EE");
    return;
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
