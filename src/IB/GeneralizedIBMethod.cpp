// Filename: GeneralizedIBMethod.cpp
// Created on 12 Dec 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <math.h>
#include <stddef.h>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "CellVariable.h"
#include "CoarsenSchedule.h"
#include "Geometry.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "boost/multi_array.hpp"
#include "ibamr/GeneralizedIBMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBKirchhoffRodForceGen.h"
#include "ibamr/IBMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of GeneralizedIBMethod restart file data.
static const int GENERALIZED_IB_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralizedIBMethod::GeneralizedIBMethod(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         bool register_for_restart)
    : IBMethod(object_name, input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object
    // name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Indicate all Lagrangian data needs ghost values to be refilled, and that
    // all intermediate data needs to be initialized.
    d_N_current_needs_ghost_fill = true;
    d_N_new_needs_ghost_fill = true;
    return;
} // GeneralizedIBMethod

GeneralizedIBMethod::~GeneralizedIBMethod()
{
    // intentionally blank
    //
    // NOTE: Parent class constructor unregisters class with the restart
    // manager.
    return;
} // ~GeneralizedIBMethod

void
GeneralizedIBMethod::registerIBKirchhoffRodForceGen(Pointer<IBKirchhoffRodForceGen> ib_force_and_torque_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ib_force_and_torque_fcn);
#endif
    d_ib_force_and_torque_fcn = ib_force_and_torque_fcn;
    return;
} // registerIBKirchhoffRodForceGen

void
GeneralizedIBMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();

    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth();
    const IntVector<NDIM> ghosts = 1;
    const IntVector<NDIM> no_ghosts = 0;

    Pointer<Variable<NDIM> > u_var = d_ib_solver->getVelocityVariable();
    Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
    Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
    if (u_cc_var)
    {
        d_f_var = new CellVariable<NDIM, double>(d_object_name + "::f", NDIM);
        d_w_var = new CellVariable<NDIM, double>(d_object_name + "::w", NDIM);
        d_n_var = new CellVariable<NDIM, double>(d_object_name + "::n", NDIM);
    }
    else if (u_sc_var)
    {
        d_f_var = new SideVariable<NDIM, double>(d_object_name + "::f");
        d_w_var = new SideVariable<NDIM, double>(d_object_name + "::w");
        d_n_var = new SideVariable<NDIM, double>(d_object_name + "::n");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::registerEulerianVariables():\n"
                                 << "  unsupported velocity data centering"
                                 << std::endl);
    }
    registerVariable(d_f_idx, d_f_var, no_ghosts, d_ib_solver->getScratchContext());
    registerVariable(d_w_idx, d_w_var, ib_ghosts, d_ib_solver->getScratchContext());
    registerVariable(d_n_idx, d_n_var, ib_ghosts, d_ib_solver->getScratchContext());
    return;
} // registerEulerianVariables

void
GeneralizedIBMethod::registerEulerianCommunicationAlgorithms()
{
    IBMethod::registerEulerianCommunicationAlgorithms();

    Pointer<Geometry<NDIM> > grid_geom = d_ib_solver->getPatchHierarchy()->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg;
    Pointer<RefineOperator<NDIM> > refine_op;

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_w_idx, d_w_idx, d_w_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name + "::w", refine_alg);

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_n_idx, d_n_idx, d_n_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name + "::n", refine_alg);
    return;
} // registerEulerianCommunicationAlgorithms

void
GeneralizedIBMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_ib_force_and_torque_fcn_needs_init = d_ib_force_fcn_needs_init || d_ib_force_and_torque_fcn_needs_init;
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double start_time = d_ib_solver->getStartTime();

    if (d_ib_force_and_torque_fcn)
    {
        if (d_ib_force_and_torque_fcn_needs_init)
        {
            const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
            resetLagrangianForceAndTorqueFunction(current_time, initial_time);
            d_ib_force_and_torque_fcn_needs_init = false;
        }
    }

    // Look-up or allocate Lagangian data.
    d_D_current_data.resize(finest_ln + 1);
    d_D_new_data.resize(finest_ln + 1);
    d_N_current_data.resize(finest_ln + 1);
    d_N_new_data.resize(finest_ln + 1);
    d_W_current_data.resize(finest_ln + 1);
    d_W_new_data.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_D_current_data[ln] = d_l_data_manager->getLData("D", ln);
        d_D_new_data[ln] = d_l_data_manager->createLData("D_new", ln, NDIM * NDIM);
        d_N_current_data[ln] = d_l_data_manager->createLData("N", ln, NDIM);
        d_N_new_data[ln] = d_l_data_manager->createLData("N_new", ln, NDIM);
        d_W_current_data[ln] = d_l_data_manager->getLData("W", ln);
        d_W_new_data[ln] = d_l_data_manager->createLData("W_new", ln, NDIM);

        // Initialize D^{n+1} to equal D^{n}, and initialize W^{n+1} to equal
        // W^{n}.
        int ierr;
        ierr = VecCopy(d_D_current_data[ln]->getVec(), d_D_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_W_current_data[ln]->getVec(), d_W_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }
    return;
} // preprocessIntegrateData

void
GeneralizedIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr;
        ierr = VecSwap(d_D_current_data[ln]->getVec(), d_D_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_W_current_data[ln]->getVec(), d_W_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }

    // Deallocate Lagrangian scratch data.
    d_D_current_data.clear();
    d_D_new_data.clear();
    d_N_current_data.clear();
    d_N_new_data.clear();
    d_W_current_data.clear();
    d_W_new_data.clear();
    return;
} // postprocessIntegrateData

void
GeneralizedIBMethod::interpolateVelocity(const int u_data_idx,
                                         const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                         const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                         const double data_time)
{
    // Interpolate the linear velocities.
    IBMethod::interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);

    // Interpolate the angular velocities.
    std::vector<Pointer<LData> >* W_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        W_data = &d_W_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        TBOX_ERROR(d_object_name << "::interpolateVelocity():\n"
                                 << "  time-stepping type MIDPOINT_RULE not supported by "
                                    "class GeneralizedIBMethod;\n"
                                 << "  use TRAPEZOIDAL_RULE instead.\n");
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        W_data = &d_W_new_data;
    }

    Pointer<Variable<NDIM> > u_var = d_ib_solver->getVelocityVariable();
    Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
    Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
    if (u_cc_var)
    {
        Pointer<CellVariable<NDIM, double> > w_cc_var = d_w_var;
        getHierarchyMathOps()->curl(d_w_idx, w_cc_var, u_data_idx, u_cc_var, NULL, data_time);
    }
    else if (u_sc_var)
    {
        Pointer<SideVariable<NDIM, double> > w_sc_var = d_w_var;
        getHierarchyMathOps()->curl(d_w_idx, w_sc_var, u_data_idx, u_sc_var, NULL, data_time);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::interpolateVelocity():\n"
                                 << "  unsupported velocity data centering"
                                 << std::endl);
    }
    std::vector<Pointer<LData> >* X_LE_data;
    bool* X_LE_needs_ghost_fill;
    getLECouplingPositionData(&X_LE_data, &X_LE_needs_ghost_fill, data_time);
    getVelocityHierarchyDataOps()->scale(d_w_idx, 0.5, d_w_idx);
    d_l_data_manager->interp(d_w_idx,
                             *W_data,
                             *X_LE_data,
                             std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                             getGhostfillRefineSchedules(d_object_name + "::w"),
                             data_time);
    resetAnchorPointValues(*W_data,
                           /*coarsest_ln*/ 0,
                           /*finest_ln*/ d_hierarchy->getFinestLevelNumber());
    return;
} // interpolateVelocity

void
GeneralizedIBMethod::eulerStep(const double current_time, const double new_time)
{
    IBMethod::eulerStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Update the value of D^{n+1} using forward Euler.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        boost::multi_array_ref<double, 2>& D_current_data = *d_D_current_data[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& W_current_data = *d_W_current_data[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& D_new_data = *d_D_new_data[ln]->getLocalFormVecArray();
        const int n_local = d_l_data_manager->getNumberOfLocalNodes(ln);
        Matrix3d R;
        Vector3d e;
        for (int l = 0; l < n_local; ++l)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                e(d) = W_current_data[l][d];
            }
            const double norm_e = e.norm();
            if (norm_e > std::numeric_limits<double>::epsilon())
            {
                const double theta = norm_e * dt;
                e /= norm_e;
                const double c_t = cos(theta);
                const double s_t = sin(theta);
                R << c_t + (1.0 - c_t) * e(0) * e(0), (1.0 - c_t) * e(0) * e(1) - s_t * e(2),
                    (1.0 - c_t) * e(0) * e(2) + s_t * e(1), (1.0 - c_t) * e(1) * e(0) + s_t * e(2),
                    c_t + (1.0 - c_t) * e(1) * e(1), (1.0 - c_t) * e(1) * e(2) - s_t * e(0),
                    (1.0 - c_t) * e(2) * e(0) - s_t * e(1), (1.0 - c_t) * e(2) * e(1) + s_t * e(0),
                    c_t + (1.0 - c_t) * e(2) * e(2);
                for (int alpha = 0; alpha < 3; ++alpha)
                {
                    Eigen::Map<const Vector3d> D_current_alpha(&D_current_data[l][3 * alpha]);
                    Eigen::Map<Vector3d> D_new_alpha(&D_new_data[l][3 * alpha]);
                    D_new_alpha = R * D_current_alpha;
                }
            }
            else
            {
                for (int alpha = 0; alpha < 3; ++alpha)
                {
                    Eigen::Map<const Vector3d> D_current_alpha(&D_current_data[l][3 * alpha]);
                    Eigen::Map<Vector3d> D_new_alpha(&D_new_data[l][3 * alpha]);
                    D_new_alpha = D_current_alpha;
                }
            }
        }
    }
    return;
} // eulerStep

void
GeneralizedIBMethod::midpointStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR(d_object_name << "::midpointStep():\n"
                             << "  time-stepping type MIDPOINT_RULE not supported by class GeneralizedIBMethod;\n"
                             << "  use TRAPEZOIDAL_RULE instead.\n");
    return;
} // midpointStep

void
GeneralizedIBMethod::trapezoidalStep(const double current_time, const double new_time)
{
    IBMethod::trapezoidalStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Update the value of D^{n+1} using the trapezoidal rule.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        boost::multi_array_ref<double, 2>& D_current_data = *d_D_current_data[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& W_current_data = *d_W_current_data[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& D_new_data = *d_D_new_data[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& W_new_data = *d_W_new_data[ln]->getLocalFormVecArray();
        const int n_local = d_l_data_manager->getNumberOfLocalNodes(ln);
        Matrix3d R;
        Vector3d e;
        for (int l = 0; l < n_local; ++l)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                e(d) = 0.5 * (W_current_data[l][d] + W_new_data[l][d]);
            }
            const double norm_e = e.norm();
            if (norm_e > std::numeric_limits<double>::epsilon())
            {
                const double theta = norm_e * dt;
                e /= norm_e;
                const double c_t = cos(theta);
                const double s_t = sin(theta);
                R << c_t + (1.0 - c_t) * e(0) * e(0), (1.0 - c_t) * e(0) * e(1) - s_t * e(2),
                    (1.0 - c_t) * e(0) * e(2) + s_t * e(1), (1.0 - c_t) * e(1) * e(0) + s_t * e(2),
                    c_t + (1.0 - c_t) * e(1) * e(1), (1.0 - c_t) * e(1) * e(2) - s_t * e(0),
                    (1.0 - c_t) * e(2) * e(0) - s_t * e(1), (1.0 - c_t) * e(2) * e(1) + s_t * e(0),
                    c_t + (1.0 - c_t) * e(2) * e(2);
                for (int alpha = 0; alpha < 3; ++alpha)
                {
                    Eigen::Map<const Vector3d> D_current_alpha(&D_current_data[l][3 * alpha]);
                    Eigen::Map<Vector3d> D_new_alpha(&D_new_data[l][3 * alpha]);
                    D_new_alpha = R * D_current_alpha;
                }
            }
            else
            {
                for (int alpha = 0; alpha < 3; ++alpha)
                {
                    Eigen::Map<const Vector3d> D_current_alpha(&D_current_data[l][3 * alpha]);
                    Eigen::Map<Vector3d> D_new_alpha(&D_new_data[l][3 * alpha]);
                    D_new_alpha = D_current_alpha;
                }
            }
        }
    }
    return;
} // trapezoidalStep

void
GeneralizedIBMethod::computeLagrangianForce(const double data_time)
{
    IBMethod::computeLagrangianForce(data_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    std::vector<Pointer<LData> >* F_data = NULL;
    std::vector<Pointer<LData> >* N_data = NULL;
    std::vector<Pointer<LData> >* X_data = NULL;
    std::vector<Pointer<LData> >* D_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        d_F_current_needs_ghost_fill = true;
        d_N_current_needs_ghost_fill = true;
        F_data = &d_F_current_data;
        N_data = &d_N_current_data;
        X_data = &d_X_current_data;
        D_data = &d_D_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        TBOX_ERROR(d_object_name << "::computeLagrangianForce():\n"
                                 << "  time-stepping type MIDPOINT_RULE not supported by "
                                    "class GeneralizedIBMethod;\n"
                                 << "  use TRAPEZOIDAL_RULE instead.\n");
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        d_F_new_needs_ghost_fill = true;
        d_N_new_needs_ghost_fill = true;
        F_data = &d_F_new_data;
        N_data = &d_N_new_data;
        X_data = &d_X_new_data;
        D_data = &d_D_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecSet((*N_data)[ln]->getVec(), 0.0);
        IBTK_CHKERRQ(ierr);
        if (d_ib_force_and_torque_fcn)
        {
            d_ib_force_and_torque_fcn->computeLagrangianForceAndTorque((*F_data)[ln],
                                                                       (*N_data)[ln],
                                                                       (*X_data)[ln],
                                                                       (*D_data)[ln],
                                                                       d_hierarchy,
                                                                       ln,
                                                                       data_time,
                                                                       d_l_data_manager);
        }
    }
    resetAnchorPointValues(*F_data, coarsest_ln, finest_ln);
    resetAnchorPointValues(*N_data, coarsest_ln, finest_ln);
    return;
} // computeLagrangianForce

void
GeneralizedIBMethod::spreadForce(const int f_data_idx,
                                 RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                 const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                                 const double data_time)
{
    IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);

    std::vector<Pointer<LData> >* N_data = NULL;
    bool* N_needs_ghost_fill = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        N_data = &d_N_current_data;
        N_needs_ghost_fill = &d_N_current_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        TBOX_ERROR(d_object_name << "::spreadForce():\n"
                                 << "  time-stepping type MIDPOINT_RULE not supported by "
                                    "class GeneralizedIBMethod;\n"
                                 << "  use TRAPEZOIDAL_RULE instead.\n");
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        N_data = &d_N_new_data;
        N_needs_ghost_fill = &d_N_new_needs_ghost_fill;
    }

    std::vector<Pointer<LData> >* X_LE_data;
    bool* X_LE_needs_ghost_fill;
    getLECouplingPositionData(&X_LE_data, &X_LE_needs_ghost_fill, data_time);
    getVelocityHierarchyDataOps()->setToScalar(d_n_idx, 0.0, false);
    d_l_data_manager->spread(d_n_idx,
                             *N_data,
                             *X_LE_data,
                             f_phys_bdry_op,
                             std::vector<Pointer<RefineSchedule<NDIM> > >(),
                             data_time,
                             *N_needs_ghost_fill,
                             *X_LE_needs_ghost_fill);
    *N_needs_ghost_fill = false;
    *X_LE_needs_ghost_fill = false;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const std::vector<Pointer<RefineSchedule<NDIM> > >& n_ghostfill_scheds =
        getGhostfillRefineSchedules(d_object_name + "::n");
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        n_ghostfill_scheds[ln]->fillData(data_time);
    }
    Pointer<Variable<NDIM> > u_var = d_ib_solver->getVelocityVariable();
    Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
    Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
    if (u_cc_var)
    {
        Pointer<CellVariable<NDIM, double> > f_cc_var = d_f_var;
        Pointer<CellVariable<NDIM, double> > n_cc_var = d_n_var;
        getHierarchyMathOps()->curl(d_f_idx, f_cc_var, d_n_idx, n_cc_var, NULL, data_time);
    }
    else if (u_sc_var)
    {
        Pointer<SideVariable<NDIM, double> > f_sc_var = d_f_var;
        Pointer<SideVariable<NDIM, double> > n_sc_var = d_n_var;
        getHierarchyMathOps()->curl(d_f_idx, f_sc_var, d_n_idx, n_sc_var, NULL, data_time);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::spreadForce():\n"
                                 << "  unsupported velocity data centering"
                                 << std::endl);
    }
    getVelocityHierarchyDataOps()->axpy(f_data_idx, 0.5, d_f_idx, f_data_idx);
    return;
} // spreadForce

void
GeneralizedIBMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                              Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                              int u_data_idx,
                                              const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                              const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                              int integrator_step,
                                              double init_data_time,
                                              bool initial_time)
{
    // Initialize various Lagrangian data objects required by the conventional
    // IB method.
    IBMethod::initializePatchHierarchy(hierarchy,
                                       gridding_alg,
                                       u_data_idx,
                                       u_synch_scheds,
                                       u_ghost_fill_scheds,
                                       integrator_step,
                                       init_data_time,
                                       initial_time);

    // Initialize various Lagrangian data objects required by the gIB method.
    if (initial_time)
    {
        // Lookup the range of hierarchy levels.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();

        // Initialize the interpolated angular velocity field.
        std::vector<Pointer<LData> > W_data(finest_ln + 1);
        std::vector<Pointer<LData> > X_data(finest_ln + 1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
            W_data[ln] = d_l_data_manager->getLData("W", ln);
        }
        Pointer<Variable<NDIM> > u_var = d_ib_solver->getVelocityVariable();
        Pointer<CellVariable<NDIM, double> > u_cc_var = u_var;
        Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
        if (u_cc_var)
        {
            Pointer<CellVariable<NDIM, double> > w_cc_var = d_w_var;
            getHierarchyMathOps()->curl(d_w_idx, w_cc_var, u_data_idx, u_cc_var, NULL, init_data_time);
        }
        else if (u_sc_var)
        {
            Pointer<SideVariable<NDIM, double> > w_sc_var = d_w_var;
            getHierarchyMathOps()->curl(d_w_idx, w_sc_var, u_data_idx, u_sc_var, NULL, init_data_time);
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializePatchHierarchy():\n"
                                     << "  unsupported velocity data centering"
                                     << std::endl);
        }
        getVelocityHierarchyDataOps()->scale(d_w_idx, 0.5, d_w_idx);
        d_l_data_manager->interp(d_w_idx,
                                 W_data,
                                 X_data,
                                 std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                 getGhostfillRefineSchedules(d_object_name + "::w"),
                                 init_data_time);
        resetAnchorPointValues(W_data, coarsest_ln, finest_ln);
    }

    // Indicate that the force-and-torque strategy needs to be re-initialized.
    d_ib_force_and_torque_fcn_needs_init = true;
    return;
} // initializePatchHierarchy

void
GeneralizedIBMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                         int level_number,
                                         double init_data_time,
                                         bool can_be_refined,
                                         bool initial_time,
                                         Pointer<BasePatchLevel<NDIM> > old_level,
                                         bool allocate_data)
{
    IBMethod::initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        // 1. Allocate LData corresponding to the curvilinear mesh node
        //    directors and angular velocities.
        Pointer<LData> D_data = d_l_data_manager->createLData("D",
                                                              level_number,
                                                              NDIM * NDIM,
                                                              /*manage_data*/ true);
        Pointer<LData> W_data = d_l_data_manager->createLData("W", level_number, NDIM, /*manage_data*/ true);

        // 2. Initialize the Lagrangian data.
        static const int global_index_offset = 0;
        static const int local_index_offset = 0;
        d_l_initializer->initializeDirectorDataOnPatchLevel(global_index_offset,
                                                            local_index_offset,
                                                            D_data,
                                                            hierarchy,
                                                            level_number,
                                                            init_data_time,
                                                            can_be_refined,
                                                            initial_time,
                                                            d_l_data_manager);

        // 3. Register data with any registered data writer.
        if (d_silo_writer)
        {
            d_silo_writer->registerVariableData("D1", D_data, 0, 3, level_number);
            d_silo_writer->registerVariableData("D2", D_data, 3, 3, level_number);
            d_silo_writer->registerVariableData("D3", D_data, 6, 3, level_number);
            d_silo_writer->registerVariableData("W", W_data, level_number);
        }
    }
    return;
} // initializeLevelData

void
GeneralizedIBMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);
    db->putInteger("GENERALIZED_IB_METHOD_VERSION", GENERALIZED_IB_METHOD_VERSION);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
GeneralizedIBMethod::resetLagrangianForceAndTorqueFunction(const double init_data_time, const bool initial_time)
{
    if (!d_ib_force_and_torque_fcn) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_ib_force_and_torque_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
} // resetLagrangianForceAndTorqueFunction

void
GeneralizedIBMethod::getFromInput(Pointer<Database> /*db*/, bool /*is_from_restart*/)
{
    // intentionally blank
    return;
} // getFromInput

void
GeneralizedIBMethod::getFromRestart()
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
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("GENERALIZED_IB_METHOD_VERSION");
    if (ver != GENERALIZED_IB_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
