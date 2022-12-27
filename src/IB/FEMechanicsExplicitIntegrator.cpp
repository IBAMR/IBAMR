// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/FEMechanicsExplicitIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/FEProjector.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/libmesh_utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"

#include "petscvec.h"

#include <algorithm>
#include <memory>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of FEMechanicsExplicitIntegrator restart file data.
const int EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION = 0;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             const Pointer<Database>& input_db,
                                                             MeshBase* mesh,
                                                             bool register_for_restart,
                                                             const std::string& restart_read_dirname,
                                                             unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, mesh, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(input_db);
}

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             const Pointer<Database>& input_db,
                                                             const std::vector<MeshBase*>& meshes,
                                                             bool register_for_restart,
                                                             const std::string& restart_read_dirname,
                                                             unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(input_db);
}

void
FEMechanicsExplicitIntegrator::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    FEMechanicsBase::preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize variables.
    d_X_vecs->copy("solution", { "current", "new", "half" });
    d_U_vecs->copy("solution", { "current", "new", "half" });
    d_F_vecs->copy("solution", { "current", "new", "half" });
    if (d_P_vecs) d_P_vecs->copy("solution", { "current", "new", "half" });
}

void
FEMechanicsExplicitIntegrator::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    std::vector<std::vector<PetscVector<double>*> > vecs{ d_X_vecs->get("new"),
                                                          d_U_vecs->get("new"),
                                                          d_F_vecs->get("new") };
    if (d_P_vecs) vecs.push_back(d_P_vecs->get("new"));
    batch_vec_ghost_update(vecs, INSERT_VALUES, SCATTER_FORWARD);

    d_X_vecs->copy("new", { "solution", "current" });
    d_U_vecs->copy("new", { "solution", "current" });
    d_F_vecs->copy("new", { "solution", "current" });
    if (d_P_vecs) d_P_vecs->copy("new", { "solution", "current" });

    FEMechanicsBase::postprocessIntegrateData(current_time, new_time, num_cycles);
}

void
FEMechanicsExplicitIntegrator::forwardEulerStep(const double current_time, const double new_time)
{
    // U^{n+1} := U^{n} + (dt/rho) F^{n}
    // X^{n+1} := X^{n} + dt       U^{n}
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        doForwardEulerStep(d_X_vecs->get("new", part),
                           d_U_vecs->get("new", part),
                           d_dynamic_pressure_part[part] ? &d_P_vecs->get("new", part) : nullptr,
                           d_X_vecs->get("current", part),
                           d_U_vecs->get("current", part),
                           partHasPressure(part) ? &d_P_vecs->get("current", part) : nullptr,
                           d_F_vecs->get("current", part),
                           d_dynamic_pressure_part[part] ? &d_P_vecs->get("dP_dt", part) : nullptr,
                           current_time,
                           new_time,
                           part);
    }

    // Copy F^{n} into F^{n+1} to make sure that it is stored in the viz file.
    d_F_vecs->copy("current", { "new" });
}

void
FEMechanicsExplicitIntegrator::modifiedEulerStep(const double current_time, const double new_time)
{
    int ierr;
    const double dt = new_time - current_time;

    // U^{n+1} := U^{n} + (dt/rho) F^{n}
    // X^{n+1} := X^{n} + dt       U^{n+1}
    computeLagrangianForce(current_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("current", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("new", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(
                d_P_vecs->get("tmp", part), d_X_vecs->get("new", part), d_U_vecs->get("new", part), new_time, part);
            ierr = VecWAXPY(d_P_vecs->get("new", part).vec(),
                            dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Copy F^{n} into F^{n+1} to make sure that it is stored in the viz file.
    d_F_vecs->copy("current", { "new" });
}

void
FEMechanicsExplicitIntegrator::backwardEulerStep(const double current_time, const double new_time)
{
    int ierr;
    const double dt = new_time - current_time;

    // U^{n+1} := U^{n} + (dt/rho) F^{n+1}
    // X^{n+1} := X^{n} + dt       U^{n+1}
    computeLagrangianForce(new_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("new", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("new", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(
                d_P_vecs->get("tmp", part), d_X_vecs->get("new", part), d_U_vecs->get("new", part), new_time, part);
            ierr = VecWAXPY(d_P_vecs->get("new", part).vec(),
                            dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }
}

void
FEMechanicsExplicitIntegrator::midpointStep(const double current_time, const double new_time)
{
    int ierr;
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;

    // Step 1:
    //    U^{n+1/2} := U^{n} + (dt/(2 rho)) F^{n}
    //    X^{n+1/2} := X^{n} + (dt/2)       U^{n}
    computeLagrangianForce(current_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("half", part).vec(),
                        0.5 * dt / d_rhos[part],
                        d_F_vecs->get("current", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("half", part).vec(),
                        0.5 * dt,
                        d_U_vecs->get("current", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(d_P_vecs->get("tmp", part),
                                               d_X_vecs->get("current", part),
                                               d_U_vecs->get("current", part),
                                               current_time,
                                               part);
            ierr = VecWAXPY(d_P_vecs->get("half", part).vec(),
                            0.5 * dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Step 2:
    //    U^{n+1} := U^{n} + (dt/rho) F^{n+1/2}
    //    X^{n+1} := X^{n} + dt       U^{n+1/2}
    computeLagrangianForce(half_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("half", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("half", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(d_P_vecs->get("tmp", part),
                                               d_X_vecs->get("half", part),
                                               d_U_vecs->get("half", part),
                                               current_time,
                                               part);
            ierr = VecWAXPY(d_P_vecs->get("new", part).vec(),
                            dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Copy F^{n+1/2} into F^{n+1} to make sure that it is stored in the viz file.
    d_F_vecs->copy("half", { "new" });
}

void
FEMechanicsExplicitIntegrator::trapezoidalStep(const double current_time, const double new_time)
{
    int ierr;
    const double dt = new_time - current_time;

    // Step 1 (forward Euler):
    //    U^{n+1,*} := U^{n} + (dt/rho) F^{n}
    //    X^{n+1,*} := X^{n} + dt       U^{n}
    computeLagrangianForce(current_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("current", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("current", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(d_P_vecs->get("tmp", part),
                                               d_X_vecs->get("current", part),
                                               d_U_vecs->get("current", part),
                                               current_time,
                                               part);
            ierr = VecWAXPY(d_P_vecs->get("new", part).vec(),
                            dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Step 2 (trapezoidal rule "corrector"):
    //    U^{n+1} := U^{n} + (dt/(2 rho)) (F^{n} + F^{n+1,*})
    //    X^{n+1} := X^{n} + (dt/2)       (U^{n} + U^{n+1,*})
    computeLagrangianForce(new_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(
                d_P_vecs->get("tmp", part), d_X_vecs->get("new", part), d_U_vecs->get("new", part), new_time, part);
        }
        ierr = VecAXPBYPCZ(d_F_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_F_vecs->get("current", part).vec(),
                           d_F_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_U_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_U_vecs->get("current", part).vec(),
                           d_U_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);

        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("half", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("half", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);

        if (d_dynamic_pressure_part[part])
        {
            ierr = VecAYPX(d_P_vecs->get("tmp", part).vec(), dt, d_P_vecs->get("new", part).vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBYPCZ(d_P_vecs->get("new", part).vec(),
                               0.5,
                               0.5,
                               0.0,
                               d_P_vecs->get("tmp", part).vec(),
                               d_P_vecs->get("current", part).vec());
        }
    }
}

void
FEMechanicsExplicitIntegrator::modifiedTrapezoidalStep(const double current_time, const double new_time)
{
    int ierr;
    const double dt = new_time - current_time;
    // Step 1 (modified Euler):
    //    U^{n+1,*} := U^{n} + (dt/rho) F^{n}
    //    X^{n+1,*} := X^{n} + dt       U^{n+1,*}
    computeLagrangianForce(current_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("current", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("new", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(d_P_vecs->get("tmp", part),
                                               d_X_vecs->get("current", part),
                                               d_U_vecs->get("current", part),
                                               current_time,
                                               part);
            ierr = VecWAXPY(d_P_vecs->get("new", part).vec(),
                            dt,
                            d_P_vecs->get("tmp", part).vec(),
                            d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Step 2 (modified trapezoidal rule "corrector"):
    //    U^{n+1} := U^{n} + (dt/(2 rho)) (F^{n} + F^{n+1,*})
    //    X^{n+1} := X^{n} + (dt/2)       (U^{n} + U^{n+1})
    computeLagrangianForce(new_time);
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        ierr = VecAXPBYPCZ(d_F_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_F_vecs->get("current", part).vec(),
                           d_F_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("half", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);

        ierr = VecAXPBYPCZ(d_U_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_U_vecs->get("current", part).vec(),
                           d_U_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("half", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);

        if (d_dynamic_pressure_part[part])
        {
            computeDynamicPressureRateOfChange(
                d_P_vecs->get("tmp", part), d_X_vecs->get("new", part), d_U_vecs->get("new", part), new_time, part);
            ierr = VecAYPX(d_P_vecs->get("tmp", part).vec(), dt, d_P_vecs->get("new", part).vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBYPCZ(d_P_vecs->get("new", part).vec(),
                               0.5,
                               0.5,
                               0.0,
                               d_P_vecs->get("tmp", part).vec(),
                               d_P_vecs->get("current", part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }
}

void
FEMechanicsExplicitIntegrator::SSPRK3Step(const double current_time, const double new_time, const unsigned int n_stages)
{
    int ierr;
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        if (n_stages == 3)
        {
            // We can avoid allocating vectors for X1, U1, and P1 because these vectors are never used with X_new,
            // U_new, or P_new.
            PetscVector<double>& U_current_vec = d_U_vecs->get("current", part);
            PetscVector<double>& U_new_vec = d_U_vecs->get("new", part);
            PetscVector<double>& U1_vec = U_new_vec;
            PetscVector<double>& U2_vec = d_U_vecs->get("U2", part);

            PetscVector<double>& X_current_vec = d_X_vecs->get("current", part);
            PetscVector<double>& X_new_vec = d_X_vecs->get("new", part);
            PetscVector<double>& X1_vec = X_new_vec;
            PetscVector<double>& X2_vec = d_X_vecs->get("X2", part);

            PetscVector<double>* P_current_vec = partHasPressure(part) ? &d_P_vecs->get("current", part) : nullptr;
            PetscVector<double>* P_new_vec = partHasPressure(part) ? &d_P_vecs->get("new", part) : nullptr;
            PetscVector<double>* P1_vec = P_new_vec;
            PetscVector<double>* P2_vec = d_dynamic_pressure_part[part] ? &d_P_vecs->get("P2", part) : nullptr;

            PetscVector<double>& F_vec = d_F_vecs->get("current", part);

            PetscVector<double>* dP_dt_vec = d_dynamic_pressure_part[part] ? &d_P_vecs->get("dP_dt", part) : nullptr;

            doForwardEulerStep(X1_vec,
                               U1_vec,
                               P1_vec,
                               X_current_vec,
                               U_current_vec,
                               P_current_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time,
                               current_time + 1.0 * dt,
                               part);

            doForwardEulerStep(X2_vec,
                               U2_vec,
                               P2_vec,
                               X1_vec,
                               U1_vec,
                               P1_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time + 1.0 * dt,
                               current_time + 2.0 * dt,
                               part);

            ierr = VecAXPBY(U2_vec.vec(), 0.75, 0.25, U_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X2_vec.vec(), 0.75, 0.25, X_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBY(P2_vec->vec(), 0.75, 0.25, P_current_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            doForwardEulerStep(X_new_vec,
                               U_new_vec,
                               P_new_vec,
                               X2_vec,
                               U2_vec,
                               P2_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time + 0.5 * dt,
                               current_time + 1.5 * dt,
                               part);

            ierr = VecAXPBY(U_new_vec.vec(), 1.0 / 3.0, 2.0 / 3.0, U_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_new_vec.vec(), 1.0 / 3.0, 2.0 / 3.0, X_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBY(P_new_vec->vec(), 1.0 / 3.0, 2.0 / 3.0, P_current_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            d_F_vecs->get("new", part) = F_vec;
        }
        else if (n_stages == 4)
        {
            // We can avoid allocating vectors for X1, U1, and P1 because these vectors are never used with X_new,
            // U_new, or P_new.
            PetscVector<double>& U_current_vec = d_U_vecs->get("current", part);
            PetscVector<double>& U_new_vec = d_U_vecs->get("new", part);
            PetscVector<double>& U1_vec = U_new_vec;
            PetscVector<double>& U2_vec = d_U_vecs->get("U2", part);
            PetscVector<double>& U3_vec = d_U_vecs->get("U3", part);

            PetscVector<double>& X_current_vec = d_X_vecs->get("current", part);
            PetscVector<double>& X_new_vec = d_X_vecs->get("new", part);
            PetscVector<double>& X1_vec = X_new_vec;
            PetscVector<double>& X2_vec = d_X_vecs->get("X2", part);
            PetscVector<double>& X3_vec = d_X_vecs->get("X3", part);

            PetscVector<double>* P_current_vec = partHasPressure(part) ? &d_P_vecs->get("current", part) : nullptr;
            PetscVector<double>* P_new_vec = partHasPressure(part) ? &d_P_vecs->get("new", part) : nullptr;
            PetscVector<double>* P1_vec = P_new_vec;
            PetscVector<double>* P2_vec = partHasPressure(part) ? &d_P_vecs->get("P2", part) : nullptr;
            PetscVector<double>* P3_vec = partHasPressure(part) ? &d_P_vecs->get("P3", part) : nullptr;

            PetscVector<double>& F_vec = d_F_vecs->get("current", part);

            PetscVector<double>* dP_dt_vec = d_dynamic_pressure_part[part] ? &d_P_vecs->get("dP_dt", part) : nullptr;

            doForwardEulerStep(X1_vec,
                               U1_vec,
                               P1_vec,
                               X_current_vec,
                               U_current_vec,
                               P_current_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time,
                               current_time + 1.0 * dt,
                               part);

            ierr = VecAXPBY(U1_vec.vec(), 0.5, 0.5, U_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X1_vec.vec(), 0.5, 0.5, X_current_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBY(P1_vec->vec(), 0.5, 0.5, P_current_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            doForwardEulerStep(X2_vec,
                               U2_vec,
                               P2_vec,
                               X1_vec,
                               U1_vec,
                               P1_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time + 0.5 * dt,
                               current_time + 1.5 * dt,
                               part);

            ierr = VecAXPBY(U2_vec.vec(), 0.5, 0.5, U1_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X2_vec.vec(), 0.5, 0.5, X1_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBY(P2_vec->vec(), 0.5, 0.5, P1_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            doForwardEulerStep(X3_vec,
                               U3_vec,
                               P3_vec,
                               X2_vec,
                               U2_vec,
                               P2_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time + 1.0 * dt,
                               current_time + 2.0 * dt,
                               part);

            ierr = VecAXPBYPCZ(U3_vec.vec(), 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, U_current_vec.vec(), U2_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBYPCZ(X3_vec.vec(), 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, X_current_vec.vec(), X2_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBYPCZ(P3_vec->vec(), 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, P_current_vec->vec(), P2_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            doForwardEulerStep(X_new_vec,
                               U_new_vec,
                               P_new_vec,
                               X3_vec,
                               U3_vec,
                               P3_vec,
                               F_vec,
                               dP_dt_vec,
                               current_time + 0.5 * dt,
                               current_time + 1.5 * dt,
                               part);

            ierr = VecAXPBY(U_new_vec.vec(), 0.5, 0.5, U3_vec.vec());
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_new_vec.vec(), 0.5, 0.5, X3_vec.vec());
            IBTK_CHKERRQ(ierr);
            if (d_dynamic_pressure_part[part])
            {
                ierr = VecAXPBY(P_new_vec->vec(), 0.5, 0.5, P3_vec->vec());
                IBTK_CHKERRQ(ierr);
            }

            d_F_vecs->get("new", part) = F_vec;
        }
        else
        {
            TBOX_ERROR("FEMechanicsExplicitIntegrator::SSPRK3Step(): unsupported number of stages: "
                       << n_stages << "\n  supported options are n_stages = 3 or 4.");
        }
    }
}

void
FEMechanicsExplicitIntegrator::computeLagrangianForce(const double data_time)
{
    const std::string data_time_str = get_data_time_str(data_time, d_current_time, d_new_time);
    batch_vec_ghost_update(d_X_vecs->get(data_time_str), INSERT_VALUES, SCATTER_FORWARD);
    d_F_vecs->zero("RHS Vector");
    d_F_vecs->zero("tmp");
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        if (d_static_pressure_part[part])
        {
            computeStaticPressure(
                d_P_vecs->get(data_time_str, part), d_X_vecs->get(data_time_str, part), data_time, part);
        }
        assembleInteriorForceDensityRHS(d_F_vecs->get("RHS Vector", part),
                                        d_X_vecs->get(data_time_str, part),
                                        d_P_vecs ? &d_P_vecs->get(data_time_str, part) : nullptr,
                                        data_time,
                                        part);
    }
    batch_vec_ghost_update(d_F_vecs->get("RHS Vector"), ADD_VALUES, SCATTER_REVERSE);
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        d_fe_projectors[part]->computeL2Projection(d_F_vecs->get("solution", part),
                                                   d_F_vecs->get("RHS Vector", part),
                                                   getForceSystemName(),
                                                   d_use_consistent_mass_matrix,
                                                   /*close_U*/ false,
                                                   /*close_F*/ false);
        const double eta = d_etas[part];
        if (eta != 0.0)
        {
            int ierr = VecAXPY(d_F_vecs->get("solution", part).vec(), -eta, d_U_vecs->get(data_time_str, part).vec());
            IBTK_CHKERRQ(ierr);
        }
    }
    d_F_vecs->copy("solution", { data_time_str });
}

void
FEMechanicsExplicitIntegrator::computeLagrangianForce(PetscVector<double>& F_vec,
                                                      PetscVector<double>& X_vec,
                                                      PetscVector<double>& U_vec,
                                                      PetscVector<double>* P_vec,
                                                      const double data_time,
                                                      const unsigned int part)
{
    batch_vec_ghost_update({ &X_vec }, INSERT_VALUES, SCATTER_FORWARD);
    auto& F_tmp_vec = d_F_vecs->get("RHS Vector", part);
    auto& F_rhs_vec = d_F_vecs->get("tmp", part);
    F_tmp_vec.zero();
    F_rhs_vec.zero();
    if (d_static_pressure_part[part])
    {
        TBOX_ASSERT(P_vec);
        computeStaticPressure(*P_vec, X_vec, data_time, part);
    }
    assembleInteriorForceDensityRHS(F_rhs_vec, X_vec, P_vec, data_time, part);
    batch_vec_ghost_update({ &F_rhs_vec }, ADD_VALUES, SCATTER_REVERSE);
    d_fe_projectors[part]->computeL2Projection(F_vec,
                                               F_rhs_vec,
                                               getForceSystemName(),
                                               d_use_consistent_mass_matrix,
                                               /*close_U*/ false,
                                               /*close_F*/ false);
    const double eta = d_etas[part];
    if (eta != 0.0)
    {
        int ierr = VecAXPY(F_vec.vec(), -eta, U_vec.vec());
        IBTK_CHKERRQ(ierr);
    }
}

void
FEMechanicsExplicitIntegrator::putToDatabase(Pointer<Database> db)
{
    FEMechanicsBase::putToDatabase(db);
    db->putInteger("EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION", EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION);
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FEMechanicsExplicitIntegrator::doInitializeFEEquationSystems()
{
    FEMechanicsBase::doInitializeFEEquationSystems();
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        EquationSystems& equation_systems = *d_equation_systems[part];
        IBTK::setup_system_vectors(&equation_systems,
                                   { getCurrentCoordinatesSystemName(), getVelocitySystemName() },
                                   { "current", "half", "new" },
                                   from_restart);
        IBTK::setup_system_vectors(&equation_systems,
                                   { getForceSystemName() },
                                   { "current", "half", "new", "tmp", "RHS Vector" },
                                   from_restart);
    }
}

void
FEMechanicsExplicitIntegrator::doInitializeFESystemVectors()
{
    FEMechanicsBase::doInitializeFESystemVectors();
    std::vector<EquationSystems*> equation_systems;
    for (const auto& es : d_equation_systems) equation_systems.push_back(es.get());
    d_X_vecs.reset(new LibMeshSystemVectors(equation_systems, getCurrentCoordinatesSystemName()));
    d_U_vecs.reset(new LibMeshSystemVectors(equation_systems, getVelocitySystemName()));
    d_F_vecs.reset(new LibMeshSystemVectors(equation_systems, getForceSystemName()));
    if (d_has_static_pressure_parts)
    {
        d_P_vecs.reset(new LibMeshSystemVectors(equation_systems, d_static_pressure_part, getPressureSystemName()));
    }
    if (d_has_dynamic_pressure_parts)
    {
        d_P_vecs.reset(new LibMeshSystemVectors(equation_systems, d_dynamic_pressure_part, getPressureSystemName()));
    }
}

void
FEMechanicsExplicitIntegrator::doInitializeFEData(const bool use_present_data)
{
    FEMechanicsBase::doInitializeFEData(use_present_data);
}

void
FEMechanicsExplicitIntegrator::doForwardEulerStep(PetscVector<double>& X_new_vec,
                                                  PetscVector<double>& U_new_vec,
                                                  PetscVector<double>* P_new_vec,
                                                  PetscVector<double>& X_current_vec,
                                                  PetscVector<double>& U_current_vec,
                                                  PetscVector<double>* P_current_vec,
                                                  PetscVector<double>& F_current_vec,
                                                  PetscVector<double>* dP_dt_current_vec,
                                                  const double start_time,
                                                  const double end_time,
                                                  const unsigned int part)
{
    TBOX_ASSERT(&X_new_vec != &X_current_vec);
    TBOX_ASSERT(&U_new_vec != &U_current_vec);
    TBOX_ASSERT((P_new_vec == nullptr && P_current_vec == nullptr) || P_new_vec != P_current_vec);

    const double dt = end_time - start_time;

    // U^{n+1} := U^{n} + (dt/rho) F^{n}
    // X^{n+1} := X^{n} + dt       U^{n}
    computeLagrangianForce(F_current_vec, X_current_vec, U_current_vec, P_current_vec, start_time, part);
    int ierr;
    ierr = VecWAXPY(U_new_vec.vec(), dt / d_rhos[part], F_current_vec.vec(), U_current_vec.vec());
    IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(X_new_vec.vec(), dt, U_current_vec.vec(), X_current_vec.vec());
    IBTK_CHKERRQ(ierr);
    if (d_dynamic_pressure_part[part])
    {
        TBOX_ASSERT(P_new_vec);
        TBOX_ASSERT(P_current_vec);
        TBOX_ASSERT(dP_dt_current_vec);
        computeDynamicPressureRateOfChange(*dP_dt_current_vec, X_current_vec, U_current_vec, start_time, part);
        ierr = VecWAXPY(P_new_vec->vec(), dt, dP_dt_current_vec->vec(), P_current_vec->vec());
        IBTK_CHKERRQ(ierr);
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsExplicitIntegrator::commonConstructor(const Pointer<Database>& input_db)
{
    // Set some default values.
    d_rhos.resize(d_meshes.size(), 1.0);
    d_etas.resize(d_meshes.size(), 0.0);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Determine how to compute weak forms.
    d_include_normal_stress_in_weak_form = false;
    d_include_tangential_stress_in_weak_form = false;
    d_include_normal_surface_forces_in_weak_form = true;
    d_include_tangential_surface_forces_in_weak_form = true;
}

void
FEMechanicsExplicitIntegrator::getFromInput(const Pointer<Database>& db, bool /*is_from_restart*/)
{
    // Problem parameters.
    if (db->isDouble("mass_density"))
        std::fill(d_rhos.begin(), d_rhos.end(), db->getDouble("mass_density"));
    else if (db->keyExists("mass_density"))
    {
        TBOX_ASSERT(static_cast<std::size_t>(db->getArraySize("mass_density")) == d_rhos.size());
        db->getDoubleArray("mass_density", d_rhos.data(), db->getArraySize("mass_density"));
    }
    if (db->isDouble("damping"))
        std::fill(d_etas.begin(), d_etas.end(), db->getDouble("damping"));
    else if (db->keyExists("damping"))
    {
        TBOX_ASSERT(static_cast<std::size_t>(db->getArraySize("damping")) == d_etas.size());
        db->getDoubleArray("damping", d_etas.data(), db->getArraySize("damping"));
    }
}

void
FEMechanicsExplicitIntegrator::getFromRestart()
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
    int ver = db->getInteger("EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION");
    if (ver != EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
