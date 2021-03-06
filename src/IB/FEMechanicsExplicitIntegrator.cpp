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

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/FEMechanicsExplicitIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/FEProjector.h"
#include "ibtk/libmesh_utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"

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
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        // U^{n+1} := U^{n} + (dt/rho) F^{n}
        // X^{n+1} := X^{n} + dt       U^{n}
        computeLagrangianForce(current_time);
        int ierr;
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
    // We copy F^{n} into F^{n+1} to make sure that it is stored in the viz
    // file.
    d_F_vecs->copy("current", { "new" });
}

void
FEMechanicsExplicitIntegrator::modifiedEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        // U^{n+1} := U^{n} + (dt/rho) F^{n}
        // X^{n+1} := X^{n} + dt       U^{n+1}
        computeLagrangianForce(current_time);
        int ierr;
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
    // We copy F^{n} into F^{n+1} to make sure that it is stored in the viz
    // file.
    d_F_vecs->copy("current", { "new" });
}

void
FEMechanicsExplicitIntegrator::backwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        // U^{n+1} := U^{n} + (dt/rho) F^{n+1}
        // X^{n+1} := X^{n} + dt       U^{n+1}
        computeLagrangianForce(new_time);
        int ierr;
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
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr;

        // Step 1:
        //    U^{n+1/2} := U^{n} + (dt/(2 rho)) F^{n}
        //    X^{n+1/2} := X^{n} + (dt/2)       U^{n}
        computeLagrangianForce(current_time);
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

        // Step 2:
        //    U^{n+1} := U^{n} + (dt/rho) F^{n+1/2}
        //    X^{n+1} := X^{n} + dt       U^{n+1/2}
        computeLagrangianForce(half_time);
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
    // We copy F^{n+1/2} into F^{n+1} to make sure that it is stored in the viz
    // file.
    d_F_vecs->copy("half", { "new" });
}

void
FEMechanicsExplicitIntegrator::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr;

        // Step 1 (forward Euler):
        //    U^{n+1,*} := U^{n} + (dt/rho) F^{n}
        //    X^{n+1,*} := X^{n} + dt       U^{n}
        computeLagrangianForce(current_time);
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

        // Step 2 (trapezoidal rule "corrector"):
        //    U^{n+1} := U^{n} + (dt/(2 rho)) (F^{n} + F^{n+1,*})
        //    X^{n+1} := X^{n} + (dt/2)       (U^{n} + U^{n+1,*})
        computeLagrangianForce(new_time);
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
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        int ierr;

        // Step 1 (modified Euler):
        //    U^{n+1,*} := U^{n} + (dt/rho) F^{n}
        //    X^{n+1,*} := X^{n} + dt       U^{n+1,*}
        computeLagrangianForce(current_time);
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

        // Step 2 (modified trapezoidal rule "corrector"):
        //    U^{n+1} := U^{n} + (dt/(2 rho)) (F^{n} + F^{n+1,*})
        //    X^{n+1} := X^{n} + (dt/2)       (U^{n} + U^{n+1})
        computeLagrangianForce(new_time);
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
                                                   FORCE_SYSTEM_NAME,
                                                   d_use_consistent_mass_matrix,
                                                   /*close_U*/ false,
                                                   /*close_F*/ false);
    }
    d_F_vecs->copy("solution", { data_time_str });
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
                                   { COORDS_SYSTEM_NAME, VELOCITY_SYSTEM_NAME },
                                   { "current", "half", "new" },
                                   from_restart);
        IBTK::setup_system_vectors(
            &equation_systems, { FORCE_SYSTEM_NAME }, { "current", "half", "new", "tmp", "RHS Vector" }, from_restart);
    }
}

void
FEMechanicsExplicitIntegrator::doInitializeFESystemVectors()
{
    FEMechanicsBase::doInitializeFESystemVectors();
    std::vector<EquationSystems*> equation_systems;
    for (const auto& es : d_equation_systems) equation_systems.push_back(es.get());
    d_X_vecs.reset(new LibMeshSystemVectors(equation_systems, COORDS_SYSTEM_NAME));
    d_U_vecs.reset(new LibMeshSystemVectors(equation_systems, VELOCITY_SYSTEM_NAME));
    d_F_vecs.reset(new LibMeshSystemVectors(equation_systems, FORCE_SYSTEM_NAME));
    if (d_has_static_pressure_parts)
    {
        d_P_vecs.reset(new LibMeshSystemVectors(equation_systems, d_static_pressure_part, PRESSURE_SYSTEM_NAME));
    }
    if (d_has_dynamic_pressure_parts)
    {
        d_P_vecs.reset(new LibMeshSystemVectors(equation_systems, d_dynamic_pressure_part, PRESSURE_SYSTEM_NAME));
    }
}

void
FEMechanicsExplicitIntegrator::doInitializeFEData(const bool use_present_data)
{
    FEMechanicsBase::doInitializeFEData(use_present_data);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsExplicitIntegrator::commonConstructor(const Pointer<Database>& input_db)
{
    // Set some default values.
    d_rhos.resize(d_meshes.size(), 1.0);

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
