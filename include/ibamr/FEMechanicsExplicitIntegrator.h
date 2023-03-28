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

#ifndef included_IBAMR_FEMechanicsExplicitIntegrator
#define included_IBAMR_FEMechanicsExplicitIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/FEMechanicsBase.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class FEMechanicsExplicitIntegrator is an implementation of the
 * abstract base class FEMechanicsBase that provides a simple explicit
 * elastodynamics time integrator with an interface that is similar to
 * IBFEMethod to facilitate model re-use.
 *
 * \see IBFEMethod
 */
class FEMechanicsExplicitIntegrator : public FEMechanicsBase
{
public:
    /*!
     * \brief Constructor for a single-part model.
     */
    FEMechanicsExplicitIntegrator(const std::string& object_name,
                                  const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                                  libMesh::MeshBase* mesh,
                                  bool register_for_restart = true,
                                  const std::string& restart_read_dirname = "",
                                  unsigned int restart_restore_number = 0);

    /*!
     * \brief Constructor for a multi-part model.
     */
    FEMechanicsExplicitIntegrator(const std::string& object_name,
                                  const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                                  const std::vector<libMesh::MeshBase*>& meshes,
                                  bool register_for_restart = true,
                                  const std::string& restart_read_dirname = "",
                                  unsigned int restart_restore_number = 0);

    /*!
     * \brief Deleted default constructor.
     */
    FEMechanicsExplicitIntegrator() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    FEMechanicsExplicitIntegrator(const FEMechanicsExplicitIntegrator& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    FEMechanicsExplicitIntegrator& operator=(const FEMechanicsExplicitIntegrator& that) = delete;

    /*!
     * \brief Defaulted destructor.
     */
    ~FEMechanicsExplicitIntegrator() override = default;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Advance the structural velocities and positions using forward Euler.
     */
    void forwardEulerStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using "modified" Euler.
     *
     * NOTE: This scheme first advances the velocity using forward Euler, and
     * then uses the updated velocity to advance the position of the structure.
     */
    void modifiedEulerStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using backward Euler.
     *
     * NOTE: This is not a true backward Euler method; we simply use the most
     * recent approximations to the updated velocity and position to compute the
     * force at the end of the time step, and then use that force to increment
     * the velocity and position via backward Euler.
     */
    void backwardEulerStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using the explicit midpoint
     * rule.
     */
    void midpointStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using the explicit SSP
     * RK2 scheme, which is the explicit trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using a modified explicit
     * trapezoidal rule.
     *
     * NOTE: In this scheme, the first stage uses a modified Euler scheme (\see
     * modifiedEulerStep), and the second stage uses a modified explicit
     * trapezoidal rule, in which the final approximation to the updated
     * velocity is used instead of the intermediate updated approximation used in
     * the standard SSP RK2 scheme.
     */
    void modifiedTrapezoidalStep(double current_time, double new_time);

    /*!
     * Advance the structural velocities and positions using the explicit SSP RK3 scheme.
     *
     * The implementation supports either 3- or 4-stage variants. The default is the 4-stage
     * version, which has a larger stability region.
     */
    void SSPRK3Step(double current_time, double new_time, unsigned int n_stages = 4);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval using the provided data vectors.
     */
    void computeLagrangianForce(libMesh::PetscVector<double>& F_vec,
                                libMesh::PetscVector<double>& X_vec,
                                libMesh::PetscVector<double>& U_vec,
                                libMesh::PetscVector<double>* P_vec,
                                double data_time,
                                unsigned int part);

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * Do the actual work in initializeFEEquationSystems.
     */
    void doInitializeFEEquationSystems() override;

    /*!
     * Do the actual work of setting up libMesh system vectors.
     */
    void doInitializeFESystemVectors() override;

    /*!
     * Do the actual work in reinitializeFEData and initializeFEData. if @p
     * use_present_data is `true` then the current content of the solution
     * vectors is used: more exactly, the coordinates and velocities (computed
     * by initializeCoordinates and initializeVelocity) are considered as
     * being up to date, as is the direct forcing kinematic data.
     */
    void doInitializeFEData(bool use_present_data) override;

    /*!
     * Perform a forward Euler step.
     */
    void doForwardEulerStep(libMesh::PetscVector<double>& X_new_vec,
                            libMesh::PetscVector<double>& U_new_vec,
                            libMesh::PetscVector<double>* P_new_vec,
                            libMesh::PetscVector<double>& X_current_vec,
                            libMesh::PetscVector<double>& U_current_vec,
                            libMesh::PetscVector<double>* P_current_vec,
                            libMesh::PetscVector<double>& F_current_vec,
                            libMesh::PetscVector<double>* dP_dt_current_vec,
                            double current_time,
                            double new_time,
                            unsigned int part);

    /// Structure mass densities.
    std::vector<double> d_rhos;

    /// Structure damping coefficients.
    std::vector<double> d_etas;

private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_FEMechanicsExplicitIntegrator
