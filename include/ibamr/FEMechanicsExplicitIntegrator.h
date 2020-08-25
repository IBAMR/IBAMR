// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibamr/FEMechanicsBase.h"

#include "ibtk/FEProjector.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class FEMechanicsExplicitIntegrator is an implementation of the abstract base class
 * FEMechanicsBase that provides a simple explicit elastodynamics time integrator with an
 * interface that is similar to IBFEMethod to facilitate model re-use.
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
     * Advance the structural velocities and positions using the explicit midpoint rule.
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
     * velocity is used instead of the intermediate updated approximation used in the
     * standard SSP RK2 scheme.
     */
    void modifiedTrapezoidalStep(double current_time, double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time);

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
     * Do the actual work in reinitializeFEData and initializeFEData. if @p
     * use_present_data is `true` then the current content of the solution
     * vectors is used: more exactly, the coordinates and velocities (computed
     * by initializeCoordinates and initializeVelocity) are considered as
     * being up to date, as is the direct forcing kinematic data.
     */
    void doInitializeFEData(bool use_present_data) override;

    /// Structure mass densities.
    std::vector<double> d_rhos;

    /// FEProjector objects.
    std::vector<std::unique_ptr<IBTK::FEProjector> > d_fe_projectors;

private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                           const std::vector<libMesh::MeshBase*>& meshes,
                           const std::string& restart_read_dirname,
                           unsigned int restart_restore_number);

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
