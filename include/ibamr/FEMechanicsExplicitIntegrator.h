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
 * IBStrategy that provides functionality required by the IB method with finite
 * element elasticity.
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
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time);

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
