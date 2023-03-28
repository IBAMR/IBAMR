// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_GeneralizedIBMethod
#define included_IBAMR_GeneralizedIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBKirchhoffRodForceGen.h"
#include "ibamr/IBMethod.h"

#include "ibtk/LData.h"

#include "Variable.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace tbox
{
class Database;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class GeneralizedIBMethod is an extension of class IBMethod that
 * provides functionality required by the generalized IB (gIB) method.
 */
class GeneralizedIBMethod : public IBMethod
{
public:
    /*!
     * \brief Constructor.
     */
    GeneralizedIBMethod(std::string object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~GeneralizedIBMethod() = default;

    /*!
     * Supply a Lagrangian force object.
     */
    void registerIBKirchhoffRodForceGen(SAMRAI::tbox::Pointer<IBKirchhoffRodForceGen> ib_force_and_torque_fcn);

    /*!
     * Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables() override;

    /*!
     * Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    void registerEulerianCommunicationAlgorithms() override;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time) override;

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time) override;

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time) override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                             bool allocate_data) override;

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_f_var, d_w_var, d_n_var;
    int d_f_idx, d_w_idx, d_n_idx;

    /*
     * Boolean values tracking whether certain quantities need to be
     * reinitialized.
     */
    bool d_N_current_needs_ghost_fill = true, d_N_new_needs_ghost_fill = true;

    /*
     * Lagrangian variables.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_D_current_data, d_D_new_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_N_current_data, d_N_new_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_W_current_data, d_W_new_data;

    /*
     * The force and torque generator.
     */
    SAMRAI::tbox::Pointer<IBKirchhoffRodForceGen> d_ib_force_and_torque_fcn;
    bool d_ib_force_and_torque_fcn_needs_init;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    GeneralizedIBMethod() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    GeneralizedIBMethod(const GeneralizedIBMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    GeneralizedIBMethod& operator=(const GeneralizedIBMethod& that) = delete;

    /*!
     * Reset the Lagrangian force function object.
     */
    void resetLagrangianForceAndTorqueFunction(double init_data_time, bool initial_time);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_GeneralizedIBMethod
