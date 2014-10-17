// Filename: GeneralizedIBMethod.h
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

#ifndef included_GeneralizedIBMethod
#define included_GeneralizedIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "Variable.h"
#include "ibamr/IBKirchhoffRodForceGen.h"
#include "ibamr/IBMethod.h"
#include "tbox/Pointer.h"

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
    GeneralizedIBMethod(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~GeneralizedIBMethod();

    /*!
     * Supply a Lagrangian force object.
     */
    void registerIBKirchhoffRodForceGen(SAMRAI::tbox::Pointer<IBKirchhoffRodForceGen> ib_force_and_torque_fcn);

    /*!
     * Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables();

    /*!
     * Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    void registerEulerianCommunicationAlgorithms();

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void eulerStep(double current_time, double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time);

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time);

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
        bool initial_time);

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
                             bool allocate_data);

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

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
    bool d_N_current_needs_ghost_fill, d_N_new_needs_ghost_fill;

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
    GeneralizedIBMethod();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    GeneralizedIBMethod(const GeneralizedIBMethod& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    GeneralizedIBMethod& operator=(const GeneralizedIBMethod& that);

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

#endif //#ifndef included_GeneralizedIBMethod
