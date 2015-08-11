// Filename: AdvDiffSemiImplicitHierarchyIntegrator.h
// Created on 22 May 2012 by Boyce Griffith
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

#ifndef included_AdvDiffSemiImplicitHierarchyIntegrator
#define included_AdvDiffSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <set>
#include <string>

#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "SAMRAI/tbox/Database.h"

namespace IBAMR
{
class ConvectiveOperator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
class PatchHierarchy;
class PatchHierarchy;
} // namespace hier
namespace mesh
{
class GriddingAlgorithm;
} // namespace mesh
namespace pdat
{
template <class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffSemiImplicitHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the advection-diffusion equation.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \kappa \f$ and damping coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Additionally, a different
 * advection velocity may be used with each quantity registered with the
 * integrator.
 *
 * This hierarchy integrator advances all levels of the patch hierarchy
 * synchronously in time.  In particular, subcycling in time is \em not
 * performed.
 *
 * Various options are available for the spatial and temporal discretizations.
 *
 * \see HierarchyIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class AdvDiffSemiImplicitHierarchyIntegrator : public AdvDiffHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffSemiImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    AdvDiffSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                           const boost::shared_ptr<SAMRAI::tbox::Database>& input_db,
                                           bool register_for_restart = true);

    /*!
     * The destructor for class AdvDiffSemiImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~AdvDiffSemiImplicitHierarchyIntegrator();

    /*!
     * Set the default convective time stepping method to use with registered
     * quantities.
     */
    void setDefaultConvectiveTimeSteppingType(TimeSteppingType default_convective_time_stepping_type);

    /*!
     * Get the default convective time stepping method being used with registered
     * quantities.
     */
    TimeSteppingType getDefaultConvectiveTimeSteppingType() const;

    /*!
     * Set the default convective time stepping method used during the initial
     * time step.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    void setDefaultInitialConvectiveTimeSteppingType(TimeSteppingType default_init_convective_time_stepping_type);

    /*!
     * Return the default convective time stepping method used during the
     * initial time step.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    TimeSteppingType getDefaultInitialConvectiveTimeSteppingType() const;

    /*!
     * \brief Set the default convective operator type to be used by the solver.
     */
    void setDefaultConvectiveOperatorType(const std::string& op_type);

    /*!
     * \brief Get the default convective operator type used by the solver.
     */
    const std::string& getDefaultConvectiveOperatorType() const;

    /*!
     * \brief Set the default convective operator input database to be used by
     * the solver.
     */
    void setDefaultConvectiveOperatorInputDatabase(const boost::shared_ptr<SAMRAI::tbox::Database>& input_db);

    /*!
     * \brief Get the default convective operator input database to be used by
     * the solver.
     */
    boost::shared_ptr<SAMRAI::tbox::Database> getDefaultConvectiveOperatorInputDatabase() const;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    void registerTransportedQuantity(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var);

    /*!
     * Set the convective time stepping method to use for a particular
     * transported quantity Q.
     */
    void setConvectiveTimeSteppingType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var,
                                       TimeSteppingType convective_time_stepping_type);

    /*!
     * Get the convective time stepping method being used for a particular
     * transported quantity Q.
     */
    TimeSteppingType
    getConvectiveTimeSteppingType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var) const;

    /*!
     * Set the convective time stepping method used during the initial time step
     * for a particular transported quantity Q.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    void setInitialConvectiveTimeSteppingType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var,
                                              TimeSteppingType init_convective_time_stepping_type);

    /*!
     * Return the convective time stepping method used during the initial time
     * step for a particular transported quantity Q.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    TimeSteppingType
    getInitialConvectiveTimeSteppingType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var) const;

    /*!
     * \brief Set the convective operator type to be used by the solver for a
     * particular transported quantity Q.
     */
    void setConvectiveOperatorType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var,
                                   const std::string& op_type);

    /*!
     * \brief Get the convective operator type used by the solver for a
     * particular transported quantity Q.
     */
    const std::string&
    getConvectiveOperatorType(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var) const;

    /*!
     * \brief Set the convective operator input database to be used by the
     * solver for a particular transported quantity Q.
     */
    void setConvectiveOperatorInputDatabase(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var,
                                            const boost::shared_ptr<SAMRAI::tbox::Database>& input_db);

    /*!
     * \brief Get the convective operator boundary input database used by the
     * solver for a particular transported quantity Q.
     */
    boost::shared_ptr<SAMRAI::tbox::Database>
    getConvectiveOperatorInputDatabase(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var) const;

    /*!
     * Register an operator to compute the convective derivative term u*grad Q
     * for a particular transported quantity Q.
     */
    void setConvectiveOperator(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var,
                               const boost::shared_ptr<ConvectiveOperator>& convective_op);

    /*!
     * Get the convective operator being used by this solver class for a
     * particular transported quantity Q.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    virtual boost::shared_ptr<ConvectiveOperator>
    getConvectiveOperator(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var);

    /*!
     * Indicate that all of the convective operators should be (re-)initialized
     * before the next time step.
     */
    void setConvectiveOperatorsNeedInit();

    /*!
     * Indicate that the convective operator should be (re-)initialized before
     * the next time step.
     */
    void setConvectiveOperatorNeedsInit(const boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>& Q_var);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                       const boost::shared_ptr<SAMRAI::mesh::GriddingAlgorithm>& gridding_alg);

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

protected:
    /*!
     * Reset cached hierarchy dependent data.
     */
    void resetHierarchyConfigurationSpecialized(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                                int coarsest_level,
                                                int finest_level);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToRestartSpecialized(const boost::shared_ptr<SAMRAI::tbox::Database>& db) const;

    /*!
     * Default convective time integration methods.
     */
    TimeSteppingType d_default_convective_time_stepping_type;
    TimeSteppingType d_default_init_convective_time_stepping_type;

    /*!
     * Default convective operator settings.
     */
    std::string d_default_convective_op_type;
    boost::shared_ptr<SAMRAI::tbox::Database> d_default_convective_op_input_db;

    /*
     * Hierarchy operations objects.
     */
    boost::shared_ptr<SAMRAI::math::HierarchyFaceDataOpsReal<double>> d_hier_fc_data_ops;

    /*!
     * Transported quantities.
     */
    std::set<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>> d_N_var, d_N_old_var;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>> d_Q_N_map,
        d_Q_N_old_map;

    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>, TimeSteppingType> d_Q_convective_time_stepping_type,
        d_Q_init_convective_time_stepping_type;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>, std::string> d_Q_convective_op_type;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>, boost::shared_ptr<SAMRAI::tbox::Database>>
        d_Q_convective_op_input_db;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>, boost::shared_ptr<ConvectiveOperator>>
        d_Q_convective_op;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double>>, bool> d_Q_convective_op_needs_init;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffSemiImplicitHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator(const AdvDiffSemiImplicitHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator& operator=(const AdvDiffSemiImplicitHierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(const boost::shared_ptr<SAMRAI::tbox::Database>& db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data are read is determined
     * by the object_name specified in the class constructor.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffSemiImplicitHierarchyIntegrator
