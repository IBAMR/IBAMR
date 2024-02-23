// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_AdvDiffSemiImplicitHierarchyIntegrator
#define included_IBAMR_AdvDiffSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "HierarchyFaceDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <set>
#include <string>

namespace IBAMR
{
class ConvectiveOperator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace pdat
{
template <int DIM, class TYPE>
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
                                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                           bool register_for_restart = true);

    /*!
     * The destructor for class AdvDiffSemiImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~AdvDiffSemiImplicitHierarchyIntegrator() = default;

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
    void setDefaultConvectiveOperatorInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Get the default convective operator input database to be used by
     * the solver.
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> getDefaultConvectiveOperatorInputDatabase() const;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator. Can optionally turn off outputting the quantity.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    void registerTransportedQuantity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                     const bool Q_output = true) override;

    /*!
     * Set the convective time stepping method to use for a particular
     * transported quantity Q.
     */
    void setConvectiveTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                       TimeSteppingType convective_time_stepping_type);

    /*!
     * Get the convective time stepping method being used for a particular
     * transported quantity Q.
     */
    TimeSteppingType
    getConvectiveTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the convective time stepping method used during the initial time step
     * for a particular transported quantity Q.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    void setInitialConvectiveTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                              TimeSteppingType init_convective_time_stepping_type);

    /*!
     * Return the convective time stepping method used during the initial time
     * step for a particular transported quantity Q.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    TimeSteppingType
    getInitialConvectiveTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * \brief Set the convective operator type to be used by the solver for a
     * particular transported quantity Q.
     */
    void setConvectiveOperatorType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                   const std::string& op_type);

    /*!
     * \brief Get the convective operator type used by the solver for a
     * particular transported quantity Q.
     */
    const std::string&
    getConvectiveOperatorType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * \brief Set the convective operator input database to be used by the
     * solver for a particular transported quantity Q.
     */
    void setConvectiveOperatorInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Get the convective operator boundary input database used by the
     * solver for a particular transported quantity Q.
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>
    getConvectiveOperatorInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Register an operator to compute the convective derivative term u*grad Q
     * for a particular transported quantity Q.
     */
    void setConvectiveOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                               SAMRAI::tbox::Pointer<ConvectiveOperator> convective_op);

    /*!
     * Get the convective operator being used by this solver class for a
     * particular transported quantity Q.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    virtual SAMRAI::tbox::Pointer<ConvectiveOperator>
    getConvectiveOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Indicate that all of the convective operators should be (re-)initialized
     * before the next time step.
     */
    void setConvectiveOperatorsNeedInit();

    /*!
     * Indicate that the convective operator should be (re-)initialized before
     * the next time step.
     */
    void setConvectiveOperatorNeedsInit(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

protected:
    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Default convective time integration methods.
     */
    TimeSteppingType d_default_convective_time_stepping_type;
    TimeSteppingType d_default_init_convective_time_stepping_type;

    /*!
     * Default convective operator settings.
     */
    std::string d_default_convective_op_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_default_convective_op_input_db;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;

    /*!
     * Transported quantities.
     */
    std::set<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_N_var, d_N_old_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > >
        d_Q_N_map, d_Q_N_old_map;

    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, TimeSteppingType>
        d_Q_convective_time_stepping_type, d_Q_init_convective_time_stepping_type;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, std::string> d_Q_convective_op_type;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> >
        d_Q_convective_op_input_db;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<ConvectiveOperator> >
        d_Q_convective_op;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, bool> d_Q_convective_op_needs_init;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffSemiImplicitHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator(const AdvDiffSemiImplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffSemiImplicitHierarchyIntegrator& operator=(const AdvDiffSemiImplicitHierarchyIntegrator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data are read is determined
     * by the object_name specified in the class constructor.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AdvDiffSemiImplicitHierarchyIntegrator
