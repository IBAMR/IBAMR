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

#ifndef included_IBAMR_IBImplicitHierarchyIntegrator
#define included_IBAMR_IBImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBImplicitStrategy.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesSolver.h"

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"

#include <limits>
#include <string>

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace solv
{
class PoissonSpecifications;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitHierarchyIntegrator is an implementation of a
 * formally second-order accurate, nonlinearly-implicit version of the immersed
 * boundary method.
 */
class IBImplicitHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    IBImplicitHierarchyIntegrator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_method_ops,
                                  SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
                                  bool register_for_restart = true);

    /*!
     * The destructor for class IBImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~IBImplicitHierarchyIntegrator() = default;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

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

protected:
    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops;

    // Whether to use "frozen" LE operators.
    bool d_use_fixed_LE_operators = false;

    // Whether to solve for the position.
    bool d_solve_for_position = false;

    // Penalty factor used in direct forcing.
    double d_eta = std::numeric_limits<double>::quiet_NaN();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitHierarchyIntegrator(const IBImplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitHierarchyIntegrator& operator=(const IBImplicitHierarchyIntegrator& that) = delete;

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /// Solver state data.
    int d_cycle_num, d_ins_cycle_num;
    double d_current_time, d_new_time;

    /*!
     * Update the solution (e.g. for fixed-point iteration) based on the current value of Y and compute the residual.
     */
    void updateSolution(Vec Y, Vec R);

    /*!
     * Static function for implicit formulation.
     */
    static PetscErrorCode IBFunction(SNES snes, Vec Y, Vec R, void* ctx);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBImplicitHierarchyIntegrator
