// Filename: AdvDiffHierarchyIntegrator.h
// Created on 21 May 2012 by Boyce Griffith
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

#ifndef included_AdvDiffHierarchyIntegrator
#define included_AdvDiffHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <string>
#include <vector>

#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class CartGridFunction;
class LaplaceOperator;
class PoissonSolver;
} // namespace IBTK
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
template <int DIM, class TYPE>
class FaceVariable;
template <int DIM, class TYPE>
class SideVariable;
} // namespace pdat
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffHierarchyIntegrator provides an abstract interface for a
 * time integrator for advection-diffusion or advection-reaction-diffusion
 * equations on an AMR grid hierarchy, along with basic data management for
 * variables defined on that hierarchy.
 */
class AdvDiffHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    /*!
     * The destructor for class AdvDiffHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~AdvDiffHierarchyIntegrator();

    /*!
     * Set the default diffusion time stepping method to use with registered
     * quantities.
     */
    void setDefaultDiffusionTimeSteppingType(TimeSteppingType default_diffusion_time_stepping_type);

    /*!
     * Get the default diffusion time stepping method being used with registered
     * quantities.
     */
    TimeSteppingType getDefaultDiffusionTimeSteppingType() const;

    /*!
     * Set the default convective differencing form to use with registered
     * quantities.
     */
    void setDefaultConvectiveDifferencingType(ConvectiveDifferencingType default_convective_difference_form);

    /*!
     * Get the default convective differencing form being used with registered
     * quantities.
     */
    ConvectiveDifferencingType getDefaultConvectiveDifferencingType() const;

    /*!
     * Register a face-centered advection velocity to be used to advect
     * cell-centered quantities by the hierarchy integrator.
     *
     * Data management for the registered advection velocity will be handled by
     * the hierarchy integrator.
     *
     * \note By default, each registered advection velocity is assumed to be
     * divergence free.
     */
    virtual void registerAdvectionVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

    /*!
     * Indicate whether a particular advection velocity is discretely divergence
     * free.
     */
    void setAdvectionVelocityIsDivergenceFree(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var,
                                              bool is_div_free);

    /*!
     * Determine whether a particular advection velocity has been indicated to
     * be discretely divergence free.
     */
    bool
    getAdvectionVelocityIsDivergenceFree(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var) const;

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular advection velocity.
     */
    void setAdvectionVelocityFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var,
                                      SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_fcn);

    /*!
     * Get the IBTK::CartGridFunction object being used to specify the value of
     * a particular advection velocity.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
    getAdvectionVelocityFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var) const;

    /*!
     * Register a cell-centered source term.
     *
     * Data management for the registered source term will be handled by the
     * hierarchy integrator.
     */
    virtual void registerSourceTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > F_var);

    /*!
     * Supply an IBTK::CartGridFunction object to specify the value of a
     * particular source term.
     */
    void setSourceTermFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > F_var,
                               SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Get the IBTK::CartGridFunction object being used to specify the value of
     * a particular source term.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
    getSourceTermFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > F_var) const;

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator.
     *
     * Data management for the registered quantity will be handled by the
     * hierarchy integrator.
     */
    virtual void registerTransportedQuantity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Set the face-centered advection velocity to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void setAdvectionVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > u_var);

    /*!
     * Get the face-centered advection velocity being used with a particular
     * cell-centered quantity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> >
    getAdvectionVelocity(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the cell-centered source term to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void setSourceTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > F_var);

    /*!
     * Get the cell-centered source term being used with a particular
     * cell-centered quantity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >
    getSourceTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the diffusion time integration scheme for a quantity that has been
     * registered with the hierarchy integrator.
     */
    void setDiffusionTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                      TimeSteppingType time_stepping_type);

    /*!
     * Get the diffusion time integration scheme for a quantity that has been
     * registered with the hierarchy integrator.
     */
    TimeSteppingType
    getDiffusionTimeSteppingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the convective differencing form for a quantity that has been
     * registered with the hierarchy integrator.
     */
    void setConvectiveDifferencingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                       ConvectiveDifferencingType difference_form);

    /*!
     * Get the convective differencing form for a quantity that has been
     * registered with the hierarchy integrator.
     */
    ConvectiveDifferencingType
    getConvectiveDifferencingType(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the constant scalar diffusion coefficient corresponding to a quantity that has
     * been registered with the hierarchy integrator.
     */
    void setDiffusionCoefficient(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var, double kappa);

    /*!
     * Get the constant scalar diffusion coefficient corresponding to a quantity that has
     * been registered with the hierarchy integrator.
     */
    double getDiffusionCoefficient(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Register a variable scalar diffusion coefficient corresponding to a quantity
     * that has been registered with the hierarchy integrator.
     */
    void registerDiffusionCoefficientVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > D_var);

    /*!
     * Supply an IBTK:CartGridFunction object to specify the value of a particular
     * variable diffusion coefficient.
     */
    void setDiffusionCoefficientFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > D_var,
                                         SAMRAI::tbox::Pointer<IBTK::CartGridFunction> D_fcn);

    /*!
     * Get the IBTK::CartGridFunction object being used to specify the value of
     * a particular variable diffusion coefficient.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
    getDiffusionCoefficientFunction(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > D_var) const;

    /*!
     * Set the cell-centered variable diffusion coefficient to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void setDiffusionCoefficientVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                         SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > D_var);

    /*!
     * Get the cell-centered variable diffusion coefficient being used with a particular
     * cell-centered quantity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> >
    getDiffusionCoefficientVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*
     * Return whether the diffusion coefficient being used with a particular
     * cell-centered quantity is variable.
     */
    bool isDiffusionCoefficientVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set the scalar linear damping coefficient corresponding to a quantity
     * that has been registered with the hierarchy integrator.
     */
    void setDampingCoefficient(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var, double lambda);

    /*!
     * Get the scalar linear damping coefficient corresponding to a quantity
     * that has been registered with the hierarchy integrator.
     */
    double getDampingCoefficient(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set a grid function to provide initial conditions for a quantity that has
     * been registered with the hierarchy integrator.
     */
    void setInitialConditions(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                              SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q_init);

    /*!
     * Get the grid function being used to provide initial conditions for a
     * quantity that has been registered with the hierarchy integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
    getInitialConditions(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Set an object to provide boundary conditions for a scalar-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoef(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                           SAMRAI::solv::RobinBcCoefStrategy<NDIM>* Q_bc_coef);

    /*!
     * Set objects to provide boundary conditions for a vector-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coef);

    /*!
     * Get objects used to provide boundary conditions for a scalar- or
     * vector-valued quantity that has been registered with the hierarchy
     * integrator.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>
    getPhysicalBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const;

    /*!
     * Register a solver for the Helmholtz equation (time-discretized diffusion
     * equation).
     */
    void setHelmholtzSolver(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                            SAMRAI::tbox::Pointer<IBTK::PoissonSolver> helmholtz_solver);

    /*!
     * Get the solver for the Helmholtz equation (time-discretized diffusion
     * equation) used by this solver class.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    getHelmholtzSolver(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Indicate that all of the Helmholtz solvers should be (re-)initialized before the
     * next time step.
     */
    void setHelmholtzSolversNeedInit();

    /*!
     * Indicate that the Helmholtz solver should be (re-)initialized before the
     * next time step.
     */
    void setHelmholtzSolverNeedsInit(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Register an operator to use to evaluate the right-hand side for the
     * Helmholtz solver (time-discretized diffusion equation).
     */
    void setHelmholtzRHSOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                 SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> helmholtz_rhs_operator);

    /*!
     * Get the operator to use to evaluate the right-hand side for the Helmholtz
     * solver (time-discretized diffusion equation).
     */
    SAMRAI::tbox::Pointer<IBTK::LaplaceOperator>
    getHelmholtzRHSOperator(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Indicate that all of the operators to evaluate the right-hand side for
     * the Helmholtz solver should be (re-)initialized before the next time
     * step.
     */
    void setHelmholtzRHSOperatorsNeedInit();

    /*!
     * Indicate that the operator to evaluate the right-hand side for the
     * Helmholtz solver should be (re-)initialized before the next time step.
     */
    void setHelmholtzRHSOperatorNeedsInit(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

protected:
    /*!
     * The constructor for class AdvDiffHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    AdvDiffHierarchyIntegrator(const std::string& object_name,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                               bool register_for_restart);

    /*!
     * Return the maximum stable time step size.
     */
    double getMaximumTimeStepSizeSpecialized();

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Standard variable registration.
     */
    void registerVariables();

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*!
     * Advective CFL condition.
     */
    double d_cfl_max;

    /*!
     * Default diffusion time integration method.
     */
    TimeSteppingType d_default_diffusion_time_stepping_type;

    /*!
     * Default convective differencing type.
     */
    ConvectiveDifferencingType d_default_convective_difference_form;

    /*!
     * Advection velocity data.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > > d_u_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> >, bool> d_u_is_div_free;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_u_fcn;

    /*!
     * Source term data.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_F_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_F_fcn;

    /*!
     * Diffusion coefficient data
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > > d_diffusion_coef_var,
        d_diffusion_coef_rhs_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_diffusion_coef_fcn;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > > d_diffusion_coef_rhs_map;

    /*!
     * Transported quantities.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_Q_var, d_Q_rhs_var;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > > d_Q_u_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_Q_F_map,
        d_Q_Q_rhs_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, TimeSteppingType>
        d_Q_diffusion_time_stepping_type;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, ConvectiveDifferencingType>
        d_Q_difference_form;

    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, double> d_Q_diffusion_coef;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > > d_Q_diffusion_coef_variable;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, bool> d_Q_is_diffusion_coef_variable;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, double> d_Q_damping_coef;

    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<IBTK::CartGridFunction> > d_Q_init;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> > d_Q_bc_coef;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;

    /*
     * Linear solvers and associated data.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_sol_vecs, d_rhs_vecs;
    std::string d_helmholtz_solver_type, d_helmholtz_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_helmholtz_solver_db, d_helmholtz_precond_db;
    std::vector<SAMRAI::tbox::Pointer<IBTK::PoissonSolver> > d_helmholtz_solvers;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LaplaceOperator> > d_helmholtz_rhs_ops;
    std::vector<bool> d_helmholtz_solvers_need_init, d_helmholtz_rhs_ops_need_init;
    int d_coarsest_reset_ln, d_finest_reset_ln;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffHierarchyIntegrator(const AdvDiffHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffHierarchyIntegrator& operator=(const AdvDiffHierarchyIntegrator& that);

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

#endif //#ifndef included_AdvDiffHierarchyIntegrator
