#ifndef included_GodunovHypPatchOps
#define included_GodunovHypPatchOps

// Filename: GodunovHypPatchOps.h
// Last modified: <17.Aug.2006 15:54:51 boyce@bigboy.nyconnect.com>
// Created on 14 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/GodunovAdvector.h>
#include <ibamr/PhysicalBCDataStrategy.h>
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <FaceData.h>
#include <FaceVariable.h>
#include <GriddingAlgorithm.h>
#include <HyperbolicLevelIntegrator.h>
#include <HyperbolicPatchStrategy.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <VariableContext.h>
#include <VisItDataWriter.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * Class GodunovHypPatchOps provides numerical routines for solving
 * the conservative advection equation
 *
 *      dQ/dt + div(uQ) = Psi - Q div(u).
 * 
 * where Q is a cell centered quantity, u is a specified face centered
 * advection velocity, and Psi is an optional source term.  If u is
 * discretely divergence free, this is equivalent to
 *
 *      dQ/dt + div(uQ) = dQ/dt + u * grad(Q) = Psi.
 *
 * The class employs a predictor-corrector method to obtain a
 * (formally) second order accurate solution.  The predicted fluxes
 * are computed using the non-conservative advection equation, namely
 *
 *       dQ/dt + u * grad(Q) = Psi.
 *
 * Consequently, the source term Psi must be the source term
 * corresponding to the non-conservative advection equation.
 */
namespace IBAMR
{
class GodunovHypPatchOps
    : public SAMRAI::algs::HyperbolicPatchStrategy<NDIM>,
      public SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for GodunovHypPatchOps sets default parameters
     * for the advection solver.  The constructor also registers this
     * object for restart with the restart manager using the object
     * name.
     *
     * After default values are set, this routine calls
     * getFromRestart() if execution from a restart file is specified.
     * Finally, getFromInput() is called to read values from the given
     * input database (potentially overriding those found in the
     * restart file).
     */
    GodunovHypPatchOps(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<GodunovAdvector> godunov_advector,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
        bool register_for_restart=true);
    
    /*!
     * The destructor for GodunovHypPatchOps unregisters the patch
     * strategy object with the restart manager when so registered.
     */
    virtual ~GodunovHypPatchOps();
    
    ///
    ///  The following routines:
    ///
    ///      registerAdvectedQuantity(),
    ///      registerAdvectedQuantityWithSourceTerm(),
    ///      registerAdvectionVelocity(),
    ///      registerVisItDataWriter()
    ///
    ///  allow the GodunovHypPatchOps to be used as a generic
    ///  advection scheme.
    ///
    
    /*!
     * Register a cell centered quantity to be advected by the
     * GodunovAdvector according to the specified advection velocity.
     *
     * Conservative differencing is employed in updating the value of
     * the quantity when conservation_form is true.  Otherwise,
     * non-conservative differencing is used to update the quantity.
     *
     * Optional concrete SetDataStrategy and PhysicalBCDataStrategy
     * objects allow for the specification of initial and boundary
     * data for the advected quantity Q.  If an initialization object
     * is not specified, Q is initialized to zero.  If a boundary
     * condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.
     * (I.e. the domain can have no "physical" boundaries.)
     *
     * When the advected quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The
     * gradient is subtracted from the predicted face centered and
     * time centered values prior to the computation of the advective
     * fluxes.
     *
     * IMPORTANT NOTE: The advection velocity must be registered with
     * the patch strategy prior to the registration of advected
     * quantities.
     */
    void registerAdvectedQuantity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<SetDataStrategy> Q_init=NULL,
        SAMRAI::tbox::Pointer<PhysicalBCDataStrategy> Q_bc=NULL,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);

    /*!
     * Register a cell centered quantity to be advected by the
     * GodunovAdvector according to the specified advection velocity
     * and source term.
     *
     * Conservative differencing is employed in updating the value of
     * the quantity when conservation_form is true.  Otherwise,
     * non-conservative differencing is used to update the quantity.
     *
     * Optional concrete SetDataStrategy and PhysicalBCDataStrategy
     * objects allow for the specification of initial and boundary
     * data for the advected quantity Q.  If an initialization object
     * is not specified, Q is initialized to zero.  If a boundary
     * condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.
     * (I.e. the domain can have no "physical" boundaries.)
     *
     * The value of the source term is determined by an (optional)
     * SetDataStrategy object.  This allows for the specification of
     * either a constant or a time-dependent source term.  If this
     * object is not provided, the source term is initialized to zero.
     *
     * When the advected quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The
     * gradient is subtracted from the predicted face centered and
     * time centered values prior to the computation of the advective
     * fluxes.
     *
     * IMPORTANT NOTE: The advection velocity must be registered with
     * the patch strategy prior to the registration of advected
     * quantities.
     */
    void registerAdvectedQuantityWithSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<SetDataStrategy> Q_init=NULL,
        SAMRAI::tbox::Pointer<PhysicalBCDataStrategy> Q_bc=NULL,
        SAMRAI::tbox::Pointer<SetDataStrategy> Psi_set=NULL,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);
    
    /*!
     * Register a face centered advection velocity, used by the
     * GodunovAdvector to advect the cell centered quantities
     * registered with the patch strategy.
     *
     * An optional SetDataStrategy object allows for the specification
     * of a constant or time-dependent advection velocity.  If this
     * object is not provided, the advection velocity is initialized
     * to zero.
     *
     * The value of u_is_div_free determines whether the patch
     * strategy is able to assume that the discrete divergence of the
     * advection velocity is zero.
     *
     * IMPORTANT NOTE: The advection velocity must be registered with
     * the patch strategy prior to the registration of advected
     * quantities.
     */
    void registerAdvectionVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        const bool u_is_div_free,
        SAMRAI::tbox::Pointer<SetDataStrategy> u_set=NULL);

#if (NDIM>1)
    /*!
     * Register a VisIt data writer so this class will write plot
     * files that may be postprocessed with the VisIt visualization
     * tool.
     */
    void registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);
#endif

    ///
    ///  The following routines:
    ///
    ///      registerModelVariables(),
    ///      initializeDataOnPatch(),
    ///      computeStableDtOnPatch(),
    ///      computeFluxesOnPatch(),
    ///      conservativeDifferenceOnPatch(),
    ///      preprocessAdvanceLevelState(),
    ///      postprocessAdvanceLevelState(),
    ///      tagGradientDetectorCells(),
    ///      tagRichardsonExtrapolationCells()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::algs::HyperbolicPatchStrategy<NDIM> abstract base class.
    ///

    /*!
     * Register GodunovHypPatchOps model variables with the
     * SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> according to the variable
     * registration function provided by the integrator.
     *
     * In other words, variables are registered according to their
     * role in the integration process (e.g. time-dependent, flux,
     * etc.).  This routine also registers variables for plotting with
     * the VisIt writer.
     */
    virtual void registerModelVariables(
        SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>* integrator);
    
    /*!
     * Set the data on the patch interior to some initial values via
     * the concrete SetDataStrategy objects registered with the patch
     * strategy when provided.  Otherwise, initialize data to zero.
     */
    virtual void initializeDataOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time);
 
    /*!
     * Compute a stable time increment for patch using an explicit CFL
     * condition and return the computed dt.
     */
    virtual double computeStableDtOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const bool initial_time,
        const double dt_time);
 
    /*!
     * Compute the time integral of the fluxes to be used in
     * conservative difference for patch integration.
     *
     * The conservative difference used to update the integrated
     * quantities is implemented in conservativeDifferenceOnPatch().
     */
    virtual void computeFluxesOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double time,
        const double dt);
    
    /*!
     * Update solution variables by performing a conservative
     * difference using the fluxes calculated by
     * computeFluxesOnPatch().
     */
    virtual void conservativeDifferenceOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double time,
        const double dt,
        bool at_synchronization);
    
    /*!
     * Compute the values of any time-dependent source terms for use
     * by the explicit predictor.
     * 
     * This routine is called after patch boundary data is filled 
     * (i.e., ghosts) and before computeFluxesOnPatch().
     *
     * Note that when this routine is called, the scratch data is
     * filled on all patches (i.e., ghost cells) and that data is the
     * same as the current level data on all patch interiors.  That
     * is, both scratch and current data correspond to current_time.
     */
    virtual void preprocessAdvanceLevelState(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
        double current_time,
        double dt,
        bool first_step,
        bool last_step,
        bool regrid_advance);        
    
    /*!
     * Add source terms to the updated solution.
     * 
     * This routine is called after conservativeDifferenceOnPatch() is
     * called and before computeStableDtOnPatch().
     * 
     * Note that when this routine is called, the scratch data is
     * filled on all patches (i.e., ghost cells) and that data is the
     * same as the new level data on all patch interiors.  That is,
     * both scratch and new data correspond to current_time + dt on
     * patch interiors.  The current data and ghost values correspond
     * to the current_time.
     */
    virtual void postprocessAdvanceLevelState(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
        double current_time,
        double dt,
        bool first_step,
        bool last_step,
        bool regrid_advance);
    
    /*!
     * Tag cells for refinement using Richardson extrapolation.
     */
    virtual void tagRichardsonExtrapolationCells(
        SAMRAI::hier::Patch<NDIM>& patch,
        const int error_level_number,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> coarsened_fine,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> advanced_coarse,
        const double regrid_time,
        const double deltat,
        const int error_coarsen_ratio,
        const bool initial_error,
        const int tag_index,
        const bool uses_gradient_detector_too);
    
    /*!
     * Tag cells for refinement using a gradient detector.
     */
    virtual void tagGradientDetectorCells(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double regrid_time,
        const bool initial_error,
        const int tag_indexindx,
        const bool uses_richardson_extrapolation_too);

    ///
    ///  The following routines:
    ///
    ///      setPhysicalBoundaryConditions()
    ///
    ///  are concrete implementations of functions declared in the
    ///  xfer::RefinePatchStrategy abstract base class.
    ///

    /*!
     * Set the data in ghost cells corresponding to physical boundary
     * conditions.
     */ 
    virtual void setPhysicalBoundaryConditions(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);
    
    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///
    
    /*!
     * Write state of GodunovHypPatchOps object to the given database
     * for restart.
     *
     * This routine is a concrete implementation of the function
     * declared in the SAMRAI::tbox::Serializable abstract base class.
     */
    virtual void putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    ///
    ///  The following routines:
    ///
    ///      printClassData()
    ///
    ///  are provided for your viewing pleasure.
    ///
    
    /*!
     * Print all data members for GodunovHypPatchOps class.
     */
    virtual void printClassData(
        ostream& os) const;

protected:
    /*
     * Get a pointer to the requested flux integral patch data on the
     * specified patch.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > getFluxIntegralData(
        const size_t l,
        SAMRAI::hier::Patch<NDIM>& patch,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context);
    
    /*
     * Get a pointer to the requested q integral patch data on the
     * specified patch.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > getQIntegralData(
        const size_t l,
        SAMRAI::hier::Patch<NDIM>& patch,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context);
    
    /*
     * The SAMRAI::algs::HyperbolicLevelIntegrator that is using the
     * patch strategy.
     */
    SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>* d_integrator;
    
    /*
     * The GodunovAdvector being used to advect the cell centered
     * quantities Q.
     */
    SAMRAI::tbox::Pointer<GodunovAdvector> d_godunov_advector;
    
    /*! 
     * Advected quantities Q, source terms Psi (possibly NULL) and the
     * optional face centered gradient terms used to enforce
     * incompressibility.
     
    */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Psi_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_grad_vars;

    /*!
     * Indicates whether conservative or non-conservative differencing
     * should be employed for a given quantity.
     */
    std::vector<bool> d_Q_conservation_form;
    
    /*!
     * When conservative differencing is employed for a quantity Q, we
     * maintain the time integral of the advective flux corresponding
     * to that quantity.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_flux_integral_vars;
    
    /*!
     * When non-conservative differencing is employed for a quantity
     * Q, we maintain the time integral of the predicted value and the
     * advection velocity.
     *
     * These values must also be maintained when the advection
     * velocity is not discretely divergence free.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_q_integral_vars;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_u_integral_var;

    /*!
     * Objects to set initial and boundary conditions as well as
     * forcing terms for each advected quantity.
     */    
    std::vector<SAMRAI::tbox::Pointer<SetDataStrategy> >        d_Q_inits;
    std::vector<SAMRAI::tbox::Pointer<PhysicalBCDataStrategy> > d_Q_bcs;
    
    std::vector<SAMRAI::tbox::Pointer<SetDataStrategy> >  d_Psi_sets;
    
    /*!
     * The advection velocity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_u_var;
    SAMRAI::tbox::Pointer<SetDataStrategy> d_u_set;
    bool d_u_is_div_free;
    bool d_compute_init_velocity, d_compute_half_velocity, d_compute_final_velocity;

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    GodunovHypPatchOps();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    GodunovHypPatchOps(
        const GodunovHypPatchOps& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    GodunovHypPatchOps& operator=(
        const GodunovHypPatchOps& that);

    /*
     * These private member functions read data from input and
     * restart.  When beginning a run from a restart file, all data
     * members are read from the restart file.  If the boolean flag is
     * true when reading from input, some restart values may be
     * overridden by those in the input file.
     *
     * An assertion results if the database pointer is null.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);
    
    void getFromRestart();
    
    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.  The boolean is
     * used to control restart file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * We cache pointers to the grid geometry and VisIt data writer
     * object to set up initial data, set physical boundary
     * conditions, and register plot variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
#if (NDIM>1)
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;
#endif
    
    /*
     *  Parameters for numerical method:
     *
     *    d_ghosts .............. number of ghost cells for cell centered
     *                            and face/side centered variables
     *    d_flux_ghosts ......... number of ghost cells for fluxes
     *    d_u_is_registered ..... indicates that an advection velocity
     *                            has been registered with the patch
     *                            strategy
     */
    SAMRAI::hier::IntVector<NDIM> d_ghosts;
    SAMRAI::hier::IntVector<NDIM> d_flux_ghosts;
    bool d_u_is_registered;
    
    /*
     * Refinement criteria parameters for gradient detection and
     * Richardson extrapolation.
     */
    SAMRAI::tbox::Array<std::string> d_refinement_criteria;
    SAMRAI::tbox::Array<double> d_dev_tol;
    SAMRAI::tbox::Array<double> d_dev;
    SAMRAI::tbox::Array<double> d_dev_time_max;
    SAMRAI::tbox::Array<double> d_dev_time_min;
    SAMRAI::tbox::Array<double> d_grad_tol;
    SAMRAI::tbox::Array<double> d_grad_time_max;
    SAMRAI::tbox::Array<double> d_grad_time_min;
    SAMRAI::tbox::Array<double> d_rich_tol;
    SAMRAI::tbox::Array<double> d_rich_time_max;
    SAMRAI::tbox::Array<double> d_rich_time_min;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifndef DEBUG_NO_INLINE
//#include "GodunovHypPatchOps.I"
//#endif

#endif //#ifndef included_GodunovHypPatchOps
