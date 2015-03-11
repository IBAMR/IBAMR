// Filename: AdvectorPredictorCorrectorHyperbolicPatchOps.h
// Created on 14 Feb 2004 by Boyce Griffith
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

#ifndef included_AdvectorPredictorCorrectorHyperbolicPatchOps
#define included_AdvectorPredictorCorrectorHyperbolicPatchOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/algs/HyperbolicPatchStrategy.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "SAMRAI/tbox/Array.h"

#include "SAMRAI/tbox/Serializable.h"

namespace IBTK
{
class CartGridFunction;
} // namespace IBTK
namespace SAMRAI
{
namespace algs
{
class HyperbolicLevelIntegrator;
} // namespace algs
namespace hier
{
class VariableContext;
class Patch;
class PatchLevel;
} // namespace hier
namespace pdat
{
template <class TYPE>
class CellVariable;
template <class TYPE>
class FaceData;
template <class TYPE>
class FaceVariable;
} // namespace pdat
namespace solv
{
class RobinBcCoefStrategy;
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
 * \brief Class AdvectorPredictorCorrectorHyperbolicPatchOps is a concrete
 * SAMRAI::algs::HyperbolicPatchStrategy that makes use of class
 *AdvectorExplicitPredictorPatchOps
 * to solve the linear advection equation.
 *
 * Class AdvectorPredictorCorrectorHyperbolicPatchOps provides numerical routines for solving
 *the
 *advection
 * equation in conservative form, \f[
 *
 *      \frac{dQ}{dt} + \nabla \cdot (\vec{u}^{\mbox{\scriptsize ADV}} Q) = F - Q \nabla \cdot
 *\vec{u}^{\mbox{\scriptsize ADV}},
 *
 * \f] where \f$ Q \f$ is a cell-centered scalar- or vector-valued quantity, \f$
 * \vec{u}^{\mbox{\scriptsize ADV}} \f$ is a specified face-centered advection
 * velocity, and \f$ F \f$ is an optional source term.  When \f$
 * \vec{u}^{\mbox{\scriptsize ADV}} \f$ is discretely divergence free, this is
 * equivalent to solving\f[
 *
 *      \frac{dQ}{dt} + \nabla \cdot (\vec{u}^{\mbox{\scriptsize ADV}} Q) = F.
 *
 * \f] The class employs a predictor-corrector method to obtain a (formally)
 * second-order accurate solution.  The predicted fluxes are computed using the
 * \em non-conservative advection equation, namely \f[
 *
 *      \frac{dQ}{dt} + (\vec{u}^{\mbox{\scriptsize ADV}} \cdot \nabla) Q = F.
 *
 * \f] Consequently, if the form of the source term \f$ F \f$ depends on whether
 * \f$ Q \f$ is being conservatively or non-conservatively differenced, it is
 * \em crucial that \f$ F \f$ correspond to the \em non-conservative form of the
 * advection equation.
 *
 * This class can also be used to solve the \em non-conservative form of the
 * advection equation, i.e., \f[
 *
 *      \frac{dQ}{dt} + (\vec{u}^{\mbox{\scriptsize ADV}} \cdot \nabla) Q = F.
 *
 * \f]
 */
class AdvectorPredictorCorrectorHyperbolicPatchOps : public SAMRAI::algs::HyperbolicPatchStrategy,
                                                     public SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for AdvectorPredictorCorrectorHyperbolicPatchOps sets default parameters
     *for
     * the advection solver.  The constructor also registers this object for
     * restart with the restart manager using the object name.
     *
     * After default values are set, this routine calls getFromRestart() if
     * execution from a restart file is specified.  Finally, getFromInput() is
     * called to read values from the given input database (potentially
     * overriding those found in the restart file).
     */
    AdvectorPredictorCorrectorHyperbolicPatchOps(
        const std::string& object_name,
        boost::shared_ptr<SAMRAI::tbox::Database> input_db,
        boost::shared_ptr<AdvectorExplicitPredictorPatchOps> explicit_predictor,
        boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geom,
        bool register_for_restart = true);

    /*!
     * The destructor for AdvectorPredictorCorrectorHyperbolicPatchOps unregisters the patch
     * strategy object with the restart manager when so registered.
     */
    ~AdvectorPredictorCorrectorHyperbolicPatchOps();

    /*!
     * Return the name of the patch operations object.
     */
    const std::string& getName() const;

    /*!
     * Register a VisIt data writer so this class will write plot files that may
     * be postprocessed with the VisIt visualization tool.
     */
    void registerVisItDataWriter(boost::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_writer);

    /*!
     * Register a face-centered advection velocity to be used to advect
     * cell-centered quantities by the hierarchy integrator.
     *
     * \note By default, each registered advection velocity is assumed to be
     * divergence free.
     */
    void registerAdvectionVelocity(boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > u_var);

    /*!
     * Indicate whether a particular advection velocity is discretely divergence
     * free.
     */
    void setAdvectionVelocityIsDivergenceFree(boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > u_var,
                                              bool is_div_free);

    /*!
     * Set an IBTK::CartGridFunction object that specifies the value of a
     * particular advection velocity.
     */
    void setAdvectionVelocityFunction(boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > u_var,
                                      boost::shared_ptr<IBTK::CartGridFunction> u_fcn);

    /*!
     * Register a cell-centered source term.
     */
    void registerSourceTerm(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > F_var);

    /*!
     * Set an IBTK::CartGridFunction object that specifies the value of a
     * particular source term.
     */
    void setSourceTermFunction(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > F_var,
                               boost::shared_ptr<IBTK::CartGridFunction> F_fcn);

    /*!
     * Register a cell-centered quantity to be advected and diffused by the
     * hierarchy integrator.
     */
    void registerTransportedQuantity(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var);

    /*!
     * Set the face-centered advection velocity to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified advection velocity must have been already registered
     * with the hierarchy integrator.
     */
    void setAdvectionVelocity(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > u_var);

    /*!
     * Set the cell-centered source term to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void setSourceTerm(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                       boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > F_var);

    /*!
     * Set the convective differencing form for a quantity that has been
     * registered with the hierarchy integrator.
     */
    void setConvectiveDifferencingType(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                                       ConvectiveDifferencingType difference_form);

    /*!
     * Set a grid function to provide initial conditions for a quantity that has
     * been registered with the hierarchy integrator.
     */
    void setInitialConditions(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                              boost::shared_ptr<IBTK::CartGridFunction> Q_init);

    /*!
     * Set an object to provide boundary conditions for a scalar-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoefs(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                            SAMRAI::solv::RobinBcCoefStrategy* Q_bc_coef);

    /*!
     * Set objects to provide boundary conditions for a vector-valued quantity
     * that has been registered with the hierarchy integrator.
     */
    void setPhysicalBcCoefs(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                            std::vector<SAMRAI::solv::RobinBcCoefStrategy*> Q_bc_coef);

    /*!
     * \brief Register AdvectorPredictorCorrectorHyperbolicPatchOps model variables with the
     * SAMRAI::algs::HyperbolicLevelIntegrator according to the variable
     * registration function provided by the integrator.
     *
     * In other words, variables are registered according to their role in the
     * integration process (e.g. time-dependent, flux, etc.).  This routine also
     * registers variables for plotting with the VisIt writer.
     */
    void registerModelVariables(SAMRAI::algs::HyperbolicLevelIntegrator* integrator);

    /*!
     * \brief Set the data on the patch interior to some initial values via the
     * concrete IBTK::CartGridFunction objects registered with the patch strategy when
     * provided.  Otherwise, initialize data to zero.
     */
    void initializeDataOnPatch(SAMRAI::hier::Patch& patch, double data_time, bool initial_time);

    /*!
     * \brief Compute a stable time increment for patch using an explicit CFL
     * condition and return the computed dt.
     */
    double computeStableDtOnPatch(SAMRAI::hier::Patch& patch, bool initial_time, double dt_time);

    /*!
     * \brief Compute the time integral of the fluxes to be used in conservative
     * difference for patch integration.
     *
     * The conservative difference used to update the integrated quantities is
     * implemented in conservativeDifferenceOnPatch().
     */
    void computeFluxesOnPatch(SAMRAI::hier::Patch& patch, double time, double dt);

    /*!
     * \brief Update solution variables by performing a conservative difference
     * using the fluxes calculated by computeFluxesOnPatch().
     */
    void conservativeDifferenceOnPatch(SAMRAI::hier::Patch& patch, double time, double dt, bool at_synchronization);

    /*!
     * \brief Compute the values of any time-dependent source terms for use by
     * the explicit predictor.
     *
     * This routine is called \em after patch boundary data is filled (i.e.,
     * ghosts) and \em before computeFluxesOnPatch().
     *
     * Note that when this routine is called, the scratch data is filled on all
     * patches (i.e., ghost cells) and that data is the same as the current
     * level data on all patch interiors.  That is, both scratch and current
     * data correspond to current_time.
     */
    void preprocessAdvanceLevelState(const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                     double current_time,
                                     double dt,
                                     bool first_step,
                                     bool last_step,
                                     bool regrid_advance);

    /*!
     * \brief Add source terms to the updated solution.
     *
     * This routine is called \em after conservativeDifferenceOnPatch() is
     * called and \em before computeStableDtOnPatch().
     *
     * Note that when this routine is called, the scratch data is filled on all
     * patches (i.e., ghost cells) and that data is the same as the new level
     * data on all patch interiors.  That is, both scratch and new data
     * correspond to current_time + dt on patch interiors.  The current data and
     * ghost values correspond to the current_time.
     */
    void postprocessAdvanceLevelState(const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                      double current_time,
                                      double dt,
                                      bool first_step,
                                      bool last_step,
                                      bool regrid_advance);

    /*!
     * \brief Tag cells for refinement using a gradient detector.
     */
    void tagGradientDetectorCells(SAMRAI::hier::Patch& patch,
                                  double regrid_time,
                                  bool initial_error,
                                  int tag_indexindx,
                                  bool uses_richardson_extrapolation_too);

    /*!
     * \brief Set the data in ghost cells corresponding to physical boundary
     * conditions.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector& ghost_width_to_fill);

    /*!
     * \brief Return maximum stencil width needed for user-defined
     * data refinement operations performed by this object.
     */
    SAMRAI::hier::IntVector getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const;

    /*!
     * \brief Perform user-defined patch data refinement operations.
     */
    void preprocessRefine(SAMRAI::hier::Patch& fine,
                          const SAMRAI::hier::Patch& coarse,
                          const SAMRAI::hier::Box& fine_box,
                          const SAMRAI::hier::IntVector& ratio);

    /*!
     * \brief Perform user-defined patch data refinement operations.
     */
    void postprocessRefine(SAMRAI::hier::Patch& fine,
                           const SAMRAI::hier::Patch& coarse,
                           const SAMRAI::hier::Box& fine_box,
                           const SAMRAI::hier::IntVector& ratio);
    /*!
     * \brief Return maximum stencil width needed for user-defined
     * data coarsening operations performed by this object.
     */
    SAMRAI::hier::IntVector getCoarsenOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const;

    /*!
     * \brief Perform user-defined patch data coarsening operations.
     */
    void preprocessCoarsen(SAMRAI::hier::Patch& coarse,
                           const SAMRAI::hier::Patch& fine,
                           const SAMRAI::hier::Box& coarse_box,
                           const SAMRAI::hier::IntVector& ratio);

    /*!
     * \brief Perform user-defined patch data coarsening operations.
     */
    void postprocessCoarsen(SAMRAI::hier::Patch& coarse,
                            const SAMRAI::hier::Patch& fine,
                            const SAMRAI::hier::Box& coarse_box,
                            const SAMRAI::hier::IntVector& ratio);

    /*!
     * \brief Write state of AdvectorPredictorCorrectorHyperbolicPatchOps object to the given
     * database for restart.
     *
     * This routine is a concrete implementation of the function declared in the
     * SAMRAI::tbox::Serializable abstract base class.
     */
    void putToRestart(const boost::shared_ptr<SAMRAI::tbox::Database>& restart_db) const;

protected:
    /*!
     * \brief Get a pointer to the requested flux integral patch data on the
     * specified patch.
     */
    boost::shared_ptr<SAMRAI::pdat::FaceData<double> >
    getFluxIntegralData(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                        SAMRAI::hier::Patch& patch,
                        boost::shared_ptr<SAMRAI::hier::VariableContext> context);

    /*!
     * \brief Get a pointer to the requested q integral patch data on the
     * specified patch.
     */
    boost::shared_ptr<SAMRAI::pdat::FaceData<double> >
    getQIntegralData(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                     SAMRAI::hier::Patch& patch,
                     boost::shared_ptr<SAMRAI::hier::VariableContext> context);

    /*!
     * \brief Get a pointer to the requested u integral patch data on the
     * specified patch.
     */
    boost::shared_ptr<SAMRAI::pdat::FaceData<double> >
    getUIntegralData(boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > Q_var,
                     SAMRAI::hier::Patch& patch,
                     boost::shared_ptr<SAMRAI::hier::VariableContext> context);

    /*
     * The SAMRAI::algs::HyperbolicLevelIntegrator that is using the patch
     * strategy.
     */
    SAMRAI::algs::HyperbolicLevelIntegrator* d_integrator;

    /*
     * The AdvectorExplicitPredictorPatchOps being used to advect the cell-centered quantities
     * Q.
     */
    boost::shared_ptr<AdvectorExplicitPredictorPatchOps> d_explicit_predictor;

    /*
     * Advection velocity data.
     */
    std::set<boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > > d_u_var;
    std::map<boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> >, bool> d_u_is_div_free;
    std::map<boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> >, boost::shared_ptr<IBTK::CartGridFunction> >
        d_u_fcn;
    bool d_compute_init_velocity, d_compute_half_velocity, d_compute_final_velocity;

    /*
     * Source term data.
     */
    std::set<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > > d_F_var;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >, boost::shared_ptr<IBTK::CartGridFunction> >
        d_F_fcn;

    /*
     * Transported quantities.
     */
    std::set<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > > d_Q_var;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >,
             boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > > d_Q_u_map;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > > d_Q_F_map;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >, ConvectiveDifferencingType> d_Q_difference_form;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >, boost::shared_ptr<IBTK::CartGridFunction> >
        d_Q_init;
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >, std::vector<SAMRAI::solv::RobinBcCoefStrategy*> >
        d_Q_bc_coef;

    /*
     * When conservative differencing is employed for a quantity Q, we maintain
     * the time integral of the advective flux corresponding to that quantity.
     */
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >,
             boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > > d_flux_integral_var;

    /*
     * When non-conservative differencing is employed for a quantity Q, we
     * maintain the time integral of the predicted value and the advection
     * velocity.
     *
     * These values must also be maintained in the case in which the advection
     * velocity is not discretely divergence free.
     */
    std::map<boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >,
             boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > > d_q_integral_var;
    std::map<boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> >,
             boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > > d_u_integral_var;

    /*
     * Boolean indicating whether or not to overwrite tag data (default is
     * true).
     */
    bool d_overwrite_tags;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvectorPredictorCorrectorHyperbolicPatchOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvectorPredictorCorrectorHyperbolicPatchOps(const AdvectorPredictorCorrectorHyperbolicPatchOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvectorPredictorCorrectorHyperbolicPatchOps& operator=(const AdvectorPredictorCorrectorHyperbolicPatchOps& that);

    /*
     * Set physical boundary conditions at inflow boundaries for predicted
     * face-centered quantities.
     */
    void setInflowBoundaryConditions(SAMRAI::hier::Patch& patch, double fill_time);

    /*
     * These private member functions read data from input and restart.  When
     * beginning a run from a restart file, all data members are read from the
     * restart file.  If the boolean flag is true when reading from input, some
     * restart values may be overridden by those in the input file.
     *
     * An assertion results if the database pointer is null.
     */
    void getFromInput(boost::shared_ptr<SAMRAI::tbox::Database> db, bool is_from_restart);
    void getFromRestart();

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * We cache pointers to the grid geometry and VisIt data writer object to
     * set up initial data, set physical boundary conditions, and register plot
     * variables.
     */
    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geometry;
    boost::shared_ptr<SAMRAI::appu::VisItDataWriter> d_visit_writer;

    /*
     * Boundary condition extrapolation helpers.
     */
    IBTK::CartExtrapPhysBdryOp d_extrap_bc_helper;

    /*
     *  Parameters for numerical method:
     *
     *    d_ghosts .............. number of ghost cells for cell-centered and
     *                            face/side-centered variables
     *    d_flux_ghosts ......... number of ghost cells for fluxes
     *
     *    d_extrap_type ......... type of extrapolation to use at outflow
     *                            boundaries (choices are: CONSTANT, LINEAR)
     */
    SAMRAI::hier::IntVector d_ghosts;
    SAMRAI::hier::IntVector d_flux_ghosts;
    std::string d_extrap_type;

    /*
     * Refinement criteria parameters for gradient detection.
     */
    std::vector<std::string> d_refinement_criteria;
    std::vector<double> d_dev_tol;
    std::vector<double> d_dev;
    std::vector<double> d_dev_time_max;
    std::vector<double> d_dev_time_min;
    std::vector<double> d_grad_tol;
    std::vector<double> d_grad_time_max;
    std::vector<double> d_grad_time_min;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvectorPredictorCorrectorHyperbolicPatchOps
