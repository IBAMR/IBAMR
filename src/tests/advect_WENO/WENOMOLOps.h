#ifndef included_QInit
#define included_QInit

// Filename: QInit.h
// Last modified: <04.Jan.2008 23:16:46 griffith@box221.cims.nyu.edu>
// Created on 04 Jan 2008 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <FaceVariable.h>
#include <GridGeometry.h>
#include <MethodOfLinesPatchStrategy.h>
#include <VisItDataWriter.h>
#include <tbox/Array.h>
#include <tbox/Database.h>

// NAMESPACE
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class to implement the WENO advection algorithm with user-specified
 * initial conditions and velocity profiles.
 */
class WENOMOLOps
    : public algs::MethodOfLinesPatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    WENOMOLOps(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~WENOMOLOps();

    /*!
     * Register a VisIt data writer so this class will write
     * plot files that may be postprocessed with the VisIt
     * visualization tool.
     */
    void
    registerVisItDataWriter(
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * Register variables specific to the problem to be solved with the
     * integrator using the registerVariable function.
     */
    virtual void
    registerModelVariables(
        algs::MethodOfLinesIntegrator<NDIM>* integrator);

    /*!
     * Set the initial data on a patch interior (i.e., NO GHOST CELLS).
     */
    virtual void
    initializeDataOnPatch(
        hier::Patch<NDIM>& patch,
        const double time,
        const bool initial_time) const;

    /*!
     * Compute the stable time increment for a patch.
     */
    virtual double
    computeStableDtOnPatch(
        hier::Patch<NDIM>& patch,
        const double time) const;

    /*!
     * Advance a single Runge Kutta step.
     *
     * \param patch patch that RK step is being applied
     * \param dt    timestep
     * \param alpha_1 first coefficient applied in the RK step
     * \param alpha_2 second coefficient
     * \param beta    third coefficient
     */
    virtual void
    singleStep(
        hier::Patch<NDIM>& patch,
        const double dt,
        const double alpha_1,
        const double alpha_2,
        const double beta) const;

    /*!
     * Set user-defined boundary conditions at the physical domain boundary.
     */
    virtual void
    setPhysicalBoundaryConditions(
        hier::Patch<NDIM>& patch,
        const double fill_time,
        const hier::IntVector<NDIM>& ghost_width_to_fill);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    WENOMOLOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    WENOMOLOps(
        const WENOMOLOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    WENOMOLOps&
    operator=(
        const WENOMOLOps& that);

    /*!
     * Read input values, indicated above, from given database.
     */
    void
    getFromInput(
        tbox::Pointer<tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The VisIt writer.
     */
    tbox::Pointer<appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*
     * The grid geometry.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * Variables.
     */
    tbox::Pointer<pdat::CellVariable<NDIM,double> > d_Q_var, d_U_var, d_F_var;
    tbox::Pointer<pdat::FaceVariable<NDIM,double> > d_flux_var;

    /*
     * The center of the initial data.
     */
    tbox::Array<double> d_X;

    /*
     * The initialization type.
     */
    string d_init_type;

    /*
     * Parameters for Gaussian initial conditions.
     */
    double d_gaussian_kappa;

    /*
     * Parameters for the Zalesak slotted cylinder.
     */
    double d_zalesak_r;
    double d_zalesak_slot_w;
    double d_zalesak_slot_l;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "WENOMOLOps.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_WENOMOLOps
