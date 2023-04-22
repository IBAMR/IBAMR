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

#ifndef included_IBAMR_IBExplicitHierarchyIntegrator
#define included_IBAMR_IBExplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBHierarchyIntegrator.h>

#include <ibtk/MarkerPatchHierarchy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <string>

namespace IBAMR
{
class IBStrategy;
class INSHierarchyIntegrator;
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
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBExplicitHierarchyIntegrator is an implementation of a formally
 * second-order accurate, semi-implicit version of the immersed boundary method.
 *
 * <h2>Working with marker points</h2>
 * - Specify the IB kernel in the input database as IB_delta_fcn
 * - Specify the output directory for H5Part data in the input database as
     viz_dump_dirname
 * - Set marker point positions with setMarkerPoints()
 */
class IBExplicitHierarchyIntegrator : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBExplicitHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    IBExplicitHierarchyIntegrator(std::string object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  SAMRAI::tbox::Pointer<IBStrategy> ib_method_ops,
                                  SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                  bool register_for_restart = true);

    /*!
     * The destructor for class IBExplicitHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~IBExplicitHierarchyIntegrator() = default;

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
     * Set the marker points stored by this class. This call is collective and
     * @p markers should be the same set of marker points on all processors.
     */
    virtual void setMarkers(const IBTK::EigenAlignedVector<IBTK::Point>& markers);

    /*!
     * Collect all marker point positions and velocities. This call is
     * collective.
     */
    virtual std::pair<IBTK::EigenAlignedVector<IBTK::Point>, IBTK::EigenAlignedVector<IBTK::Vector> >
    collectAllMarkers() const;

    /*!
     * Write marker plot data. This function does nothing if there are no markers.
     */
    virtual void writeMarkerPlotData(const int time_step,
                                     const double simulation_time = 0.0,
                                     const bool save_velocites = true) const;

protected:
    /*!
     * Perform necessary data movement, workload estimation, and logging prior
     * to regridding.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Perform necessary data movement and logging after regridding.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Velocity patch index at the halfway point with IB ghosting. Used only for
     * advecting marker points.
     */
    int d_u_half_idx = IBTK::invalid_index;

    /*!
     * Implementation of marker points.
     */
    SAMRAI::tbox::Pointer<IBTK::MarkerPatchHierarchy> d_markers;

    /*!
     * Boolean indicating whether or not we have set marker point velocities.
     * Since markers can be added at any time (including before the hierarchy is
     * set up) we must check at the beginning of each time step to see if we
     * need to call MarkerPatchHierarchy::setVelocities().
     */
    bool d_marker_velocities_set = false;

    /*!
     * IB kernel used for interpolation with d_marker_points.
     */
    std::string d_marker_kernel;

    /*!
     * Directory for marker point output data.
     */
    std::string d_viz_dump_dirname;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBExplicitHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBExplicitHierarchyIntegrator(const IBExplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBExplicitHierarchyIntegrator& operator=(const IBExplicitHierarchyIntegrator& that) = delete;

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Structure describing data temporarily stored by this class during a regrid.
     */
    struct RegridData
    {
        IBTK::EigenAlignedVector<IBTK::Point> d_marker_positions;

        IBTK::EigenAlignedVector<IBTK::Vector> d_marker_velocities;
    };

    /*!
     * Pointer to data only stored during regrids.
     */
    SAMRAI::tbox::Pointer<RegridData> d_regrid_temporary_data;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBExplicitHierarchyIntegrator
