// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_StandardTagAndInitStrategySet
#define included_IBTK_StandardTagAndInitStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "BasePatchLevel.h"
#include "IntVector.h"
#include "StandardTagAndInitStrategy.h"
#include "tbox/Pointer.h"

#include <ostream>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class StandardTagAndInitStrategySet is a utility class that allows
 * multiple SAMRAI::mesh::StandardTagAndInitStrategy objects to be employed by a
 * single SAMRAI::mesh::StandardTagAndInitialize object.
 *
 * This class is primarily intended to provide a mechanism for allowing
 * user-define tagging criteria.
 *
 * \note When the return value of a method is a timestep size, class
 * StandardTagAndInitStrategySet returns the \em minimum timestep size computed
 * by the collection of concrete SAMRAI::mesh::StandardTagAndInitStrategy
 * objects.  Consequently, user-defined tagging classes \em must provide
 * implementations of getLevelDt() and advanceLevel() which return a large value
 * such as \a MAX_DOUBLE as the timestep size.
 *
 * \warning Not all concrete implementations of class
 * SAMRAI::mesh::StandardTagAndInitStrategy will work properly with class
 * StandardTagAndInitStrategySet.
 */
class StandardTagAndInitStrategySet : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    template <typename InputIterator>
    StandardTagAndInitStrategySet(InputIterator first, InputIterator last, bool managed = true)
        : d_strategy_set(first, last), d_managed(managed)
    {
        // intentionally blank
        return;
    } // StandardTagAndInitStrategySet

    /*!
     * \brief Destructor.
     */
    ~StandardTagAndInitStrategySet();

    /*!
     * Determine time increment to advance data on level.
     */
    double getLevelDt(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
                      double dt_time,
                      bool initial_time) override;

    /*!
     * Advance data on all patches on specified patch level from current time
     * (current_time) to new time (new_time).  This routine is called only
     * during time-dependent regridding procedures, such as Richardson
     * extrapolation.  It is virtual with an empty implementation here (rather
     * than pure virtual) so that users are not required to provide an
     * implementation when the function is not needed.  The boolean arguments
     * are used to determine the state of the algorithm and the data when the
     * advance routine is called.  Note that this advance function is also used
     * during normal time integration steps.
     *
     * When this function is called, the level data required to begin the
     * advance must be allocated and be defined appropriately.  Typically, this
     * is equivalent to what is needed to initialize a new level after
     * regridding.  Upon exiting this routine, both current and new data may
     * exist on the level.  This data is needed until level synchronization
     * occurs, in general.  Current and new data may be reset by calling the
     * member function resetTimeDependentData().
     *
     * This routine is called from two different points within the Richardson
     * extrapolation process: to advance a temporary level that is coarser than
     * the hierarchy level on which error estimation is performed, and to
     * advance the hierarchy level itself.  In the first case, the values of the
     * boolean flags are:
     *


     *    - \b  first_step
     *        = true.
     *    - \b  last_step
     *        = true.
     *    - \b  regrid_advance
     *        = true.
     *


     * In the second case, the values of the boolean flags are:
     *


     *    - \b  first_step
     *      (when regridding during time integration sequence)
     *        = true when the level is not coarsest level to synchronize
     *          immediately before the regridding process; else, false.
     *      (when generating initial hierarchy construction)
     *        = true, even though there may be multiple advance steps.
     *    - \b  last_step
     *        = true when the advance is the last in the Richardson
     *          extrapolation step sequence; else false.
     *    - \b  regrid_advance
     *        = true.
     *


     */
    double advanceLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                        double current_time,
                        double new_time,
                        bool first_step,
                        bool last_step,
                        bool regrid_advance = false) override;

    /*!
     * Reset time-dependent data storage for the specified patch level.
     */
    void resetTimeDependentData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
                                double new_time,
                                bool can_be_refined) override;

    /*!
     * Reset data on the patch level by destroying all patch data other than
     * that which is needed to initialize the solution on that level.  In other
     * words, this is the data needed to begin a time integration step on the
     * level.
     */
    void resetDataToPreadvanceState(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level) override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.  The level number indicates that of
     * the new level.
     *
     * Generally, when data is set, it is interpolated from coarser levels in
     * the hierarchy.  If the old level pointer in the argument list is
     * non-null, then data is copied from the old level to the new level on
     * regions of intersection between those levels before interpolation occurs.
     * In this case, the level number must match that of the old level.  The
     * specific operations that occur when initializing level data are
     * determined by the particular solution methods in use; i.e., in the
     * subclass of this abstract base class.
     *
     * The boolean argument initial_time indicates whether the level is being
     * introduced for the first time (i.e., at initialization time), or after
     * some regrid process during the calculation beyond the initial hierarchy
     * construction.  This information is provided since the initialization of
     * the data may be different in each of those circumstances.  The
     * can_be_refined boolean argument indicates whether the level is the finest
     * allowable level in the hierarchy.
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level =
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
                             bool allocate_data = true) override;

    /*!
     * After hierarchy levels have changed and data has been initialized on the
     * new levels, this routine can be used to reset any information needed by
     * the solution method that is particular to the hierarchy configuration.
     * For example, the solution procedure may cache communication schedules to
     * amortize the cost of data movement on the AMR patch hierarchy.  This
     * function will be called by the gridding algorithm after the
     * initialization occurs so that the algorithm-specific subclass can reset
     * such things.  Also, if the solution method must make the solution
     * consistent across multiple levels after the hierarchy is changed, this
     * process may be invoked by this routine.  Of course the details of these
     * processes are determined by the particular solution methods in use.
     *
     * The level number arguments indicate the coarsest and finest levels in the
     * current hierarchy configuration that have changed.  It should be assumed
     * that all intermediate levels have changed as well.
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to some user-supplied gradient criteria.  The
     * double time argument is the regrid time.  The integer "tag_index"
     * argument is the patch descriptor index of the cell-centered integer tag
     * array on each patch in the hierarchy.  The boolean argument initial_time
     * indicates whether the level is being subject to refinement at the initial
     * simulation time.  If it is false, then the error estimation process is
     * being invoked at some later time after the AMR hierarchy was initially
     * constructed.  Typically, this information is passed to the user's patch
     * tagging routines since the error estimator or gradient detector may be
     * different in each case.
     *
     * The boolean uses_richardson_extrapolation_too is true when Richardson
     * extrapolation error estimation is used in addition to the gradient
     * detector, and false otherwise.  This argument helps the user to manage
     * multiple regridding criteria.
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to some user-supplied Richardson extrapolation
     * criteria.  The "error_data_time" argument is the regrid time.  The
     * "deltat" argument is the time increment to advance the solution on the
     * level to be refined.  Note that that level is finer than the level in the
     * argument list, in general.  The ratio between the argument level and the
     * actual hierarchy level is given by the integer "coarsen ratio".
     *
     * The integer "tag_index" argument is the patch descriptor index of the
     * cell-centered integer tag array on each patch in the hierarchy.
     *
     * The boolean argument initial_time indicates whether the level is being
     * subject to refinement at the initial simulation time.  If it is false,
     * then the error estimation process is being invoked at some later time
     * after the AMR hierarchy was initially constructed.  Typically, this
     * information is passed to the user's patch tagging routines since the
     * application of the Richardson extrapolation process may be different in
     * each case.
     *
     * The boolean uses_gradient_detector_too is true when a gradient detector
     * procedure is used in addition to Richardson extrapolation, and false
     * otherwise.  This argument helps the user to manage multiple regridding
     * criteria.
     */
    void applyRichardsonExtrapolation(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                                      double error_data_time,
                                      int tag_index,
                                      double deltat,
                                      int error_coarsen_ratio,
                                      bool initial_time,
                                      bool uses_gradient_detector_too) override;

    /*!
     * Coarsen solution data from level to coarse_level for Richardson
     * extrapolation.  Note that this routine will be called twice during the
     * Richardson extrapolation error estimation process, once to set data on
     * the coarser level and once to coarsen data from after advancing the fine
     * level.  The init_coarse_level boolean argument indicates whether data is
     * set on the coarse level by coarsening the "old" time level solution or by
     * coarsening the "new" solution on the fine level (i.e., after it has been
     * advanced).
     */
    void coarsenDataForRichardsonExtrapolation(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                               int level_number,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level,
                                               double coarsen_data_time,
                                               bool before_advance) override;

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StandardTagAndInitStrategySet() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StandardTagAndInitStrategySet(const StandardTagAndInitStrategySet& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StandardTagAndInitStrategySet& operator=(const StandardTagAndInitStrategySet& that) = delete;

    /*!
     * \brief The set of SAMRAI::xfer:StandardTagAndInitStrategy objects.
     */
    std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> d_strategy_set;

    /*!
     * \brief Boolean value that indicates whether this class should provide
     * memory management for the strategy objects.
     */
    const bool d_managed;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_StandardTagAndInitStrategySet
