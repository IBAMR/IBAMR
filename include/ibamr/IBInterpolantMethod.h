// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_IBInterpolantMethod
#define included_IBAMR_IBInterpolantMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBStrategy.h"

#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core"
#include "Eigen/Geometry"
IBTK_ENABLE_EXTRA_WARNINGS

#include <set>
#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
class HierarchyIntegrator;
} // namespace IBTK

namespace IBTK
{
class LData;
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace tbox
{
class Database;
template <class TYPE>
class Array;
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
 * \brief Class IBInterpolantMethod is an implementation of the abstract base class
 * IBStrategy that provides a simple functionality of interpolating and spreading scalar
 * or vector-valued quantities to and from the moving IB structure on the background Cartesian
 * grid.
 *
 * This class does not provide direct support for fluid-structure interaction, but rather
 * facilitates Eulerian-Lagrangian data management for LMesh.
 */
class IBInterpolantMethod : public IBStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBInterpolantMethod(std::string object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        int no_structures = 1,
                        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~IBInterpolantMethod();

    /*!
     * Register Eulerian variables with the parent IBHierarchyIntegrator with
     * the VariableDatabase. They are used for interpolating the Eulerian quantities
     * on the Lagrangian mesh.
     */
    void registerEulerianVariables() override;

    /*!
     * Register ghost cell filling algrorithms for the Eulerian variables.
     */
    void registerEulerianCommunicationAlgorithms() override;

    /*!
     * \brief Register a variable and the HierarchyIntegrator managing the variable
     * with this class.
     */
    virtual void
    registerVariableAndHierarchyIntegrator(const std::string& var_name,
                                           const int var_depth,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> hier_integrator);

    /*!
     * Supply a Lagrangian initialization object.
     */
    void registerLInitStrategy(SAMRAI::tbox::Pointer<IBTK::LInitStrategy> l_initializer);

    /*!
     * Free references to Lagrangian initialization objects.
     */
    void freeLInitStrategy();

    /*!
     * Return a pointer to the Lagrangian data manager object.
     */
    IBTK::LDataManager* getLDataManager() const;

    /*!
     * \brief Get the level on which the structures reside.
     */
    int getStructuresLevelNumber() const;

    /*!
     * \brief Get the structure handle to which this Lagrangian index belongs.
     */
    int getStructureHandle(const int lag_idx) const;

    /*!
     * \brief Get the number of nodes for this structure.
     */
    unsigned int getNumberOfNodes(const unsigned int struct_no) const
    {
        std::pair<int, int> lag_idx_range = d_struct_lag_idx_range[struct_no];
        return (lag_idx_range.second - lag_idx_range.first);

    } // getNumberOfStructuresNodes

    /*!
     * Register a Lagrangian Silo data writer so this class will write plot
     * files that may be postprocessed with the VisIt visualization tool.
     */
    void registerLSiloDataWriter(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer);

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>& getMinimumGhostCellWidth() const override;

    /*!
     * Setup the tag buffer.
     */
    void setupTagBuffer(SAMRAI::tbox::Array<int>& tag_buffer,
                        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const override;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) override;

    /*!
     * Interpolate the Eulerian quantities to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateQ(double data_time);

    /*!
     * Interpolate the Eulerian quantities to the curvilinear mesh at both current and new time.
     */
    void interpolateQ();

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the backward Euler
     * method.
     */
    void backwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the trapezoidal
     * rule.
     */
    void trapezoidalStep(double current_time, double new_time) override;

    /*!
     * Update the position of the mesh at new time.
     *
     * X^n+1 = R[theta]*X^n + U*dt
     */
    void updateMeshPosition(double current_time,
                            double new_time,
                            const IBTK::EigenAlignedVector<Eigen::Vector3d>& U,
                            const IBTK::EigenAlignedVector<Eigen::Vector3d>& W);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time) override;

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * Spread Q on the Eulerian grid.
     */
    void spreadQ(double data_time);

    /*!
     * Get the patch data index of Q.
     */
    int getEulerianQPatchDataIndex(const std::string& var_name)
    {
        return d_q_interp_idx[var_name];
    } // getEulerianQPatchDataIndex

    /*!
     * Copy the spreaded data back to the integrator.
     */
    void copyEulerianDataToIntegrator(double data_time);

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
        bool initial_time) override;

    /*!
     * Add the estimated computational work from the current object per cell
     * into the specified <code>workload_data_idx</code>.
     */
    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const int workload_data_idx) override;

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void beginDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void endDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

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
                             bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too) override;

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * Get the structure position data.
     */
    void getPositionData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** X_data, double data_time);

    /*!
     * Get the q data defined on Eulerian grid for interpolation.
     */
    void copyEulerianDataFromIntegrator(const std::string& var_name, int q_data_idx, double data_time);

    /*!
     * Get the q data defined on Eulerian grid for spreading.
     */
    void zeroOutEulerianData(const std::string& var_name, int q_data_idx);

    /*!
     * Copy the q_data_idx to suitable intergrator.
     */
    void copyEulerianDataToIntegrator(const std::string& name, int q_data_idx, double data_time);

    /*!
     * Get the Q data defined on Lagrangian mesh.
     */
    void
    getQData(const std::string& var_name, std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** Q_data, double data_time);

    /*!
     * \brief Compute the center of mass of the mesh.
     */
    void computeCenterOfMass(IBTK::EigenAlignedVector<Eigen::Vector3d>& center_of_mass,
                             std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data);

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log = false;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Vector of Lagrangian indices of all structures.
     */
    std::vector<std::pair<int, int> > d_struct_lag_idx_range;

    // Number of rigid parts.
    unsigned int d_num_rigid_parts;

    /*
     * Rigid body motion of the Lagrangian mesh.
     */
    IBTK::EigenAlignedVector<Eigen::Vector3d> d_center_of_mass_initial, d_center_of_mass_current, d_center_of_mass_new;
    IBTK::EigenAlignedVector<Eigen::Quaterniond> d_quaternion_current, d_quaternion_new;

    /*
     * The LDataManager is used to coordinate the distribution of Lagrangian
     * data on the patch hierarchy.
     */
    IBTK::LDataManager* d_l_data_manager;
    std::string d_interp_kernel_fcn = "IB_4", d_spread_kernel_fcn = "IB_4";
    bool d_error_if_points_leave_domain = false;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*!
     * \brief Eulerian variables and their HierarchyIntegrators.
     */
    std::map<std::string, SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_q_var;
    std::map<std::string, int> d_q_depth;
    std::map<std::string, int> d_q_interp_idx;
    std::map<std::string, SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> > d_q_hier_integrator;

    /*
     * Lagrangian variables.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_X_current_data, d_X_new_data;
    std::map<std::string, std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > > d_Q_current_data, d_Q_new_data;

    /*
     * The specification and initialization information for the Lagrangian data
     * used by the integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::LInitStrategy> d_l_initializer;

    /*
     * Visualization data writers.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBInterpolantMethod() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInterpolantMethod(const IBInterpolantMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInterpolantMethod& operator=(const IBInterpolantMethod& that) = delete;

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

#endif //#ifndef included_IBAMR_IBInterpolantMethod
