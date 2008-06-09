// Filename: IBHierarchyIntegrator.C
// Last modified: <08.Jun.2008 20:45:06 griffith@box230.cims.nyu.edu>
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBInstrumentationSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeIndexData2.h>
#include <ibtk/LagMarkerCoarsenOperator.h>
#include <ibtk/LagMarkerRefineOperator.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellIndex.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <IndexData.h>
#include <Patch.h>
#include <PatchCellDataOpsReal.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_integrator;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_advance_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Version of IBHierarchyIntegrator restart file data.
static const int IB_HIERARCHY_INTEGRATOR_VERSION = 1;

inline std::string
discard_comments(
    const std::string& input_string)
{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();

    return output_string;
}// discard_comments

inline double
cos_delta(
    const double x,
    const double eps)
{
    if (fabs(x) > eps)
    {
        return 0.0;
    }
    else
    {
        return 0.5*(1.0+cos(M_PI*x/eps))/eps;
    }
}// cos_delta
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> force_strategy,
    SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> source_strategy,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_delta_fcn("IB_4"),
      d_ghosts(-1),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_visit_writer(NULL),
      d_silo_writer(NULL),
      d_load_balancer(NULL),
      d_ins_hier_integrator(ins_hier_integrator),
      d_lag_data_manager(NULL),
      d_instrument_panel(NULL),
      d_total_flow_volume(),
      d_U_init(NULL),
      d_P_init(NULL),
      d_U_bc_coefs(),
      d_lag_init(NULL),
      d_body_force_set(NULL),
      d_eulerian_force_set(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_eulerian_source_set(NULL),
      d_source_strategy(source_strategy),
      d_source_strategy_needs_init(true),
      d_X_src(),
      d_r_src(),
      d_P_src(),
      d_Q_src(),
      d_n_src(),
      d_using_pIB_method(false),
      d_gravitational_acceleration(NDIM,0.0),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_num_cycles(1),
      d_num_init_cycles(5),
      d_regrid_interval(1),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_dt_max(std::numeric_limits<double>::max()),
      d_dt_max_time_max(std::numeric_limits<double>::max()),
      d_dt_max_time_min(-(d_dt_max_time_max-std::numeric_limits<double>::epsilon())),
      d_is_initialized(false),
      d_do_log(false),
      d_mark_input_file_name(""),
      d_mark_init_posns(),
      d_reinterpolate_after_regrid(false),
      d_hier_cc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_force_current_ralg(),
      d_force_new_ralg(),
      d_source_ralg(),
      d_force_current_rstrategy(),
      d_force_new_rstrategy(),
      d_source_rstrategy(),
      d_force_current_rscheds(),
      d_force_new_rscheds(),
      d_source_rscheds(),
      d_V_var(NULL),
      d_W_var(NULL),
      d_F_var(NULL),
      d_Q_var(NULL),
      d_mark_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_W_idx(-1),
      d_F_idx(-1),
      d_F_scratch1_idx(-1),
      d_F_scratch2_idx(-1),
      d_mark_current_idx(-1),
      d_mark_scratch_idx(-1),
      d_Q_idx(-1),
      d_Q_scratch_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!ins_hier_integrator.isNull());
    TBOX_ASSERT(!force_strategy.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);

    // Read in the initial marker positions.
    if (!from_restart && !d_mark_input_file_name.empty())
    {
        const int mpi_rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
        const int mpi_size = SAMRAI::tbox::SAMRAI_MPI::getNodes();

        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        const double* const grid_xLower = grid_geom->getXLower();
        const double* const grid_xUpper = grid_geom->getXUpper();

        for (int rank = 0; rank < mpi_size; ++rank)
        {
            if (rank == mpi_rank)
            {
                std::string line_string;
                std::ifstream file_stream(d_mark_input_file_name.c_str(), std::ios::in);

                // The first entry in the file is the number of markers.
                int num_marks;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << d_mark_input_file_name << "\n");
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_marks))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << d_mark_input_file_name << "\n");
                    }
                }

                if (num_marks <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << d_mark_input_file_name << "\n");
                }

                // Each successive line provides the initial position of each
                // marker in the input file.
                d_mark_init_posns.resize(NDIM*num_marks,0.0);
                for (int k = 0; k < num_marks; ++k)
                {
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << d_mark_input_file_name << "\n");
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            if (!(line_stream >> d_mark_init_posns[NDIM*k+d]))
                            {
                                TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_mark_input_file_name << "\n");
                            }
                        }

                        // Ensure the initial marker position lies within the
                        // physical domain.
                        const double* const X = &d_mark_init_posns[NDIM*k];
                        for (int d = 0; d < NDIM; ++d)
                        {
                            if (SAMRAI::tbox::MathUtilities<double>::equalEps(X[d],grid_xLower[d]))
                            {
                                TBOX_ERROR(d_object_name << "::IBHierarchyIntegrator():\n"
                                           << "  encountered marker intersecting lower physical boundary.\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                            else if (X[d] <= grid_xLower[d])
                            {
                                TBOX_ERROR(d_object_name << "::IBHierarchyIntegrator():\n"
                                           << "  encountered marker below lower physical boundary\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }

                            if (SAMRAI::tbox::MathUtilities<double>::equalEps(X[d],grid_xUpper[d]))
                            {
                                TBOX_ERROR(d_object_name << "::IBHierarchyIntegrator():\n"
                                           << "  encountered marker intersecting upper physical boundary.\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                            else if (X[d] >= grid_xUpper[d])
                            {
                                TBOX_ERROR(d_object_name << "::IBHierarchyIntegrator():\n"
                                           << "  encountered marker above upper physical boundary\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                        }
                    }
                }
            }
        }
    }

    // Read in the constraint force info.
    d_using_constraint_forces.resize(128,false);  // XXXX
    d_constraint_lag_coarse_idxs.resize(128);
    d_constraint_lag_fine_idxs.resize(128);
    d_constraint_petsc_coarse_idxs.resize(128);
    d_constraint_petsc_fine_idxs.resize(128);
    if (!d_constraint_forces_file_name.empty())
    {
        SAMRAI::tbox::pout << "add restart capability so this doesn't have to be re-read each time!\n";

        const int mpi_rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
        const int mpi_size = SAMRAI::tbox::SAMRAI_MPI::getNodes();

        for (int rank = 0; rank < mpi_size; ++rank)
        {
            if (rank == mpi_rank)
            {
                std::string line_string;
                std::ifstream file_stream(d_constraint_forces_file_name.c_str(), std::ios::in);

                // The first entry in the file is the number of force specs.
                int num_specs;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << d_constraint_forces_file_name << "\n");
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_specs))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << d_constraint_forces_file_name << "\n");
                    }
                }

                // Each successive line provides a force specification for a
                // constraint force.
                for (int k = 0; k < num_specs; ++k)
                {
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        int coarse_ln;
                        if (!(line_stream >> coarse_ln))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                        }
                        int fine_ln;
                        if (!(line_stream >> fine_ln))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                        }
                        int coarse_idx;
                        if (!(line_stream >> coarse_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                        }
                        int fine_idx;
                        if (!(line_stream >> fine_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                        }

                        if (fine_ln != coarse_ln+1)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << d_constraint_forces_file_name << "\n");
                        }

                        d_using_constraint_forces   [coarse_ln] = true;
                        d_constraint_lag_coarse_idxs[coarse_ln].push_back(coarse_idx);
                        d_constraint_lag_fine_idxs  [coarse_ln].push_back(  fine_idx);
                    }
                }
            }
        }
    }

    // Determine the ghost cell width required for cell-centered spreading and
    // interpolating.
    const int stencil_size = IBTK::LEInteractor::getStencilSize(d_delta_fcn);
    d_ghosts = int(floor(0.5*double(stencil_size)))+1;

    // Get the Lagrangian Data Manager.
    d_lag_data_manager = IBTK::LDataManager::getManager(
        d_object_name+"::LDataManager", d_ghosts, d_registered_for_restart);

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(
        d_object_name+"::IBInstrumentPanel",
        (input_db->isDatabase("IBInstrumentPanel")
         ? input_db->getDatabase("IBInstrumentPanel")
         : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL)));

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager = SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->getTimer(
            "IBAMR::IBHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBHierarchyIntegrator

IBHierarchyIntegrator::~IBHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin();
         it != d_rstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin();
         it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        int ierr;
        if (d_constraint_force_src_is[ln] != static_cast<IS>(NULL))
        {
            ierr = ISDestroy(d_constraint_force_src_is[ln]); IBTK_CHKERRQ(ierr);
        }
        if (d_constraint_force_dst_is[ln] != static_cast<IS>(NULL))
        {
            ierr = ISDestroy(d_constraint_force_dst_is[ln]); IBTK_CHKERRQ(ierr);
        }
        if (d_constraint_force_vec_scatter[ln] != static_cast<VecScatter>(NULL))
        {
            ierr = VecScatterDestroy(d_constraint_force_vec_scatter[ln]); IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~IBHierarchyIntegrator

const std::string&
IBHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_U_init = U_init;
    d_ins_hier_integrator->registerVelocityInitialConditions(d_U_init);
    return;
}// registerVelocityInitialConditions

void
IBHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::registerVelocityPhysicalBcCoefs()\n"
                   << "  velocity boundary conditions must be registered prior to initialization\n"
                   << "  of the hierarchy integrator object.\n");
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < U_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(U_bc_coefs[l] != NULL);
    }
#endif
    d_U_bc_coefs = U_bc_coefs;
    d_ins_hier_integrator->registerVelocityPhysicalBcCoefs(d_U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBHierarchyIntegrator::registerPressureInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> P_init)
{
    d_P_init = P_init;
    d_ins_hier_integrator->registerPressureInitialConditions(d_P_init);
    return;
}// registerPressureInitialConditions

void
IBHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> body_force_set)
{
    d_body_force_set = body_force_set;
    return;
}// registerBodyForceSpecification

void
IBHierarchyIntegrator::registerLNodeInitStrategy(
    SAMRAI::tbox::Pointer<IBTK::LNodeInitStrategy> lag_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!lag_init.isNull());
#endif
    d_lag_init = lag_init;
    d_lag_data_manager->registerLNodeInitStrategy(d_lag_init);
    return;
}// registerLNodeInitStrategy

void
IBHierarchyIntegrator::freeLNodeInitStrategy()
{
    d_lag_init.setNull();
    d_lag_data_manager->freeLNodeInitStrategy();
    return;
}// freeLNodeInitStrategy

void
IBHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    d_ins_hier_integrator->registerVisItDataWriter(d_visit_writer);
    d_lag_data_manager->registerVisItDataWriter(d_visit_writer);
    return;
}// registerVisItDataWriter

void
IBHierarchyIntegrator::registerLagSiloDataWriter(
    SAMRAI::tbox::Pointer<IBTK::LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_lag_data_manager->registerLagSiloDataWriter(d_silo_writer);
    return;
}// registerLagSiloDataWriter

#if (NDIM == 3)
void
IBHierarchyIntegrator::registerLagM3DDataWriter(
    SAMRAI::tbox::Pointer<IBTK::LagM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    d_lag_data_manager->registerLagM3DDataWriter(d_m3D_writer);
    return;
}// registerLagM3DDataWriter
#endif

void
IBHierarchyIntegrator::registerLoadBalancer(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_lag_data_manager->registerLoadBalancer(d_load_balancer);
    return;
}// registerLoadBalancer

///
///  The following routines:
///
///      initializeHierarchyIntegrator(),
///      initializeHierarchy(),
///      advanceHierarchy(),
///      atRegridPoint(),
///      getIntegratorTime(),
///      getStartTime(),
///      getEndTime(),
///      getIntegratorStep(),
///      getMaxIntegratorSteps(),
///      stepsRemaining(),
///      getPatchHierarchy(),
///      getGriddingAlgorithm(),
///      getLDataManager(),
///      getIBInstrumentPanel()
///
///  allow the IBHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");
    const SAMRAI::hier::IntVector<NDIM> ghosts = d_ghosts;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    d_V_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::V",NDIM);
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_current, ghosts);

    d_W_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::W",NDIM);
    d_W_idx = var_db->registerVariableAndContext(d_W_var, d_current, ghosts);

    d_F_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::F",NDIM);
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_current, no_ghosts);
    d_F_scratch1_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);
    d_F_scratch2_idx = var_db->registerClonedPatchDataIndex(d_F_var, d_F_scratch1_idx);

    d_mark_var = new SAMRAI::pdat::IndexVariable<NDIM,IBTK::LagMarker>(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, d_scratch, ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

    if (!d_source_strategy.isNull())
    {
        const SAMRAI::hier::IntVector<NDIM> source_ghosts = 1;
        d_Q_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Q",1);
        d_Q_idx = var_db->registerVariableAndContext(d_Q_var, d_current, no_ghosts);
        d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_var, d_scratch, source_ghosts);
    }

    // Initialize the objects used to manage Lagragian-Eulerian interaction.
    //
    // NOTE: The IBEulerianForceSetter only has to set the new Cartesian grid
    // force.  The current Cartesian grid force is set manually by
    // IBHierarchyIntegrator::advanceHierarchy().
    d_eulerian_force_set = new IBEulerianForceSetter(
        d_object_name+"::IBEulerianForceSetter", -1, d_F_idx, -1);
    d_ins_hier_integrator->registerBodyForceSpecification(d_eulerian_force_set);

    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_set = new IBEulerianSourceSetter(
            d_object_name+"::IBEulerianSourceSetter", d_Q_idx, d_Q_idx, d_Q_idx);
        d_ins_hier_integrator->registerDivergenceSpecification(
            d_eulerian_source_set);
    }

    // Initialize the INSHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(
        d_V_idx,       // destination
        U_current_idx, // source
        d_V_idx,       // temporary work space
        refine_operator);
    d_rstrategies["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new IBTK::CartExtrapPhysBdryOp(d_V_idx, "QUADRATIC");

    const int U_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getNewContext());

    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(
        d_W_idx,   // destination
        U_new_idx, // source
        d_W_idx,   // temporary work space
        refine_operator);
    d_rstrategies["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"] = new IBTK::CartExtrapPhysBdryOp(d_W_idx, "QUADRATIC");

    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getScratchContext());
    const int P_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getNewContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getScratchContext());

    d_ralgs["INSTRUMENTATION_DATA_FILL"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getVelocityVar(),
        "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(
        U_scratch_idx,  // destination
        U_new_idx,      // source
        U_scratch_idx,  // temporary work space
        refine_operator);

    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getPressureVar(),
        "LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(
        P_scratch_idx,  // destination
        P_new_idx,      // source
        P_scratch_idx,  // temporary work space
        refine_operator);

    SAMRAI::hier::ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(U_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(P_scratch_idx);
    d_rstrategies["INSTRUMENTATION_DATA_FILL"] = new IBTK::CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "QUADRATIC");

    // NOTE: When using conservative averaging to coarsen the velocity from
    // finer levels to coarser levels, the appropriate prolongation operator for
    // the force is constant refinement.  This choice results in IB spreading
    // and interpolation being adjoint.
    const int F_current_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getForceVar(),
        d_ins_hier_integrator->getCurrentContext());

    d_force_current_ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getForceVar(), "CONSTANT_REFINE");
    d_force_current_ralg->registerRefine(F_current_idx,     // destination
                                         F_current_idx,     // source
                                         d_F_scratch1_idx,  // temporary work space
                                         refine_operator);
    d_force_current_rstrategy = new IBTK::CartExtrapPhysBdryOp(
        d_F_scratch1_idx, "QUADRATIC");

    d_force_new_ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(
        d_ins_hier_integrator->getForceVar(), "CONSTANT_REFINE");
    d_force_new_ralg->registerRefine(d_F_idx,           // destination
                                     d_F_idx,           // source
                                     d_F_scratch2_idx,  // temporary work space
                                     refine_operator);
    d_force_new_rstrategy = new IBTK::CartExtrapPhysBdryOp(
        d_F_scratch2_idx, "QUADRATIC");

    if (!d_source_strategy.isNull())
    {
        d_source_ralg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        refine_operator = grid_geom->lookupRefineOperator(
            d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
        d_source_ralg->registerRefine(d_Q_idx,          // destination
                                      d_Q_idx,          // source
                                      d_Q_scratch_idx,  // temporary work space
                                      refine_operator);
        d_source_rstrategy = new IBTK::CartExtrapPhysBdryOp(
            d_Q_scratch_idx, "QUADRATIC");
    }

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(
        U_current_idx, // destination
        U_current_idx, // source
        coarsen_operator);

    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_W_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::N->N::CONSERVATIVE_COARSEN"]->registerCoarsen(
        U_new_idx, // destination
        U_new_idx, // source
        coarsen_operator);

    d_calgs["F->F::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_ins_hier_integrator->getForceVar(), "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->registerCoarsen(
        F_current_idx, // destination
        F_current_idx, // source
        coarsen_operator);
    coarsen_operator = grid_geom->lookupCoarsenOperator(
        d_F_var, "CONSERVATIVE_COARSEN");
    d_calgs["F->F::CONSERVATIVE_COARSEN"]->registerCoarsen(
        d_F_idx, // destination
        d_F_idx, // source
        coarsen_operator);

    if (!d_source_strategy.isNull())
    {
        d_calgs["Q->Q::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
        coarsen_operator = grid_geom->lookupCoarsenOperator(
            d_Q_var, "CONSERVATIVE_COARSEN");
        d_calgs["Q->Q::CONSERVATIVE_COARSEN"]->registerCoarsen(
            d_Q_idx, // destination
            d_Q_idx, // source
            coarsen_operator);
    }

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
IBHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    // Use the INSHierarchyIntegrator to initialize the patch hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(dt_next, d_dt_max);
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Begin Lagrangian data movement.
    d_lag_data_manager->beginDataRedistribution();

    // Finish Lagrangian data movement.
    d_lag_data_manager->endDataRedistribution();

    // Update the workload.
    d_lag_data_manager->updateWorkloadData(
        coarsest_ln, finest_ln);

    // Initialize the instrumentation data.
    d_instrument_panel->initializeHierarchyIndependentData(
        d_hierarchy, d_lag_data_manager);
    if (d_instrument_panel->isInstrumented())
    {
        d_instrument_panel->initializeHierarchyDependentData(
            d_hierarchy, d_lag_data_manager, d_integrator_step, d_integrator_time);
        if (d_total_flow_volume.empty())
        {
            d_total_flow_volume.resize(d_instrument_panel->getFlowValues().size(),0.0);
        }
    }

    // Indicate that the force and source strategies need to be re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
#endif

    const double current_time = d_integrator_time;
    const double new_time     = d_integrator_time+dt;
    const bool initial_time   = SAMRAI::tbox::MathUtilities<double>::equalEps(current_time,d_start_time);

    // Regrid the patch hierarchy.
    const bool do_regrid = (d_regrid_interval == 0
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Set the current time interval in the force specification objects.
    d_eulerian_force_set->setTimeInterval(current_time, new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);
    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_set->setTimeInterval(current_time, new_time);
        d_source_strategy->setTimeInterval(current_time, new_time);
    }

    // (Re)initialize the force and source strategies.
    if (d_force_strategy_needs_init)
    {
        resetLagrangianForceStrategy(current_time, initial_time);
        d_force_strategy_needs_init = false;
    }
    if (d_source_strategy_needs_init && !d_source_strategy.isNull())
    {
        resetLagrangianSourceStrategy(current_time, initial_time);
        d_source_strategy_needs_init = false;
    }

    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_new_data(finest_ln+1);

    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > K_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > M_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > Y_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > dY_dt_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_K_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > Y_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > dY_dt_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_K_new_data(finest_ln+1);

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy.
    //
    // NOTE: Since we are maintaining the Lagrangian velocity data, this step is
    // skipped for each timestep following the initial one, except for those
    // timesteps that immediately follow a regridding.
    if (initial_time || d_reinterpolate_after_regrid)
    {
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            d_cscheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_V_idx, current_time);
            d_rscheds["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
        }
    }

    // Compute the initial updated marker positions.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > v_data = patch->getPatchData(d_V_idx);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);

            // Collect the positions and velocities of all markers in the patch.
            std::vector<double> U_mark, X_mark;
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const IBTK::LagMarker& mark = it();
                    const std::vector<double>& X = mark.getPositions();
                    const std::vector<double>& U = mark.getVelocities();
                    X_mark.insert(X_mark.end(), X.begin(), X.end());
                    U_mark.insert(U_mark.end(), U.begin(), U.end());
                }
            }

            // When necessary, interpolate the velocity field.
            if (initial_time || d_reinterpolate_after_regrid)
            {
                IBTK::LEInteractor::interpolate(
                    U_mark, NDIM, X_mark, NDIM, v_data,
                    patch, patch_box, d_delta_fcn);
            }

            // Update the marker positions.
            for (size_t k = 0; k < X_mark.size()/NDIM; ++k)
            {
                double* const X = &X_mark[NDIM*k];
                const double* const U = &U_mark[NDIM*k];
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] += dt*U[d];
                }
            }

            // Store the marker positions and velocities in the scratch marker
            // patch data.
            mark_scratch_data->copy(*mark_data);
            int marker_offset = 0;
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_scratch_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    IBTK::LagMarker& mark = it();
                    const int nmarks = mark.getNumberOfMarkers();
                    std::vector<double> X(X_mark.begin()+NDIM*marker_offset,
                                          X_mark.begin()+NDIM*(marker_offset+nmarks));
                    std::vector<double> U(U_mark.begin()+NDIM*marker_offset,
                                          U_mark.begin()+NDIM*(marker_offset+nmarks));
                    mark.setPositions(X);
                    mark.setVelocities(U);
                    marker_offset += nmarks;
                }
            }
        }
    }

    // Initialize the various LNodeLevelData objects, including X_data, F_data,
    // U_data, and X_new_data, on each level of the patch hierarchy.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::VEL_DATA_NAME,ln);
            F_data[ln] = d_lag_data_manager->createLNodeLevelData("F",ln,NDIM);

            X_data[ln]->restoreLocalFormVec();
            U_data[ln]->restoreLocalFormVec();
            F_data[ln]->restoreLocalFormVec();

            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);
            U_new_data[ln] = d_lag_data_manager->createLNodeLevelData("U_new",ln,NDIM);
            F_new_data[ln] = d_lag_data_manager->createLNodeLevelData("F_new",ln,NDIM);

            X_new_data[ln]->restoreLocalFormVec();
            U_new_data[ln]->restoreLocalFormVec();
            F_new_data[ln]->restoreLocalFormVec();

            if (d_using_pIB_method)
            {
                K_data[ln]         = d_lag_data_manager->getLNodeLevelData("K",ln);
                M_data[ln]         = d_lag_data_manager->getLNodeLevelData("M",ln);
                Y_data    [ln]     = d_lag_data_manager->getLNodeLevelData("Y",ln);
                dY_dt_data[ln]     = d_lag_data_manager->getLNodeLevelData("dY_dt",ln);
                F_K_data  [ln]     = d_lag_data_manager->createLNodeLevelData("F_K",ln,NDIM);
                Y_new_data    [ln] = d_lag_data_manager->createLNodeLevelData("Y_new",ln,NDIM);
                dY_dt_new_data[ln] = d_lag_data_manager->createLNodeLevelData("dY_dt_new",ln,NDIM);
                F_K_new_data  [ln] = d_lag_data_manager->createLNodeLevelData("F_K_new",ln,NDIM);

                K_data[ln]        ->restoreLocalFormVec();
                M_data[ln]        ->restoreLocalFormVec();
                Y_data    [ln]    ->restoreLocalFormVec();
                dY_dt_data[ln]    ->restoreLocalFormVec();
                F_K_data  [ln]    ->restoreLocalFormVec();
                Y_new_data    [ln]->restoreLocalFormVec();
                dY_dt_new_data[ln]->restoreLocalFormVec();
                F_K_new_data  [ln]->restoreLocalFormVec();
            }
        }
    }

    // 1. Interpolate u(n) from the Cartesian grid onto the Lagrangian mesh.
    //
    // NOTE: Since we are maintaining the Lagrangian velocity data, this step is
    // skipped for each timestep following the initial one execpt immediately
    // following a regridding operation.
    if (initial_time || d_reinterpolate_after_regrid)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n) to U(n) on level number " << ln << "\n";
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > v_data = patch->getPatchData(d_V_idx);
                    const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                        d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                    IBTK::LEInteractor::interpolate(
                        U_data[ln], X_data[ln], idx_data, v_data,
                        patch, patch_box, periodic_shift,
                        d_delta_fcn);
                }
            }
        }
    }

    // 2. Compute F(n) = F(X(n),n), the Lagrangian force corresponding to
    //    configuration X(n) at time t_{n}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F(n) on level number " << ln << "\n";
            Vec F_vec = F_data[ln]->getGlobalVec();
            int ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_data[ln], X_data[ln],
                d_hierarchy, ln, current_time, d_lag_data_manager);
            if (d_using_pIB_method)
            {
                // Add the penalty force associated with the massive ghost
                // particles, i.e., F_K = K*(Y-X).
                Vec K_vec = K_data[ln]->getGlobalVec();
                Vec F_vec = F_data[ln]->getGlobalVec();
                Vec X_vec = X_data[ln]->getGlobalVec();
                Vec Y_vec = Y_data[ln]->getGlobalVec();
                Vec F_K_vec = F_K_data[ln]->getGlobalVec();

                if (initial_time)
                {
                    // The initial positions of the massive ghost particles
                    // should be the same as the initial positions of the nodes
                    // of the Lagrangian mesh.
                    ierr = VecCopy(X_vec, Y_vec);  IBTK_CHKERRQ(ierr);
                }

                int n_local = 0;
                ierr = VecGetLocalSize(K_vec, &n_local);  IBTK_CHKERRQ(ierr);

                double* K_arr, * X_arr, * Y_arr, * F_K_arr;
                ierr = VecGetArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(F_K_vec, &F_K_arr);  IBTK_CHKERRQ(ierr);

                static double max_displacement = 0.0;
                double max_config_displacement = 0.0;
                for (int i = 0; i < n_local; ++i)
                {
                    double displacement = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_K_arr[NDIM*i+d] = K_arr[i]*(Y_arr[NDIM*i+d] - X_arr[NDIM*i+d]);
                        displacement += pow(Y_arr[NDIM*i+d]-X_arr[NDIM*i+d],2.0);
                    }
                    displacement = sqrt(displacement);
                    if (displacement > max_config_displacement)
                    {
                        max_config_displacement = displacement;
                    }
                }
                max_config_displacement = SAMRAI::tbox::SAMRAI_MPI::maxReduction(max_config_displacement);
                if (max_config_displacement > max_displacement)
                {
                    max_displacement = max_config_displacement;
                }
                if (d_do_log && !SAMRAI::tbox::MathUtilities<double>::equalEps(max_config_displacement,0.0))
                {
                    SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy():\n";
                    SAMRAI::tbox::plog << "  maximum massive boundary point displacement [present configuration] = " << max_config_displacement << "\n";
                    SAMRAI::tbox::plog << "  maximum massive boundary point displacement [entire simulation] = " << max_displacement << "\n";
                }

                ierr = VecRestoreArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(F_K_vec, &F_K_arr);  IBTK_CHKERRQ(ierr);

                ierr = VecAXPY(F_vec, 1.0, F_K_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    computeConstraintForces(X_data, F_data, coarsest_ln, finest_ln, current_time, initial_time);

    // 3. Spread F(n) from the Lagrangian mesh onto the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        level->allocatePatchData(d_F_scratch1_idx, current_time);

        // On the coarsest level of the patch hierarchy, simply initialize the
        // Cartesian force density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy, conservatively
        // interpolate the Cartesian force density from coarser levels in the
        // patch hierarchy.
        if (ln == coarsest_ln)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch->getPatchData(
                    d_ins_hier_integrator->getForceVar(), d_ins_hier_integrator->getCurrentContext());
                f_current_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser level in
            // the hierarchy.
            //
            // Note that these refine schedules initialize both the force data
            // maintained by the IBHierarchyIntegrator and the current context
            // of the force data maintained by the Navier-Stokes solver.
            d_force_current_rscheds[ln]->fillData(current_time);
        }

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln]->beginGhostUpdate();
            F_data[ln]->beginGhostUpdate();
            X_data[ln]->endGhostUpdate();
            F_data[ln]->endGhostUpdate();
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F(n) to f(n) on level number " << ln << "\n";
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch->getPatchData(
                    d_ins_hier_integrator->getForceVar(), d_ins_hier_integrator->getCurrentContext());
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                IBTK::LEInteractor::spread(
                    f_current_data, F_data[ln], X_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts), periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // 4. Compute X~(n+1), the preliminary structure configuration at time
    //    t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X~(n+1) on level number " << ln << "\n";
            Vec U_vec = U_data[ln]->getGlobalVec();
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            int ierr = VecWAXPY(X_new_vec, dt, U_vec, X_vec);  IBTK_CHKERRQ(ierr);
            if (d_using_pIB_method)
            {
                // Advance the positions and velocities of the massive ghost
                // particles forward in time via forward Euler.
                Vec M_vec = M_data[ln]->getGlobalVec();

                Vec Y_vec = Y_data[ln]->getGlobalVec();
                Vec Y_new_vec = Y_new_data[ln]->getGlobalVec();

                Vec dY_dt_vec = dY_dt_data[ln]->getGlobalVec();
                Vec dY_dt_new_vec = dY_dt_new_data[ln]->getGlobalVec();

                Vec F_K_vec = F_K_data[ln]->getGlobalVec();

                if (initial_time)
                {
                    // The initial velocities of the massive ghost particles
                    // should be the same as the interpolated velocity field
                    // evaluated at the nodes of the Lagrangian mesh.
                    ierr = VecCopy(U_vec, dY_dt_vec);  IBTK_CHKERRQ(ierr);
                }

                int n_local = 0;
                ierr = VecGetLocalSize(M_vec, &n_local);  IBTK_CHKERRQ(ierr);

                double* M_arr, * Y_arr, * Y_new_arr, * dY_dt_arr, * dY_dt_new_arr, * F_K_arr;
                ierr = VecGetArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(dY_dt_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(F_K_vec, &F_K_arr);  IBTK_CHKERRQ(ierr);

                for (int i = 0; i < n_local; ++i)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        Y_new_arr[NDIM*i+d] = Y_arr[NDIM*i+d] + dt*dY_dt_arr[NDIM*i+d];
                        dY_dt_new_arr[NDIM*i+d] = dY_dt_arr[NDIM*i+d] - (dt/M_arr[i])*F_K_arr[NDIM*i+d] + dt*d_gravitational_acceleration[d];
                    }
                }

                ierr = VecRestoreArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(Y_vec, &Y_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(dY_dt_vec, &dY_dt_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(F_K_vec, &F_K_arr);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    computeConstraintForces(X_new_data, F_new_data, coarsest_ln, finest_ln, new_time, false);

    // 5. Compute F~(n+1) = F(X~(n+1),n+1), the Lagrangian force corresponding
    //    to configuration X~(n+1) at time t_{n+1}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing F~(n+1) on level number " << ln << "\n";
            Vec F_new_vec = F_new_data[ln]->getGlobalVec();
            int ierr = VecSet(F_new_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(
                F_new_data[ln], X_new_data[ln],
                d_hierarchy, ln, new_time, d_lag_data_manager);
            if (d_using_pIB_method)
            {
                // Add the penalty force associated with the massive ghost
                // particles, i.e., F_K = K*(Y-X).
                Vec K_vec = K_data[ln]->getGlobalVec();
                Vec F_new_vec = F_new_data[ln]->getGlobalVec();
                Vec X_new_vec = X_new_data[ln]->getGlobalVec();
                Vec Y_new_vec = Y_new_data[ln]->getGlobalVec();
                Vec F_K_new_vec = F_K_new_data[ln]->getGlobalVec();

                int n_local = 0;
                ierr = VecGetLocalSize(K_vec, &n_local);  IBTK_CHKERRQ(ierr);

                double* K_arr, * X_new_arr, * Y_new_arr, * F_K_new_arr;
                ierr = VecGetArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(X_new_vec, &X_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(F_K_new_vec, &F_K_new_arr);  IBTK_CHKERRQ(ierr);

                for (int i = 0; i < n_local; ++i)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_K_new_arr[NDIM*i+d] = K_arr[i]*(Y_new_arr[NDIM*i+d] - X_new_arr[NDIM*i+d]);
                    }
                }

                ierr = VecRestoreArray(K_vec, &K_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(X_new_vec, &X_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(F_K_new_vec, &F_K_new_arr);  IBTK_CHKERRQ(ierr);

                ierr = VecAXPY(F_new_vec, 1.0, F_K_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // 6. Spread F~(n+1) from the Lagrangian mesh onto the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        level->allocatePatchData(d_F_idx, current_time);
        level->allocatePatchData(d_F_scratch2_idx, current_time);

        // On the coarsest level of the patch hierarchy, simply initialize the
        // Cartesian force density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy, conservatively
        // interpolate the Cartesian force density from coarser levels in the
        // patch hierarchy.
        if (ln == coarsest_ln)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->getPatchData(d_F_idx);
                f_new_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian force from the next coarser level in
            // the hierarchy.
            //
            // Note that these refine schedules initialize both the force data
            // maintained by the IBHierarchyIntegrator and the current context
            // of the force data maintained by the Navier-Stokes solver.
            d_force_new_rscheds[ln]->fillData(current_time);
        }

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_new_data[ln]->beginGhostUpdate();
            F_new_data[ln]->beginGhostUpdate();
            X_new_data[ln]->endGhostUpdate();
            F_new_data[ln]->endGhostUpdate();
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading F~(n+1) to f~(n+1) on level number " << ln << "\n";
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch->getPatchData(d_F_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                IBTK::LEInteractor::spread(
                    f_new_data, F_new_data[ln], X_new_data[ln], idx_data,
                    patch, SAMRAI::hier::Box<NDIM>::grow(patch_box,d_ghosts), periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // 7. If an additional body force specification object is provided, compute
    //    the body force at the beginning and end of the time step, and add
    //    those values to f(n) and f(n+1).
    if (!d_body_force_set.isNull())
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int F_current_idx = var_db->mapVariableAndContextToIndex(
            d_ins_hier_integrator->getForceVar(),
            d_ins_hier_integrator->getCurrentContext());

        d_body_force_set->setDataOnPatchHierarchy(
            d_F_scratch1_idx, d_F_var, d_hierarchy,
            current_time, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(F_current_idx, F_current_idx, d_F_scratch1_idx);

        d_body_force_set->setDataOnPatchHierarchy(
            d_F_scratch2_idx, d_F_var, d_hierarchy,
            current_time+dt, false, coarsest_ln, finest_ln);

        d_hier_cc_data_ops->add(d_F_idx, d_F_idx, d_F_scratch2_idx);
    }

    // Synchronize the Cartesian grid force field on the patch hierarchy.
    //
    // Note that these coarsen schedules synchronize both the force data
    // maintained by the IBHierarchyIntegrator and the current context of the
    // force data maintained by the Navier-Stokes solver.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["F->F::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }

    // Deallocate unneeded scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (initial_time || d_reinterpolate_after_regrid) level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_scratch1_idx);
        level->deallocatePatchData(d_F_scratch2_idx);
    }

    // Compute the source/sink strengths corresponding to any distributed
    // internal fluid sources or sinks.
    if (!d_source_strategy.isNull())
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_Q_idx, current_time);
            level->allocatePatchData(d_Q_scratch_idx, current_time);
        }
        computeSourceStrengths(coarsest_ln, finest_ln, current_time, X_data);
    }

    // Reset the reinterpolation flag.
    d_reinterpolate_after_regrid = false;

    // Solve the incompressible Navier-Stokes equations for U(n+1) and P(n+1/2).
    //
    // Note that the pressure at start_time is not an initial value for the
    // incompressible Navier-Stokes equations, so we must solve for it by
    // cycling the solution for the first timestep:
    //
    // Solve the Navier-Stokes equations with initial guess P(n=0)=0, solve for
    // P(n=1/2), discard U(n=1) and u(n=1), rinse, wash, repeat.
    //
    // For all other timesteps, we just use the previous value of P as the guess
    // for P(n+1/2).
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing u(n+1), p(n+1/2)\n";

    int cycle = 0;
    const bool performing_init_cycles = initial_time;
    const int num_cycles = (performing_init_cycles
                            ? d_num_init_cycles
                            : d_num_cycles);
    for (cycle = 0; cycle < num_cycles; ++cycle)
    {
        if (d_do_log && performing_init_cycles)
        {
            SAMRAI::tbox::plog << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "+ Performing cycle " << cycle+1 << " of " << d_num_init_cycles << " to initialize P(n=1/2)\n";
            SAMRAI::tbox::plog << "+\n";
            SAMRAI::tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
        }

        // Solve the Navier-Stokes equations for U(n+1), u(n+1), P(n+1/2).  For
        // better or worse, each of the major algorithmic steps of the
        // projection method is separated into its own member function.
        d_ins_hier_integrator->predictAdvectionVelocity(current_time, new_time);

        d_ins_hier_integrator->integrateAdvDiff(current_time, new_time);

        d_ins_hier_integrator->projectVelocity(current_time, new_time);

        d_ins_hier_integrator->updatePressure(current_time, new_time, true);

        d_ins_hier_integrator->synchronizeHierarchy();

        if (cycle < (num_cycles-1))
        {
            // Reset data to the preadvance state.
            d_ins_hier_integrator->resetHierDataToPreadvanceState();

            // Reset time of all current data.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(current_time);
            }
        }
        else
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                level->setTime(new_time);
            }
        }
    }

    // Synchronize the Cartesian grid velocity u(n+1) on the patch hierarchy.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::N->N::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_W_idx, current_time);
        d_rscheds["U->W::N->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }

    // In the following loop, we perform the following operations:
    //
    // 1. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh
    //    using X~(n+1).
    //
    // 2. Compute X(n+1), the final structure configuration at time t_{n+1}.
    //
    // 3. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian mesh
    //    using X(n+1).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(
            level->getRatio());

        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;

            // 1. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian
            //    mesh.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data = patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                IBTK::LEInteractor::interpolate(
                    U_new_data[ln], X_new_data[ln], idx_data, w_data,
                    patch, patch_box, periodic_shift,
                    d_delta_fcn);
            }

            // 2. Compute X(n+1), the final structure configuration at time
            //    t_{n+1}.
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing X(n+1) on level number " << ln << "\n";

            Vec U_new_vec = U_new_data[ln]->getGlobalVec();
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            ierr = VecAXPY(X_new_vec, dt, U_new_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPBY(X_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);

            if (d_using_pIB_method)
            {
                // Advance the positions and velocities of the massive ghost
                // particles forward in time via 2nd order SSP Runge-Kutta.
                Vec M_vec = M_data[ln]->getGlobalVec();

                Vec Y_vec = Y_data[ln]->getGlobalVec();
                Vec Y_new_vec = Y_new_data[ln]->getGlobalVec();

                Vec dY_dt_vec = dY_dt_data[ln]->getGlobalVec();
                Vec dY_dt_new_vec = dY_dt_new_data[ln]->getGlobalVec();

                Vec F_K_new_vec = F_K_new_data[ln]->getGlobalVec();

                int n_local = 0;
                ierr = VecGetLocalSize(M_vec, &n_local);  IBTK_CHKERRQ(ierr);

                double* M_arr, * Y_new_arr, * dY_dt_new_arr, * F_K_new_arr;
                ierr = VecGetArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(F_K_new_vec, &F_K_new_arr);  IBTK_CHKERRQ(ierr);

                for (int i = 0; i < n_local; ++i)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        Y_new_arr[NDIM*i+d] = Y_new_arr[NDIM*i+d] + dt*dY_dt_new_arr[NDIM*i+d];
                        dY_dt_new_arr[NDIM*i+d] = dY_dt_new_arr[NDIM*i+d] - (dt/M_arr[i])*F_K_new_arr[NDIM*i+d] + dt*d_gravitational_acceleration[d];
                    }
                }

                ierr = VecRestoreArray(M_vec, &M_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(Y_new_vec, &Y_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(dY_dt_new_vec, &dY_dt_new_arr);  IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArray(F_K_new_vec, &F_K_new_arr);  IBTK_CHKERRQ(ierr);

                ierr = VecAXPBY(Y_vec, 0.5, 0.5, Y_new_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecAXPBY(dY_dt_vec, 0.5, 0.5, dY_dt_new_vec);  IBTK_CHKERRQ(ierr);
            }

            // 3. Interpolate u(n+1) from the Cartesian grid onto the Lagrangian
            //    mesh using X(n+1).
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): interpolating u(n+1) to U(n+1) on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data = patch->getPatchData(d_W_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(
                    d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());

                IBTK::LEInteractor::interpolate(
                    U_data[ln], X_data[ln], idx_data, w_data,
                    patch, patch_box, periodic_shift,
                    d_delta_fcn);
            }
        }
    }

    // Compute the final updated marker positions.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > w_data = patch->getPatchData(d_W_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_data = patch->getPatchData(d_mark_current_idx);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);

            // Collect the current and predicted new positions of all markers in
            // the patch.
            std::vector<double> X_mark_current;
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const IBTK::LagMarker& mark = it();
                    const std::vector<double>& X = mark.getPositions();
                    X_mark_current.insert(X_mark_current.end(), X.begin(), X.end());
                }
            }

            std::vector<double> X_mark_new;
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_scratch_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const IBTK::LagMarker& mark = it();
                    const std::vector<double>& X = mark.getPositions();
                    X_mark_new.insert(X_mark_new.end(), X.begin(), X.end());
                }
            }

            // Interpolate the velocity field at the predicted new marker
            // positions.
            std::vector<double> U_mark_new(X_mark_new.size(),0.0);
            IBTK::LEInteractor::interpolate(
                U_mark_new, NDIM, X_mark_new, NDIM, w_data,
                patch, patch_box, d_delta_fcn);

            // Update the marker positions.
            for (size_t k = 0; k < X_mark_new.size()/NDIM; ++k)
            {
                double* const X_current = &X_mark_current[NDIM*k];
                double* const X_new = &X_mark_new[NDIM*k];
                const double* const U_new = &U_mark_new[NDIM*k];
                for (int d = 0; d < NDIM; ++d)
                {
                    X_new[d] += dt*U_new[d];
                    X_new[d] = 0.5*(X_current[d]+X_new[d]);
                }
            }

            // Re-interpolate the velocity field at the final marker positions.
            IBTK::LEInteractor::interpolate(
                U_mark_new, NDIM, X_mark_new, NDIM, w_data,
                patch, patch_box, d_delta_fcn);

            // Prevent markers from leaving the computational domain.
            static const double eps = sqrt(std::numeric_limits<double>::epsilon());
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
            if (periodic_shift.min() == 0)
            {
                const double* const xLower = grid_geom->getXLower();
                const double* const xUpper = grid_geom->getXUpper();
                for (size_t k = 0; k < X_mark_new.size()/NDIM; ++k)
                {
                    double* const X = &X_mark_new[NDIM*k];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (periodic_shift[d] == 0) X[d] = std::min(std::max(X[d],xLower[d]+eps),xUpper[d]-eps);
                    }
                }
            }

            // Store the marker positions and velocities in the current marker
            // patch data.
            int marker_offset = 0;
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    IBTK::LagMarker& mark = it();
                    const int nmarks = mark.getNumberOfMarkers();
                    std::vector<double> X(X_mark_new.begin()+NDIM*marker_offset,
                                          X_mark_new.begin()+NDIM*(marker_offset+nmarks));
                    std::vector<double> U(U_mark_new.begin()+NDIM*marker_offset,
                                          U_mark_new.begin()+NDIM*(marker_offset+nmarks));
                    mark.setPositions(X);
                    mark.setVelocities(U);
                    marker_offset += nmarks;
                }
            }
        }
        level->deallocatePatchData(d_mark_scratch_idx);
    }

    // Update the instrumentation data.
    updateIBInstrumentationData(d_integrator_step+1,new_time);
    if (d_instrument_panel->isInstrumented())
    {
        const std::vector<std::string>& instrument_name = d_instrument_panel->getInstrumentNames();
        const std::vector<double>& flow_data = d_instrument_panel->getFlowValues();
        for (unsigned m = 0; m < flow_data.size(); ++m)
        {
            d_total_flow_volume[m] += flow_data[m]*dt;
            SAMRAI::tbox::plog << "flow volume through " << instrument_name[m] << ":\t " << d_total_flow_volume[m] << "\n";
        }
        SAMRAI::tbox::plog << "NOTE: flow volume in default units\n";
    }

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (!d_source_strategy.isNull())
    {
        computeSourcePressures(coarsest_ln, finest_ln, new_time, X_data);
    }

    // Deallocate any remaining scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_W_idx);
        level->deallocatePatchData(d_F_idx);
        if (!d_source_strategy.isNull())
        {
            level->deallocatePatchData(d_Q_idx);
            level->deallocatePatchData(d_Q_scratch_idx);
        }
    }

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep.
    double dt_next = d_ins_hier_integrator->getStableTimestep();

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min &&
        d_integrator_time <= d_dt_max_time_max)
    {
        dt_next = std::min(d_dt_max,dt_next);
    }

    if (!initial_time)
    {
        dt_next = std::min(dt_next,d_grow_dt*d_old_dt);
    }

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

bool
IBHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;

    return ((d_integrator_step>0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBHierarchyIntegrator::stepsRemaining() const
{
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

IBTK::LDataManager*
IBHierarchyIntegrator::getLDataManager() const
{
    return d_lag_data_manager;
}// getLDataManager

SAMRAI::tbox::Pointer<IBInstrumentPanel>
IBHierarchyIntegrator::getIBInstrumentPanel() const
{
    return d_instrument_panel;
}// getIBInstrumentPanel

///
///  The following routines:
///
///      regridHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the IBHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);

    // Update the marker data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): resetting markers particles.\n";
    const int num_marks = countMarkers(0,d_hierarchy->getFinestLevelNumber());
    resetMarkersOnPatchHierarchy();
    collectMarkersOnPatchHierarchy();

    // Update the workload pre-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Before regriding, begin Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): calling INSHierarchyIntegrator::regridHierarchy().\n";
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Ensure that we haven't misplaced any of the markers.
    const int num_marks_after_regrid = countMarkers(0,d_hierarchy->getFinestLevelNumber());
    if (num_marks != num_marks_after_regrid)
    {
        TBOX_ERROR(d_object_name << "::regridHierarchy()\n"
                   << "  number of marker particles changed during regrid\n"
                   << "  number of markers before regrid = " << num_marks << "\n"
                   << "  number of markers after  regrid = " << num_marks_after_regrid << "\n");
    }

    // Indicate that the force and source strategies need to be re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;

    // Update any constraint force information.
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(d_hierarchy->getFinestLevelNumber()+1);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
        }
    }
    d_constraint_force_src_is.resize(d_hierarchy->getFinestLevelNumber()+1,static_cast<IS>(NULL));
    d_constraint_force_dst_is.resize(d_hierarchy->getFinestLevelNumber()+1,static_cast<IS>(NULL));
    d_constraint_force_vec_scatter.resize(d_hierarchy->getFinestLevelNumber()+1,static_cast<VecScatter>(NULL));
    computeConstraintForceDataStructures(X_data,0,d_hierarchy->getFinestLevelNumber(), d_integrator_time, initial_time);

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->synchronizeNewLevels(
        hierarchy, coarsest_level, finest_level,
        sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetHierDataToPreadvanceState();

    t_reset_data_to_preadvance_state->stop();
    return;
}// resetHierDataToPreadvanceState

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
IBHierarchyIntegrator::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->initializeLevelData(
        hierarchy, level_number, init_data_time,
        can_be_refined, initial_time, old_level,
        allocate_data);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_lag_data_manager->setPatchHierarchy(hierarchy);
    d_lag_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());

    d_lag_data_manager->initializeLevelData(
        hierarchy, level_number, init_data_time,
        can_be_refined, initial_time, old_level,
        allocate_data);

    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    if (allocate_data)
    {
        level->allocatePatchData(d_mark_current_idx, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_mark_current_idx);
    }

    // Copy data on the coarsest level of the patch hierarchy; otherwise, fill
    // data from coarsesr levels when available.
    if (!old_level.isNull() && level_number == 0)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(old_level->getLevelNumber() == level_number);
#endif
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > copy_mark_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        copy_mark_alg->registerRefine(d_mark_current_idx, // destination
                                      d_mark_current_idx, // source
                                      d_mark_scratch_idx, // temporary work space
                                      NULL);

        level->allocatePatchData(d_mark_scratch_idx, init_data_time);
        copy_mark_alg->createSchedule(level, old_level)->fillData(init_data_time);
        level->deallocatePatchData(d_mark_scratch_idx);
    }
    else if (level_number > 0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_mark_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        refine_mark_alg->registerRefine(d_mark_current_idx, // destination
                                        d_mark_current_idx, // source
                                        d_mark_scratch_idx, // temporary work space
                                        new IBTK::LagMarkerRefineOperator());

        level->allocatePatchData(d_mark_scratch_idx, init_data_time);

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > dst_level = level;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > src_level = old_level;
        SAMRAI::xfer::RefinePatchStrategy<NDIM>* refine_mark_op = NULL;
        refine_mark_alg->createSchedule(dst_level, src_level, level_number-1, hierarchy, refine_mark_op)->fillData(d_integrator_time);

        level->deallocatePatchData(d_mark_scratch_idx);
    }

    // Setup the pIB data at the inital time only.
    if (initial_time && d_lag_init->getLevelHasLagrangianData(level_number, can_be_refined))
    {
        static const bool manage_data = true;
        if (d_using_pIB_method)
        {
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> M_data = d_lag_data_manager->createLNodeLevelData("M",level_number,1,manage_data);
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> K_data = d_lag_data_manager->createLNodeLevelData("K",level_number,1,manage_data);
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> Y_data = d_lag_data_manager->createLNodeLevelData("Y",level_number,NDIM,manage_data);
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> dY_dt_data = d_lag_data_manager->createLNodeLevelData("dY_dt",level_number,NDIM,manage_data);

            static const int global_index_offset = 0;
            static const int local_index_offset = 0;
            d_lag_init->initializeMassDataOnPatchLevel(
                global_index_offset, local_index_offset,
                M_data, K_data,
                hierarchy, level_number,
                init_data_time, can_be_refined, initial_time, d_lag_data_manager);
        }
    }

    // Initialize the marker data, but do so only on the coarsest level of the
    // patch hierarchy.  Marker data on finer levels is initialized via
    // refinement.
    if (initial_time && !d_mark_input_file_name.empty() && level_number == 0)
    {
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            const double* const patchDx = patch_geom->getDx();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_data = patch->getPatchData(d_mark_current_idx);
            for (size_t k = 0; k < d_mark_init_posns.size()/NDIM; ++k)
            {
                const double* const X = &d_mark_init_posns[NDIM*k];
                static const std::vector<double> U(NDIM,0.0);
                const bool patch_owns_node_at_loc =
                    ((  patchXLower[0] <= X[0])&&(X[0] < patchXUpper[0]))
#if (NDIM > 1)
                    &&((patchXLower[1] <= X[1])&&(X[1] < patchXUpper[1]))
#if (NDIM > 2)
                    &&((patchXLower[2] <= X[2])&&(X[2] < patchXUpper[2]))
#endif
#endif
                    ;
                if (patch_owns_node_at_loc)
                {
                    const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(
                        X, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                    if (!mark_data->isElement(i))
                    {
                        mark_data->appendItem(i, IBTK::LagMarker());
                    }
                    IBTK::LagMarker& new_mark = *(mark_data->getItem(i));
                    std::vector<double>& new_X = new_mark.getPositions();
                    std::vector<double>& new_U = new_mark.getVelocities();
                    std::vector<int>& new_idx = new_mark.getIndices();

                    new_X.insert(new_X.end(),X,X+NDIM);
                    new_U.insert(new_U.end(),U.begin(),U.end());
                    new_idx.push_back(k);
                }
            }
        }
    }

    // Determine the initial source/sink locations.
    if (initial_time && !d_source_strategy.isNull())
    {
        d_X_src.resize(std::max(int(d_X_src.size()),level_number+1));
        d_r_src.resize(std::max(int(d_r_src.size()),level_number+1));
        d_P_src.resize(std::max(int(d_P_src.size()),level_number+1));
        d_Q_src.resize(std::max(int(d_Q_src.size()),level_number+1));
        d_n_src.resize(std::max(int(d_n_src.size()),level_number+1),0);

        d_n_src[level_number] = d_source_strategy->getNumSources(hierarchy, level_number, d_integrator_time, d_lag_data_manager);

        d_X_src[level_number].resize(d_n_src[level_number], std::vector<double>(NDIM,std::numeric_limits<double>::quiet_NaN()));
        d_r_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
        d_P_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());
        d_Q_src[level_number].resize(d_n_src[level_number], std::numeric_limits<double>::quiet_NaN());

        if (d_n_src[level_number] > 0)
        {
            d_source_strategy->getSourceLocations(
                d_X_src[level_number], d_r_src[level_number],
                (d_lag_data_manager->levelContainsLagrangianData(level_number)
                 ? d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,level_number)
                 : SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>(NULL)),
                hierarchy, level_number, d_integrator_time, d_lag_data_manager);
        }
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // We use the INSHierarchyIntegrator to handle as much structured data
    // management as possible.
    d_ins_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_lag_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Indicate that the velocity field needs to be re-interpolated (but only in
    // the multi-level case).
    d_reinterpolate_after_regrid = d_reinterpolate_after_regrid || (finest_level>0);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the source/sink data vectors.
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    d_force_current_rscheds.resize(finest_hier_level+1);
    d_force_new_rscheds    .resize(finest_hier_level+1);
    if (!d_source_strategy.isNull())
    {
        d_source_rscheds.resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // Prune duplicate markers following regridding.
    pruneDuplicateMarkers(0,finest_hier_level);
    const bool initial_time = SAMRAI::tbox::MathUtilities<double>::equalEps(d_integrator_time,d_start_time);
    if (initial_time)
    {
        const unsigned int num_marks = countMarkers(0,d_hierarchy->getFinestLevelNumber());
        if (num_marks != d_mark_init_posns.size()/NDIM)
        {
            TBOX_ERROR(d_object_name << "::resetHierarchyConfiguration()\n"
                       << "  number of marker particles at initial time is incorrect\n"
                       << "  expected number of markers = " << d_mark_init_posns.size()/NDIM << "\n"
                       << "  actual   number of markers = " << num_marks << "\n");
        }
    }

    // (Re)build generic refine communication schedules.  These are created for
    // all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it!= d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[(*it).first][ln] = (*it).second->createSchedule(
                level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build specialized refine communication schedules used to compute the
    // Cartesian force and source densities.  These are set only for levels >=
    // 1.
    for (int ln = std::max(1,coarsest_level); ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        d_force_current_rscheds[ln] = d_force_current_ralg->createSchedule(
            level, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(),
            ln-1, hierarchy, d_force_current_rstrategy);

        d_force_new_rscheds[ln] = d_force_new_ralg->createSchedule(
            level, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(),
            ln-1, hierarchy, d_force_new_rstrategy);

        if (!d_source_strategy.isNull())
        {
            d_source_rscheds[ln] = d_source_ralg->createSchedule(
                level, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(),
                ln-1, hierarchy, d_source_rstrategy);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(
                coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBHierarchyIntegrator::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // It is necessary to untag all cells prior to tagging.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells for refinement according to the criteria specified by the
    // INSHierarchyIntegrator.
    d_ins_hier_integrator->applyGradientDetector(
        hierarchy, level_number, error_data_time,
        tag_index, initial_time,
        uses_richardson_extrapolation_too);

    // Tag cells which contain Lagrangian nodes.
    d_lag_data_manager->applyGradientDetector(
        hierarchy, level_number, error_data_time,
        tag_index, initial_time,
        uses_richardson_extrapolation_too);

    // Tag cells for refinement where the Cartesian source/sink strength is
    // nonzero.
    if (!d_source_strategy.isNull() && !initial_time &&
        hierarchy->finerLevelExists(level_number))
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

        const SAMRAI::hier::Index<NDIM>& lower = grid_geom->getPhysicalDomain()[0].lower();
        const SAMRAI::hier::Index<NDIM>& upper = grid_geom->getPhysicalDomain()[0].upper();
        const double* const xLower = grid_geom->getXLower();
        const double* const xUpper = grid_geom->getXUpper();
        const double* const dx = grid_geom->getDx();

        const int finer_level_number = level_number+1;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = hierarchy->getPatchLevel(finer_level_number);
        for (int n = 0; n < d_n_src[finer_level_number]; ++n)
        {
            double dx_finer[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                dx_finer[d] = dx[d]/double(finer_level->getRatio()(d));
            }

            // The source radius must be an integer multiple of the grid
            // spacing.
            double r[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                r[d] = floor(d_r_src[finer_level_number][n]/dx_finer[d])*dx_finer[d];
                r[d] = std::max(r[d],2.0*dx_finer[d]);
            }

            // Determine the approximate source stencil box.
            const SAMRAI::hier::Index<NDIM> i_center = IBTK::IndexUtilities::getCellIndex(
                d_X_src[finer_level_number][n], xLower, xUpper, dx_finer, lower, upper);
            SAMRAI::hier::Box<NDIM> stencil_box(i_center,i_center);
            for (int d = 0; d < NDIM; ++d)
            {
                stencil_box.grow(d, int(ceil(r[d]/dx_finer[d])));
            }

            const SAMRAI::hier::Box<NDIM> coarsened_stencil_box = SAMRAI::hier::Box<NDIM>::coarsen(
                stencil_box, finer_level->getRatioToCoarserLevel());

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                tags_data->fillAll(1, coarsened_stencil_box);
            }
        }
    }

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getLagMarkerVar()
///
///  allows access to the various state variables maintained by the integrator.
///

SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<NDIM,IBTK::LagMarker> >
IBHierarchyIntegrator::getLagMarkerVar() const
{
    return d_mark_var;
}// getLagMarkerVar

///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getOldContext(),
///      getScratchContext(),
///      getPlotContext()
///
///  allow access to the various variable contexts maintained by the integrator.
///

///
/// We simply reuse the SAMRAI::hier::VariableContext objects defined in the
/// INSHierarchyIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getOldContext() const
{
    return d_ins_hier_integrator->getOldContext();
}// getOldContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBHierarchyIntegrator::getPlotContext() const
{
    return d_ins_hier_integrator->getPlotContext();
}// getPlotContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
IBHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_HIERARCHY_INTEGRATOR_VERSION",
                   IB_HIERARCHY_INTEGRATOR_VERSION);

    db->putString("d_delta_fcn", d_delta_fcn);
    db->putInteger("d_total_flow_volume_sz", d_total_flow_volume.size());
    if (!d_total_flow_volume.empty())
    {
        db->putDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    db->putBool("d_using_pIB_method", d_using_pIB_method);
    db->putDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putDouble("d_old_dt", d_old_dt);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);
    db->putDouble("d_dt_max", d_dt_max);
    db->putDouble("d_dt_max_time_max", d_dt_max_time_max);
    db->putDouble("d_dt_max_time_min", d_dt_max_time_min);

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    db->putInteger("finest_hier_level", finest_hier_level);
    db->putIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->putDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            db->putDouble("d_r_src_"+id_string, d_r_src[ln][n]);
            db->putDouble("d_P_src_"+id_string, d_P_src[ln][n]);
            db->putDouble("d_Q_src_"+id_string, d_Q_src[ln][n]);
        }
    }

    const std::vector<std::string> instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (!instrument_names.empty())
    {
        db->putInteger("instrument_names_sz", instrument_names.size());
        db->putStringArray("instrument_names", &instrument_names[0], instrument_names.size());
    }

    t_put_to_database->stop();
    return;
}// putToDatabase

///
///  The following routines:
///
///      printClassData()
///
///  are provided for your viewing pleasure.
///

void
IBHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nIBHierarchyIntegrator::printClassData...\n";
    os << "this = " << const_cast<IBHierarchyIntegrator*>(this) << "\n";
    os << "d_object_name = " << d_object_name << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart << "\n";
    os << "d_delta_fcn = " << d_delta_fcn << "\n"
       << "d_ghosts = " << d_ghosts << "\n";
    os << "d_hierarchy = " << d_hierarchy.getPointer() << "\n"
       << "d_gridding_alg = " << d_gridding_alg.getPointer() << "\n";
    os << "d_visit_writer = " << d_visit_writer.getPointer() << "\n"
       << "d_silo_writer = " << d_silo_writer.getPointer() << "\n";
    os << "d_load_balancer = " << d_load_balancer.getPointer() << "\n";
    os << "d_ins_hier_integrator = " << d_ins_hier_integrator.getPointer() << "\n";
    os << "d_lag_data_manager = " << d_lag_data_manager << "\n";
    os << "d_lag_init = " << d_lag_init.getPointer() << "\n";
    os << "d_body_force_set = " << d_body_force_set.getPointer() << "\n"
       << "d_eluerian_force_set = " << d_eulerian_force_set.getPointer() << "\n"
       << "d_force_strategy = " << d_force_strategy.getPointer() << "\n"
       << "d_force_strategy_needs_init = " << d_force_strategy_needs_init << "\n";
    os << "d_eluerian_source_set = " << d_eulerian_source_set.getPointer() << "\n"
       << "d_source_strategy = " << d_source_strategy.getPointer() << "\n"
       << "d_source_strategy_needs_init = " << d_source_strategy_needs_init << "\n";
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            os << "d_X_src_" << ln << "_" << n << " = ";
            std::copy(d_X_src[ln][n].begin(), d_X_src[ln][n].end(), std::ostream_iterator<double>(os," , "));
            os << " ]\n"
               << "d_r_src_" << ln << "_" << n << d_r_src[ln][n] << "\n"
               << "d_P_src_" << ln << "_" << n << d_P_src[ln][n] << "\n"
               << "d_Q_src_" << ln << "_" << n << d_Q_src[ln][n] << "\n";
        }
        os << "d_n_src_" << ln << " = " << d_n_src[ln] << "\n";
    }
    os << "d_using_pIB_method = " << d_using_pIB_method << "\n"
       << "d_gravitational_acceleration = [ ";
    std::copy(d_gravitational_acceleration.begin(), d_gravitational_acceleration.end(), std::ostream_iterator<double>(os," , "));
    os << " ]\n";
    os << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_grow_dt = " << d_grow_dt << "\n"
       << "d_max_integrator_steps = " << d_max_integrator_steps << "\n";
    os << "d_num_cycles = " << d_num_cycles << "\n";
    os << "d_num_init_cycles = " << d_num_init_cycles << "\n";
    os << "d_regrid_interval = " << d_regrid_interval << "\n";
    os << "d_old_dt = " << d_old_dt << "\n"
       << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_integrator_step = " << d_integrator_step << "\n";
    os << "d_dt_max = " << d_dt_max << "\n"
       << "d_dt_max_time_max = " << d_dt_max_time_max << "\n"
       << "d_dt_max_time_min = " << d_dt_max_time_min << "\n";
    os << "d_is_initialized = " << d_is_initialized << "\n";
    os << "d_do_log = " << d_do_log << "\n";
    os << "d_reinterpolate_after_regrid = " << d_reinterpolate_after_regrid << "\n";
    os << "d_hier_cc_data_ops = " << d_hier_cc_data_ops.getPointer() << "\n";
    os << "Skipping variables, patch data descriptors, communications algorithms, etc.\n";
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_force_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_lag_data_manager);
        }
    }

    return;
}// resetLagrangianForceStrategy

void
IBHierarchyIntegrator::resetLagrangianSourceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_source_strategy->initializeLevelData(
                d_hierarchy, ln, init_data_time, initial_time,
                d_lag_data_manager);
        }
    }

    return;
}// resetLagrangianSourceStrategy

void
IBHierarchyIntegrator::updateIBInstrumentationData(
    const int timestep_num,
    const double data_time)
{
    if (!d_instrument_panel->isInstrumented()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute the positions of the flow meter nets.
    d_instrument_panel->initializeHierarchyDependentData(
        d_hierarchy, d_lag_data_manager, timestep_num, data_time);

    // Compute the flow rates and pressures.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getVelocityVar(),
        d_ins_hier_integrator->getScratchContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getScratchContext());

    std::vector<bool> deallocate_U_scratch_data(finest_ln+1,false);
    std::vector<bool> deallocate_P_scratch_data(finest_ln+1,false);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(U_scratch_idx))
        {
            deallocate_U_scratch_data[ln] = true;
            level->allocatePatchData(U_scratch_idx, data_time);
        }
        if (!level->checkAllocated(P_scratch_idx))
        {
            deallocate_P_scratch_data[ln] = true;
            level->allocatePatchData(P_scratch_idx, data_time);
        }
        d_rscheds["INSTRUMENTATION_DATA_FILL"][ln]->fillData(data_time);
    }

    d_instrument_panel->readInstrumentData(
        U_scratch_idx, P_scratch_idx,
        d_hierarchy, d_lag_data_manager,
        timestep_num, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (deallocate_U_scratch_data[ln]) level->deallocatePatchData(U_scratch_idx);
        if (deallocate_P_scratch_data[ln]) level->deallocatePatchData(P_scratch_idx);
    }
    return;
}// updateIBInstrumentationData

void
IBHierarchyIntegrator::resetMarkersOnPatchHierarchy()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset the assignment of markers to Cartesian grid cells on each level of
    // the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > mark_level_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    mark_level_fill_alg->registerRefine(d_mark_current_idx, // destination
                                        d_mark_current_idx, // source
                                        d_mark_scratch_idx, // temporary work space
                                        new IBTK::LagMarkerRefineOperator());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
        mark_level_fill_alg->createSchedule(level, ln-1, d_hierarchy, NULL)->fillData(d_integrator_time);
        level->deallocatePatchData(d_mark_scratch_idx);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            const double* const patchDx = patch_geom->getDx();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > current_mark_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > new_mark_data = new SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>(
                current_mark_data->getBox(), current_mark_data->getGhostCellWidth());

            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*current_mark_data); it; it++)
            {
                const IBTK::LagMarker& old_mark = it();
                const std::vector<double>& old_X = old_mark.getPositions();
                const std::vector<double>& old_U = old_mark.getVelocities();
                const std::vector<int>& old_idx = old_mark.getIndices();
                const SAMRAI::hier::IntVector<NDIM>& offset = old_mark.getPeriodicOffset();
                double X_shifted[NDIM];

                for (int k = 0; k < old_mark.getNumberOfMarkers(); ++k)
                {
                    const double* const X = &old_X[NDIM*k];
                    const double* const U = &old_U[NDIM*k];
                    const int& idx = old_idx[k];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X_shifted[d] = X[d] + double(offset(d))*patchDx[d];
                    }

                    const bool patch_owns_node_at_new_loc =
                        ((  patchXLower[0] <= X_shifted[0])&&(X_shifted[0] < patchXUpper[0]))
#if (NDIM > 1)
                        &&((patchXLower[1] <= X_shifted[1])&&(X_shifted[1] < patchXUpper[1]))
#if (NDIM > 2)
                        &&((patchXLower[2] <= X_shifted[2])&&(X_shifted[2] < patchXUpper[2]))
#endif
#endif
                        ;

                    if (patch_owns_node_at_new_loc)
                    {
                        const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(
                            X_shifted, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                        if (!new_mark_data->isElement(i))
                        {
                            new_mark_data->appendItem(i, IBTK::LagMarker());
                        }

                        IBTK::LagMarker& new_mark = *(new_mark_data->getItem(i));
                        std::vector<double>& new_X = new_mark.getPositions();
                        std::vector<double>& new_U = new_mark.getVelocities();
                        std::vector<int>& new_idx = new_mark.getIndices();

                        new_X.insert(new_X.end(),X_shifted,X_shifted+NDIM);
                        new_U.insert(new_U.end(),U,U+NDIM);
                        new_idx.push_back(idx);
                    }
                }
            }

            // Swap the old and new patch data pointers.
            patch->setPatchData(d_mark_current_idx, new_mark_data);
        }
    }
    return;
}// resetMarkersOnPatchHierarchy

void
IBHierarchyIntegrator::collectMarkersOnPatchHierarchy()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Collect all marker data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > mark_coarsen_alg = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    mark_coarsen_alg->registerCoarsen(d_mark_scratch_idx, // destination
                                      d_mark_current_idx, // source
                                      new IBTK::LagMarkerCoarsenOperator());
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln-1);

        // Allocate scratch data.
        coarser_level->allocatePatchData(d_mark_scratch_idx);

        // Coarsen fine data onto coarser level.
        SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* mark_coarsen_op = NULL;
        mark_coarsen_alg->createSchedule(coarser_level, level, mark_coarsen_op)->coarsenData();

        // Merge the coarsened fine data with the coarse data.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(coarser_level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = coarser_level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_current_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_scratch_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (!mark_current_data->isElement(i))
                {
                    mark_current_data->appendItem(i,IBTK::LagMarker());
                }

                const IBTK::LagMarker& src_mark = it();
                const std::vector<double>& src_X = src_mark.getPositions();
                const std::vector<double>& src_U = src_mark.getVelocities();
                const std::vector<int>& src_idx = src_mark.getIndices();

                IBTK::LagMarker& dst_mark = *(mark_current_data->getItem(i));
                std::vector<double>& dst_X = dst_mark.getPositions();
                std::vector<double>& dst_U = dst_mark.getVelocities();
                std::vector<int>& dst_idx = dst_mark.getIndices();

                dst_X.insert(dst_X.end(),src_X.begin(),src_X.end());
                dst_U.insert(dst_U.end(),src_U.begin(),src_U.end());
                dst_idx.insert(dst_idx.end(),src_idx.begin(),src_idx.end());
            }
        }

        // Deallocate scratch data.
        coarser_level->deallocatePatchData(d_mark_scratch_idx);
    }
    return;
}// collectMarkersOnPatchHierarchy

void
IBHierarchyIntegrator::pruneDuplicateMarkers(
    const int coarsest_ln,
    const int finest_ln)
{
    const int finest_hier_level_number = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        const bool at_finest_hier_level = ln == finest_hier_level_number;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = (at_finest_hier_level ? SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL) : d_hierarchy->getPatchLevel(ln+1));
        SAMRAI::hier::BoxArray<NDIM> refined_region_boxes = (at_finest_hier_level ? SAMRAI::hier::BoxArray<NDIM>() : finer_level->getBoxes());
        const SAMRAI::hier::IntVector<NDIM>& ratio = (at_finest_hier_level ? SAMRAI::hier::IntVector<NDIM>(1) : finer_level->getRatioToCoarserLevel());
        refined_region_boxes.coarsen(ratio);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_data = patch->getPatchData(d_mark_current_idx);
            const SAMRAI::hier::Box<NDIM>& ghost_box = mark_data->getGhostBox();
            for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
            {
                const SAMRAI::hier::Box<NDIM>& refined_box = refined_region_boxes[i];
                const SAMRAI::hier::Box<NDIM> intersection = ghost_box * refined_box;
                if (!intersection.empty())
                {
                    mark_data->removeInsideBox(intersection);
                }
            }
        }
    }
    return;
}// pruneDuplicateMarkers

int
IBHierarchyIntegrator::countMarkers(
    const int coarsest_ln,
    const int finest_ln)
{
    std::vector<int> num_marks(finest_ln+1,0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker> > mark_data = patch->getPatchData(d_mark_current_idx);
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker>::Iterator it(*mark_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (patch_box.contains(i))
                {
                    const IBTK::LagMarker& mark = it();
                    num_marks[ln] += mark.getNumberOfMarkers();
                }
            }
        }
    }
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&num_marks[0],finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::plog << "number of markers on level " << ln << ": " << num_marks[ln] << "\n";
    }
    return std::accumulate(num_marks.begin(), num_marks.end(), 0);
}// countMarkers

void
IBHierarchyIntegrator::computeConstraintForceDataStructures(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    const int coarsest_ln,
    const int finest_ln,
    const double data_time,
    const bool initial_time)
{
    (void) data_time;

    int ierr;

    // Compute the PETSc VecScatters required to compute the constraint penalty
    // forces.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_using_constraint_forces[ln])
        {
            if (d_constraint_force_src_is[ln] != static_cast<IS>(NULL))
            {
                ierr = ISDestroy(d_constraint_force_src_is[ln]); IBTK_CHKERRQ(ierr);
            }
            if (d_constraint_force_dst_is[ln] != static_cast<IS>(NULL))
            {
                ierr = ISDestroy(d_constraint_force_dst_is[ln]); IBTK_CHKERRQ(ierr);
            }
            if (d_constraint_force_vec_scatter[ln] != static_cast<VecScatter>(NULL))
            {
                ierr = VecScatterDestroy(d_constraint_force_vec_scatter[ln]); IBTK_CHKERRQ(ierr);
            }

            d_constraint_petsc_coarse_idxs[ln] = d_constraint_lag_coarse_idxs[ln];
            d_constraint_petsc_fine_idxs  [ln] = d_constraint_lag_fine_idxs  [ln];

            d_lag_data_manager->mapLagrangianToPETSc(d_constraint_petsc_coarse_idxs[ln], ln  );
            d_lag_data_manager->mapLagrangianToPETSc(d_constraint_petsc_fine_idxs  [ln], ln+1);

            Vec X_coarse_vec = X_data[ln]->getGlobalVec();
            Vec X_fine_vec = X_data[ln+1]->getGlobalVec();

            int coarse_idx_lo, coarse_idx_hi;
            ierr = VecGetOwnershipRange(X_coarse_vec, &coarse_idx_lo, &coarse_idx_hi); IBTK_CHKERRQ(ierr);
            coarse_idx_lo /= NDIM;
            coarse_idx_hi /= NDIM;

            std::vector<int> src_idxs, dst_idxs;
            for (int k = 0; k < int(d_constraint_petsc_coarse_idxs[ln].size()); ++k)
            {
                const int& coarse_idx = d_constraint_petsc_coarse_idxs[ln][k];
                const int&   fine_idx = d_constraint_petsc_fine_idxs  [ln][k];
                if (coarse_idx >= coarse_idx_lo && coarse_idx < coarse_idx_hi)
                {
                    src_idxs.push_back(NDIM*  fine_idx);
                    dst_idxs.push_back(NDIM*coarse_idx);
                }
            }
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(int(d_constraint_petsc_coarse_idxs[ln].size()) == SAMRAI::tbox::SAMRAI_MPI::sumReduction(int(dst_idxs.size())));
            TBOX_ASSERT(int(d_constraint_petsc_fine_idxs  [ln].size()) == SAMRAI::tbox::SAMRAI_MPI::sumReduction(int(src_idxs.size())));
#endif
            ierr = ISCreateBlock(PETSC_COMM_WORLD, NDIM, src_idxs.size(), &src_idxs[0], &d_constraint_force_src_is[ln]); IBTK_CHKERRQ(ierr);
            ierr = ISCreateBlock(PETSC_COMM_WORLD, NDIM, dst_idxs.size(), &dst_idxs[0], &d_constraint_force_dst_is[ln]); IBTK_CHKERRQ(ierr);
            ierr = VecScatterCreate(X_fine_vec  , d_constraint_force_src_is[ln],
                                    X_coarse_vec, d_constraint_force_dst_is[ln],
                                    &d_constraint_force_vec_scatter[ln]); IBTK_CHKERRQ(ierr);
        }
    }

    // Ensure the initial constrained positions agree, assuming that the fine
    // grid positions are the "true" positions.
    if (initial_time && d_reset_constrained_initial_posns)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_using_constraint_forces[ln])
            {
                Vec X_coarse_vec = X_data[ln]->getGlobalVec();
                Vec X_fine_vec = X_data[ln+1]->getGlobalVec();

                ierr = VecScatterBegin(d_constraint_force_vec_scatter[ln],
                                       X_fine_vec, X_coarse_vec,
                                       INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
                ierr = VecScatterEnd(d_constraint_force_vec_scatter[ln],
                                     X_fine_vec, X_coarse_vec,
                                     INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
            }
        }
    }
    return;
}// computeConstraintForceDataStructures

void
IBHierarchyIntegrator::computeConstraintForces(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data,
    const int coarsest_ln,
    const int finest_ln,
    const double data_time,
    const bool initial_time)
{
    (void) data_time;
    (void) initial_time;

    int ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_using_constraint_forces[ln])
        {
            Vec X_coarse_vec = X_data[ln]->getGlobalVec();
            Vec F_coarse_vec = F_data[ln]->getGlobalVec();
            Vec X_fine_vec = X_data[ln+1]->getGlobalVec();
            Vec F_fine_vec = F_data[ln+1]->getGlobalVec();

            int coarse_idx_lo, coarse_idx_hi;
            ierr = VecGetOwnershipRange(X_coarse_vec, &coarse_idx_lo, &coarse_idx_hi); IBTK_CHKERRQ(ierr);
            coarse_idx_lo /= NDIM;
            coarse_idx_hi /= NDIM;

            Vec X_coarsened_fine_vec;
            ierr = VecDuplicate(X_coarse_vec, &X_coarsened_fine_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecSet(X_coarsened_fine_vec, 0.0);  IBTK_CHKERRQ(ierr);

            ierr = VecScatterBegin(d_constraint_force_vec_scatter[ln],
                                   X_fine_vec, X_coarsened_fine_vec,
                                   INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(d_constraint_force_vec_scatter[ln],
                                 X_fine_vec, X_coarsened_fine_vec,
                                 INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);

            Vec F_coarse_constraint_vec;
            ierr = VecDuplicate(F_coarse_vec, &F_coarse_constraint_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecSet(F_coarse_constraint_vec, 0.0);  IBTK_CHKERRQ(ierr);

            double* X_coarse_arr;
            ierr = VecGetArray(X_coarse_vec, &X_coarse_arr);  IBTK_CHKERRQ(ierr);
            double* X_coarsened_fine_arr;
            ierr = VecGetArray(X_coarsened_fine_vec, &X_coarsened_fine_arr);  IBTK_CHKERRQ(ierr);
            double* F_coarse_constraint_arr;
            ierr = VecGetArray(F_coarse_constraint_vec, &F_coarse_constraint_arr);  IBTK_CHKERRQ(ierr);

            double r_max = 0.0;
            for (int k = 0; k < int(d_constraint_petsc_coarse_idxs[ln].size()); ++k)
            {
                const int& coarse_idx = d_constraint_petsc_coarse_idxs[ln][k];
                if (coarse_idx >= coarse_idx_lo && coarse_idx < coarse_idx_hi)
                {
                    const int local_idx = coarse_idx-coarse_idx_lo;
                    const double* const X_coarse         = &X_coarse_arr        [NDIM*local_idx];
                    const double* const X_coarsened_fine = &X_coarsened_fine_arr[NDIM*local_idx];
                    double* const F_coarse_constraint = &F_coarse_constraint_arr[NDIM*local_idx];
                    double D[NDIM];
                    double r_sq = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        D[d] = X_coarsened_fine[d] - X_coarse[d];
                        r_sq += D[d]*D[d];
                        F_coarse_constraint[d] = d_constraint_kappa*D[d];
                    }
                    r_max = std::max(r_max,sqrt(r_sq));
                }
            }
            r_max = SAMRAI::tbox::SAMRAI_MPI::maxReduction(r_max);

            if (d_do_log)
            {
                SAMRAI::tbox::plog << d_object_name << "::computeConstraintForces():\n";
                SAMRAI::tbox::plog << "  maximum displacement between constrained points on levels " << ln << " and " << ln+1 << " = " << r_max << "\n";
            }

            ierr = VecRestoreArray(X_coarse_vec, &X_coarse_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(X_coarsened_fine_vec, &X_coarsened_fine_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArray(F_coarse_constraint_vec, &F_coarse_constraint_arr);  IBTK_CHKERRQ(ierr);

            Vec F_fine_constraint_vec;
            ierr = VecDuplicate(F_fine_vec, &F_fine_constraint_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecSet(F_fine_constraint_vec, 0.0);  IBTK_CHKERRQ(ierr);

            ierr = VecScatterBegin(d_constraint_force_vec_scatter[ln],
                                   F_coarse_constraint_vec, F_fine_constraint_vec,
                                   INSERT_VALUES, SCATTER_REVERSE);  IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(d_constraint_force_vec_scatter[ln],
                                 F_coarse_constraint_vec, F_fine_constraint_vec,
                                 INSERT_VALUES, SCATTER_REVERSE);  IBTK_CHKERRQ(ierr);

            ierr = VecAXPY(F_coarse_vec, +1.0, F_coarse_constraint_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(F_fine_vec, -1.0, F_fine_constraint_vec);  IBTK_CHKERRQ(ierr);

            ierr = VecDestroy(X_coarsened_fine_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(F_coarse_constraint_vec);  IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(F_fine_constraint_vec);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// computeConstraintForces

void
IBHierarchyIntegrator::computeSourceStrengths(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> >& X_data)
{
    if (d_source_strategy.isNull()) return;

    // Reset the values of Q_src.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        d_Q_src[ln] = std::vector<double>(d_n_src[ln],0.0);
    }

    // Get the present source locations.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source locations and radii.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source locations on level number " << ln << "\n";

            d_source_strategy->getSourceLocations(
                d_X_src[ln], d_r_src[ln], X_data[ln],
                d_hierarchy, ln, data_time, d_lag_data_manager);
        }
    }

    // Compute the source strengths for each of the sources/sinks.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source strengths.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing fluid source strengths on level number " << ln << "\n";

            d_source_strategy->computeSourceStrengths(
                d_Q_src[ln], d_hierarchy, ln, data_time, d_lag_data_manager);
        }
    }

    // Spread the sources/sinks onto the Cartesian grid.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // On the coarsest level of the patch hierarchy, simply initialize the
        // Cartesian source density to equal zero.
        //
        // For each of the finer levels in the patch hierarchy, conservatively
        // interpolate the Cartesian source density from coarser levels in the
        // patch hierarchy.
        if (ln == coarsest_level)
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch->getPatchData(d_Q_idx);
                q_data->fillAll(0.0);
            }
        }
        else
        {
            // Interpolate the Cartesian source density from the next coarser
            // level in the hierarchy.
            d_source_rscheds[ln]->fillData(data_time);
        }

        if (d_n_src[ln] > 0)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading fluid source strengths to the Cartesian grid on level number " << ln << "\n";

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
                const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                const double* const xLower = pgeom->getXLower();
                const double* const xUpper = pgeom->getXUpper();
                const double* const dx = pgeom->getDx();

                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch->getPatchData(d_Q_idx);
                for (int n = 0; n < d_n_src[ln]; ++n)
                {
                    // The source radius must be an integer multiple of the grid
                    // spacing.
                    double r[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        r[d] = floor(d_r_src[ln][n]/dx[d])*dx[d];
                        r[d] = std::max(r[d],2.0*dx[d]);
                    }

                    // Determine the approximate source stencil box.
                    const SAMRAI::hier::Index<NDIM> i_center = IBTK::IndexUtilities::getCellIndex(
                        d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    SAMRAI::hier::Box<NDIM> stencil_box(i_center,i_center);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, int(ceil(r[d]/dx[d])));
                    }

                    // Spread the source strength onto the Cartesian grid.
                    for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                    {
                        const SAMRAI::hier::Index<NDIM>& i = b();

                        double wgt = 1.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            const double X_center = xLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                            wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d]);
                        }
                        (*q_data)(i) += d_Q_src[ln][n]*wgt;
                    }
                }
            }
        }
    }

    // Compute the net inflow into the computational domain.
    const int wgt_idx = d_ins_hier_integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
    SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
    double Q_sum = 0.0;
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        Q_sum = std::accumulate(d_Q_src[ln].begin(), d_Q_src[ln].end(), Q_sum);
    }
    const double q_total = d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx);

    if (d_do_log)
    {
        SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy():\n";
#if (NDIM == 2)
        SAMRAI::tbox::plog << "    Sum_{i,j} q_{i,j} h^2     = " << q_total << "\n"
                           << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum << "\n";
#endif
#if (NDIM == 3)
        SAMRAI::tbox::plog << "    Sum_{i,j,k} q_{i,j,k} h^3 = " << q_total << "\n"
                           << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum <<  "\n";
#endif
    }

    if (!SAMRAI::tbox::MathUtilities<double>::equalEps(q_total, Q_sum))
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Lagrangian and Eulerian source/sink strengths are inconsistent.");
    }

    // Balance the net inflow with outflow along the upper/lower boundaries.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const SAMRAI::hier::Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];
    const double* const dx_coarsest = grid_geom->getDx();

    SAMRAI::hier::Box<NDIM> interior_box = domain_box;
    for (int d = 0; d < NDIM-1; ++d)
    {
        interior_box.grow(d,-1);
    }

    SAMRAI::hier::BoxList<NDIM> bdry_boxes;
    bdry_boxes.removeIntersections(domain_box,interior_box);
    double vol = double(bdry_boxes.getTotalSizeOfBoxes());
    for (int d = 0; d < NDIM; ++d)
    {
        vol *= dx_coarsest[d];
    }

    const double q_norm = -q_total/vol;
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::hier::BoxList<NDIM> level_bdry_boxes(bdry_boxes);
        level_bdry_boxes.refine(level->getRatio());
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch->getPatchData(d_Q_idx);
            for (SAMRAI::hier::BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
            {
                for (SAMRAI::hier::Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                {
                    (*q_data)(b()) += q_norm;
                }
            }
        }
    }

    SAMRAI::tbox::plog << "    q_integral = " << d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx) << "\n";

    // Synchronize the Cartesian grid source density on the patch hierarchy.
    for (int ln = finest_level; ln > coarsest_level; --ln)
    {
        d_cscheds["Q->Q::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    return;
}// computeSourceStrengths

void
IBHierarchyIntegrator::computeSourcePressures(
    const int coarsest_level,
    const int finest_level,
    const double data_time,
    const std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> >& X_data)
{
    if (d_source_strategy.isNull()) return;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int P_new_idx = var_db->mapVariableAndContextToIndex(
        d_ins_hier_integrator->getPressureVar(),
        d_ins_hier_integrator->getNewContext());
    const int wgt_idx = d_ins_hier_integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();

    // Compute the normalization pressure.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const SAMRAI::hier::Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    SAMRAI::hier::Box<NDIM> interior_box = domain_box;
    for (int d = 0; d < NDIM-1; ++d)
    {
        interior_box.grow(d,-1);
    }

    SAMRAI::hier::BoxList<NDIM> bdry_boxes;
    bdry_boxes.removeIntersections(domain_box,interior_box);

    double p_norm = 0.0;
    double vol = 0.0;
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::hier::BoxList<NDIM> level_bdry_boxes(bdry_boxes);
        level_bdry_boxes.refine(level->getRatio());
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > p_data = patch->getPatchData(P_new_idx);
            const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > wgt_data = patch->getPatchData(wgt_idx);
            for (SAMRAI::hier::BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
            {
                for (SAMRAI::hier::Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = b();
                    p_norm += (*p_data)(i)*(*wgt_data)(i);
                    vol += (*wgt_data)(i);
                }
            }
        }
    }

    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&p_norm,1);
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&vol,1);

    p_norm /= vol;

    // Reset the values of P_src.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        d_P_src[ln] = std::vector<double>(d_n_src[ln],0.0);
    }

    // Get the present source locations.
    //
    // IMPORTANT NOTE: Here, we require that each MPI process is assigned the
    // same source locations and radii.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source locations on level number " << ln << "\n";

            d_source_strategy->getSourceLocations(
                d_X_src[ln], d_r_src[ln], X_data[ln],
                d_hierarchy, ln, data_time, d_lag_data_manager);
        }
    }

    // Compute the mean pressure at the sources/sinks associated with each level
    // of the Cartesian grid.
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): computing source pressures on level number " << ln << "\n";

            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
                const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                const double* const xLower = pgeom->getXLower();
                const double* const xUpper = pgeom->getXUpper();
                const double* const dx = pgeom->getDx();

                const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > p_data = patch->getPatchData(P_new_idx);
                for (int n = 0; n < d_n_src[ln]; ++n)
                {
                    // The source radius must be an integer multiple of the grid
                    // spacing.
                    double r[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        r[d] = floor(d_r_src[ln][n]/dx[d])*dx[d];
                        r[d] = std::max(r[d],2.0*dx[d]);
                    }

                    // Determine the approximate source stencil box.
                    const SAMRAI::hier::Index<NDIM> i_center = IBTK::IndexUtilities::getCellIndex(
                        d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                    SAMRAI::hier::Box<NDIM> stencil_box(i_center,i_center);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        stencil_box.grow(d, int(ceil(r[d]/dx[d])));
                    }

                    // Interpolate the pressure from the Cartesian grid.
                    for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                    {
                        const SAMRAI::hier::Index<NDIM>& i = b();

                        double wgt = 1.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            const double X_center = xLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                            wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d])*dx[d];
                        }
                        d_P_src[ln][n] += (*p_data)(i)*wgt;
                    }
                }
            }

            SAMRAI::tbox::SAMRAI_MPI::sumReduction(&d_P_src[ln][0], d_P_src[ln].size());
            std::transform(d_P_src[ln].begin(), d_P_src[ln].end(), d_P_src[ln].begin(),
                           std::bind2nd(std::plus<double>(),-p_norm));

            // Set the pressures for the Lagrangian source strategy.
            d_source_strategy->setSourcePressures(
                d_P_src[ln], d_hierarchy, ln, data_time, d_lag_data_manager);
        }
    }
    return;
}// computeSourcePressures

void
IBHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault(
        "dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault(
        "dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        d_delta_fcn = db->getStringWithDefault("delta_fcn", d_delta_fcn);

        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);

        d_num_init_cycles = db->getIntegerWithDefault(
            "num_init_cycles", d_num_init_cycles);

        d_using_pIB_method = db->getBoolWithDefault(
            "using_pIB_method", d_using_pIB_method);

        if (d_using_pIB_method)
        {
            if (db->keyExists("gravitational_acceleration"))
            {
                db->getDoubleArray("gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
            }
            else
            {
                TBOX_WARNING(d_object_name << ":  "
                             << "Using penalty-IB method but key data `gravitational_acceleration' not found in input.");
            }
        }

        d_mark_input_file_name = db->getStringWithDefault("marker_input_file_name", d_mark_input_file_name);
    }

    d_constraint_forces_file_name = db->getStringWithDefault("constraint_forces_file_name", d_constraint_forces_file_name);
    if (!d_constraint_forces_file_name.empty())
    {
        d_constraint_kappa = db->getDouble("constraint_kappa");  // XXXX
        d_reset_constrained_initial_posns = db->getBool("reset_constrained_initial_posns");  // XXXX
    }
    return;
}// getFromInput

void
IBHierarchyIntegrator::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db = SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("IB_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_delta_fcn = db->getString("d_delta_fcn");
    const int total_flow_volume_sz = db->getInteger("d_total_flow_volume_sz");
    d_total_flow_volume.resize(total_flow_volume_sz, std::numeric_limits<double>::quiet_NaN());
    if (!d_total_flow_volume.empty())
    {
        db->getDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    d_using_pIB_method = db->getBool("d_using_pIB_method");
    db->getDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_old_dt = db->getDouble("d_old_dt");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");
    d_dt_max = db->getDouble("d_dt_max");
    d_dt_max_time_max = db->getDouble("d_dt_max_time_max");
    d_dt_max_time_min = db->getDouble("d_dt_max_time_min");

    const int finest_hier_level = db->getInteger("finest_hier_level");
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);

    db->getIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        d_X_src[ln].resize(d_n_src[ln],std::vector<double>(NDIM,std::numeric_limits<double>::quiet_NaN()));
        d_r_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_P_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_Q_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->getDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            d_r_src[ln][n] = db->getDouble("d_r_src_"+id_string);
            d_P_src[ln][n] = db->getDouble("d_P_src_"+id_string);
            d_Q_src[ln][n] = db->getDouble("d_Q_src_"+id_string);
        }
    }

    if (db->keyExists("instrument_names"))
    {
        const int sz = db->getInteger("instrument_names_sz");
        std::vector<std::string> instrument_names(sz);
        db->getStringArray("instrument_names", &instrument_names[0], sz);
        IBInstrumentationSpec::setInstrumentNames(instrument_names);
    }

    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
