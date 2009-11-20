// Filename: IBStaggeredHierarchyIntegrator.C
// Last modified: <20.Nov.2009 11:21:54 griffith@boyce-griffiths-mac-pro.local>
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBStaggeredHierarchyIntegrator.h"

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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/IBInstrumentationSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeIndexData2.h>
#include <ibtk/LagMarkerCoarsen.h>
#include <ibtk/LagMarkerRefine.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <IndexData.h>
#include <Patch.h>
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

// Version of IBStaggeredHierarchyIntegrator restart file data.
static const int IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;

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
    if (std::abs(x) > eps)
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

IBStaggeredHierarchyIntegrator::IBStaggeredHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> force_strategy,
    SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> source_strategy,
    SAMRAI::tbox::Pointer<IBDataPostProcessor> post_processor,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_interp_delta_fcn("IB_4"),
      d_spread_delta_fcn("IB_4"),
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
      d_lag_init(NULL),
      d_body_force_setter(NULL),
      d_eulerian_force_setter(NULL),
      d_force_strategy(force_strategy),
      d_force_strategy_needs_init(true),
      d_eulerian_source_setter(NULL),
      d_source_strategy(source_strategy),
      d_source_strategy_needs_init(true),
      d_post_processor(post_processor),
      d_post_processor_needs_init(true),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_num_cycles(1),
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
      d_hier_cc_data_ops(),
      d_hier_sc_data_ops(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_V_var(NULL),
      d_F_var(NULL),
      d_mark_var(NULL),
      d_current(NULL),
      d_scratch(NULL),
      d_V_idx(-1),
      d_F_idx(-1),
      d_mark_current_idx(-1),
      d_mark_scratch_idx(-1)
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

    // Check the choices for the delta function.
    if (d_interp_delta_fcn != d_spread_delta_fcn)
    {
        SAMRAI::tbox::pout << "WARNING: different delta functions are being used for velocity interpolation and force spreading.\n"
                           << "         recommended usage is to employ the same delta functions for both interpolation and spreading.\n";
    }

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
                                TBOX_ERROR(d_object_name << "::IBStaggeredHierarchyIntegrator():\n"
                                           << "  encountered marker intersecting lower physical boundary.\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                            else if (X[d] <= grid_xLower[d])
                            {
                                TBOX_ERROR(d_object_name << "::IBStaggeredHierarchyIntegrator():\n"
                                           << "  encountered marker below lower physical boundary\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }

                            if (SAMRAI::tbox::MathUtilities<double>::equalEps(X[d],grid_xUpper[d]))
                            {
                                TBOX_ERROR(d_object_name << "::IBStaggeredHierarchyIntegrator():\n"
                                           << "  encountered marker intersecting upper physical boundary.\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                            else if (X[d] >= grid_xUpper[d])
                            {
                                TBOX_ERROR(d_object_name << "::IBStaggeredHierarchyIntegrator():\n"
                                           << "  encountered marker above upper physical boundary\n"
                                           << "  please ensure that all markers are within the computational domain."<< std::endl);
                            }
                        }
                    }
                }
            }
        }
    }

    // Determine the ghost cell width required for side-centered spreading and
    // interpolating.
    const int stencil_size = std::max(IBTK::LEInteractor::getStencilSize(d_interp_delta_fcn),
                                      IBTK::LEInteractor::getStencilSize(d_spread_delta_fcn));
    d_ghosts = int(floor(0.5*double(stencil_size)))+1;

    // Get the Lagrangian Data Manager.
    d_lag_data_manager = IBTK::LDataManager::getManager(d_object_name+"::LDataManager", d_ghosts, d_registered_for_restart);

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(d_object_name+"::IBInstrumentPanel", (input_db->isDatabase("IBInstrumentPanel") ? input_db->getDatabase("IBInstrumentPanel") : SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL)));

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager = SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var = new SAMRAI::pdat::SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy);

    // Initialize all variable contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_current = var_db->getContext(d_object_name+"::CURRENT");
    d_scratch = var_db->getContext(d_object_name+"::SCRATCH");

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy            = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy               = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy                = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::regridHierarchy()");
        t_synchronize_hierarchy           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels          = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_data       = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_data_to_preadvance_state  = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration   = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector         = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                 = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBStaggeredHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// IBStaggeredHierarchyIntegrator

IBStaggeredHierarchyIntegrator::~IBStaggeredHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin(); it != d_rstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin(); it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }
    return;
}// ~IBStaggeredHierarchyIntegrator

const std::string&
IBStaggeredHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
IBStaggeredHierarchyIntegrator::registerVelocityInitialConditions(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> U_init)
{
    d_ins_hier_integrator->registerVelocityInitialConditions(U_init);
    return;
}// registerVelocityInitialConditions

void
IBStaggeredHierarchyIntegrator::registerVelocityPhysicalBcCoefs(
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
    d_ins_hier_integrator->registerVelocityPhysicalBcCoefs(U_bc_coefs);
    return;
}// registerVelocityPhysicalBcCoefs

void
IBStaggeredHierarchyIntegrator::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> body_force_setter)
{
    d_body_force_setter = body_force_setter;
    return;
}// registerBodyForceSpecification

void
IBStaggeredHierarchyIntegrator::registerLNodeInitStrategy(
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
IBStaggeredHierarchyIntegrator::freeLNodeInitStrategy()
{
    d_lag_init.setNull();
    d_lag_data_manager->freeLNodeInitStrategy();
    return;
}// freeLNodeInitStrategy

void
IBStaggeredHierarchyIntegrator::registerVisItDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLagSiloDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLagM3DDataWriter(
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
IBStaggeredHierarchyIntegrator::registerLoadBalancer(
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
///      postProcessData(),
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
///  allow the IBStaggeredHierarchyIntegrator to be used as a hierarchy integrator.
///

void
IBStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Initialize all variables.
    const SAMRAI::hier::IntVector<NDIM> ghosts = d_ghosts;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_V_var = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::V");
    d_V_idx = var_db->registerVariableAndContext(d_V_var, d_scratch, ghosts);

    d_F_var = new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::F");
    d_F_idx = var_db->registerVariableAndContext(d_F_var, d_scratch, no_ghosts);

    if (!d_source_strategy.isNull())
    {
        d_Q_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::Q",1);
        d_Q_idx = var_db->registerVariableAndContext(d_Q_var, d_scratch, no_ghosts);
    }

    d_mark_var = new SAMRAI::pdat::IndexVariable<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >(d_object_name+"::mark");
    d_mark_current_idx = var_db->registerVariableAndContext(d_mark_var, getCurrentContext(), ghosts);
    d_mark_scratch_idx = var_db->registerVariableAndContext(d_mark_var, getScratchContext(), ghosts);
    if (d_registered_for_restart)
    {
        var_db->registerPatchDataForRestart(d_mark_current_idx);
    }

    // Initialize the objects used to manage Lagragian-Eulerian interaction.
    d_eulerian_force_setter = new IBEulerianForceSetter(d_object_name+"::IBEulerianForceSetter", -1, -1, d_F_idx);
    d_ins_hier_integrator->registerBodyForceSpecification(d_eulerian_force_setter);

    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_setter = new IBEulerianSourceSetter(d_object_name+"::IBEulerianSourceSetter", d_Q_idx, d_Q_idx, d_Q_idx);
        d_ins_hier_integrator->registerSourceSpecification(d_eulerian_source_setter);
    }

    // Initialize the INSStaggeredHierarchyIntegrator.
    //
    // NOTE: This must be done after all variables created by the
    // IBStaggeredHierarchyIntegrator class have been registered.
    d_ins_hier_integrator->initializeHierarchyIntegrator(d_gridding_alg);

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator;

    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getScratchContext());
    const int P_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getNewContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getScratchContext());

    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(d_V_idx, U_current_idx, d_V_idx, refine_operator);
    d_rstrategies["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"] = new IBTK::CartExtrapPhysBdryOp(d_V_idx, "QUADRATIC");

    d_ralgs["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"]->registerRefine(d_V_idx, d_V_idx, d_V_idx, refine_operator);

    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    coarsen_operator = grid_geom->lookupCoarsenOperator(d_V_var, "CONSERVATIVE_COARSEN");
    d_calgs["U->U::C->C::CONSERVATIVE_COARSEN"]->registerCoarsen(U_current_idx, U_current_idx, coarsen_operator);

    d_ralgs["INSTRUMENTATION_DATA_FILL"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getVelocityVar(), "CONSERVATIVE_LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(U_scratch_idx, U_new_idx, U_scratch_idx, refine_operator);
    refine_operator = grid_geom->lookupRefineOperator(d_ins_hier_integrator->getPressureVar(), "LINEAR_REFINE");
    d_ralgs["INSTRUMENTATION_DATA_FILL"]->registerRefine(P_scratch_idx, P_new_idx, P_scratch_idx, refine_operator);
    SAMRAI::hier::ComponentSelector instrumentation_data_fill_bc_idxs;
    instrumentation_data_fill_bc_idxs.setFlag(U_scratch_idx);
    instrumentation_data_fill_bc_idxs.setFlag(P_scratch_idx);
    d_rstrategies["INSTRUMENTATION_DATA_FILL"] = new IBTK::CartExtrapPhysBdryOp(instrumentation_data_fill_bc_idxs, "QUADRATIC");

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
IBStaggeredHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    // Use the INSStaggeredHierarchyIntegrator to initialize the patch
    // hierarchy.
    double dt_next = d_ins_hier_integrator->initializeHierarchy();

    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
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
    d_lag_data_manager->updateWorkloadData(coarsest_ln, finest_ln);

    // Prune duplicate markers following initialization.
    pruneDuplicateMarkers(0,d_hierarchy->getFinestLevelNumber());

    // Ensure that we haven't misplaced any of the markers.
    const int num_marks = d_mark_init_posns.size()/NDIM;
    const int num_marks_after_init = countMarkers(0,d_hierarchy->getFinestLevelNumber(),true);
    if (num_marks != num_marks_after_init)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy()\n"
                   << "  number of marker particles is incorrect\n"
                   << "  expected number of markers = " << num_marks << "\n"
                   << "  actual   number of markers = " << num_marks_after_init << "\n");
    }

    // Initialize the instrumentation data.
    d_instrument_panel->initializeHierarchyIndependentData(d_hierarchy, d_lag_data_manager);
    if (d_instrument_panel->isInstrumented())
    {
        d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_lag_data_manager, d_integrator_step, d_integrator_time);
        if (d_total_flow_volume.empty())
        {
            d_total_flow_volume.resize(d_instrument_panel->getFlowValues().size(),0.0);
        }
    }

    // Indicate that the force strategy and post processor need to be
    // re-initialized.
    d_force_strategy_needs_init = true;
    d_source_strategy_needs_init = true;
    d_post_processor_needs_init = true;

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
IBStaggeredHierarchyIntegrator::advanceHierarchy(
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
    const bool do_regrid = (d_regrid_interval == 0 ? false : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid)
    {
        if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): regridding prior to timestep " << d_integrator_step << "\n";
        regridHierarchy();
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Set the current time interval in the force and (optional) source
    // specification objects.
    d_eulerian_force_setter->registerBodyForceSpecification(d_body_force_setter);
    d_eulerian_force_setter->setTimeInterval(current_time, new_time);
    d_force_strategy->setTimeInterval(current_time, new_time);
    if (!d_source_strategy.isNull())
    {
        d_eulerian_source_setter->setTimeInterval(current_time, new_time);
        d_source_strategy->setTimeInterval(current_time, new_time);
    }

    // (Re)initialize the force and (optional) source strategies and the
    // post-processor.
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
    if (d_post_processor_needs_init && !d_post_processor.isNull())
    {
        resetPostProcessor(current_time, initial_time);
        d_post_processor_needs_init = false;
    }

    // Allocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_V_idx, current_time);
        level->allocatePatchData(d_F_idx, current_time);
        if (!d_source_strategy.isNull())
        {
            level->allocatePatchData(d_Q_idx, current_time);
        }
    }

    // Initialize the various LNodeLevelData objects on each level of the patch hierarchy.
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_new_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_half_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_half_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_half_data(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::VEL_DATA_NAME,ln);
            X_new_data[ln] = d_lag_data_manager->createLNodeLevelData("X_new",ln,NDIM);
            X_half_data[ln] = d_lag_data_manager->createLNodeLevelData("X_half",ln,NDIM);
            U_half_data[ln] = d_lag_data_manager->createLNodeLevelData("U_half",ln,NDIM);
            F_half_data[ln] = d_lag_data_manager->createLNodeLevelData("F_half",ln,NDIM);

            X_data[ln]->restoreLocalFormVec();
            U_data[ln]->restoreLocalFormVec();
            X_new_data[ln]->restoreLocalFormVec();
            X_half_data[ln]->restoreLocalFormVec();
            U_half_data[ln]->restoreLocalFormVec();
            F_half_data[ln]->restoreLocalFormVec();
        }
    }

    // Get patch data descriptors for the current and new velocity data.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int U_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getNewContext());

    // Synchronize the Cartesian grid velocity u(n) on the patch hierarchy, then
    // interpolate the velocity field.
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["U->U::C->C::CONSERVATIVE_COARSEN"][ln]->coarsenData();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_rscheds["U->V::C->S::CONSERVATIVE_LINEAR_REFINE"][ln]->fillData(current_time);
    }
    interp(U_data, d_V_idx, false, X_data, false, coarsest_ln, finest_ln);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Initialize X(n+1/2) to equal X(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_half_vec = X_half_data[ln]->getGlobalVec();
            int ierr = VecCopy(X_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize X(n+1) to equal X(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            int ierr = VecCopy(X_vec, X_new_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize X_mark(n+1) to equal X_mark(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_scratch_data->copy(*mark_data);
        }
    }

    // Initialize U(n+1/2) to equal U(n).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec U_vec = U_data[ln]->getGlobalVec();
            Vec U_half_vec = U_half_data[ln]->getGlobalVec();
            int ierr = VecCopy(U_vec, U_half_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Perform one or more cycles of fixed point iteration to compute the
    // updated configuration of the coupled fluid-structure system.
    d_ins_hier_integrator->integrateHierarchy_initialize(current_time, new_time);
    for (int cycle = 0; cycle < d_num_cycles; ++cycle)
    {
        // Compute F(n+1/2) = F(X(n+1/2),U(n+1/2),t_{n+1/2}).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                Vec F_half_vec = F_half_data[ln]->getGlobalVec();
                int ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
                d_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], d_hierarchy, ln, current_time+0.5*dt, d_lag_data_manager);
            }
        }
        resetAnchorPointValues(F_half_data, coarsest_ln, finest_ln);

        // Spread F(n+1/2) to f(n+1/2).
        spread(d_F_idx, F_half_data, true, X_half_data, true, coarsest_ln, finest_ln);

        // Compute the source/sink strengths corresponding to any distributed
        // internal fluid sources or sinks.
        if (!d_source_strategy.isNull())
        {
            computeSourceStrengths(coarsest_ln, finest_ln, current_time+0.5*dt, X_half_data);
        }

        // Solve the incompressible Navier-Stokes equations.
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time);

        // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
        d_hier_sc_data_ops->linearSum(d_V_idx, 0.5, U_current_idx, 0.5, U_new_idx);

        // Interpolate u(n+1/2) to U(n+1/2).
        interp(U_half_data, d_V_idx, true, X_half_data, false, coarsest_ln, finest_ln);
        resetAnchorPointValues(U_half_data, coarsest_ln, finest_ln);

        // Set X(n+1) = X(n) + dt*U(n+1/2).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                int ierr;
                Vec X_vec = X_data[ln]->getGlobalVec();
                Vec X_new_vec = X_new_data[ln]->getGlobalVec();
                Vec U_half_vec = U_half_data[ln]->getGlobalVec();
                ierr = VecWAXPY(X_new_vec, dt, U_half_vec, X_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Set X_mark(n+1) = X_mark(n) + dt*U(n+1/2).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > v_data = patch->getPatchData(d_V_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);

                // Collect into a single vector the current positions of all
                // markers in the patch.
                std::vector<double> X_mark;
#ifdef DEBUG_CHECK_ASSERTIONS
                std::vector<int> idx_mark;
#endif
                for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*mark_data); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                    if (patch_box.contains(i))
                    {
                        const IBTK::LagMarker& mark = it();
                        const std::vector<double>& X = mark.getPositions();
                        X_mark.insert(X_mark.end(), X.begin(), X.end());
#ifdef DEBUG_CHECK_ASSERTIONS
                        const std::vector<int>& idx = mark.getIndices();
                        idx_mark.insert(idx_mark.end(), idx.begin(), idx.end());
#endif
                    }
                }

                // Collect into a single vector the predicted new positions of
                // all markers in the patch.
                std::vector<double> X_mark_new;
#ifdef DEBUG_CHECK_ASSERTIONS
                std::vector<int> idx_mark_new;
#endif
                for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*mark_scratch_data); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                    if (patch_box.contains(i))
                    {
                        const IBTK::LagMarker& mark_new = it();
                        const std::vector<double>& X_new = mark_new.getPositions();
                        X_mark_new.insert(X_mark_new.end(), X_new.begin(), X_new.end());
#ifdef DEBUG_CHECK_ASSERTIONS
                        const std::vector<int>& idx_new = mark_new.getIndices();
                        idx_mark_new.insert(idx_mark_new.end(), idx_new.begin(), idx_new.end());
#endif
                    }
                }

                // Compute X_mark(n+1/2) = 0.5*(X_mark(n+1) + X_mark(n)).
                //
                // NOTE: It is important here that mark_scratch_data is
                // initialized as a copy of mark_data.  Otherwise, the ordering
                // of the current and new markers are not guaranteed to be the
                // same.
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(X_mark.size() == X_mark_new.size());
                TBOX_ASSERT(idx_mark.size() == idx_mark_new.size());
                for (unsigned k = 0; k < idx_mark.size(); ++k)
                {
                    TBOX_ASSERT(idx_mark[k] == idx_mark_new[k]);
                }
#endif
                std::vector<double> X_mark_half(X_mark.size());
                for (unsigned k = 0; k < X_mark.size(); ++k)
                {
                    X_mark_half[k] = 0.5*(X_mark_new[k]+X_mark[k]);
                }

                // Compute U_mark(n+1/2) = u(X_mark(n+1/2),n+1/2).
                std::vector<double> U_mark_half(X_mark.size());
                IBTK::LEInteractor::interpolate(U_mark_half, NDIM, X_mark_half, NDIM, v_data, patch, patch_box, d_interp_delta_fcn);

                // Compute X_mark(n+1) = X_mark(n) + dt*U_mark(n+1/2).
                for (unsigned k = 0; k < X_mark.size(); ++k)
                {
                    X_mark_new[k] = X_mark[k] + dt*U_mark_half[k];
                }

                // Prevent markers from leaving the computational domain through
                // physical boundaries (but *not* through periodic boundaries).
                static const double edge_tol = 1.0e-6;
                const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
                if (periodic_shift.min() == 0)
                {
                    const double* const xLower = grid_geom->getXLower();
                    const double* const xUpper = grid_geom->getXUpper();
                    for (unsigned k = 0; k < X_mark_new.size()/NDIM; ++k)
                    {
                        double* const X = &X_mark_new[NDIM*k];
                        for (int d = 0; d < NDIM; ++d)
                        {
                            if (periodic_shift[d] == 0)
                            {
                                X[d] = std::max(X[d],xLower[d]+edge_tol);
                                X[d] = std::min(X[d],xUpper[d]-edge_tol);
                            }
                        }
                    }
                }

                // Store the updated marker positions.
                int marker_offset = 0;
                for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*mark_scratch_data); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                    if (patch_box.contains(i))
                    {
                        IBTK::LagMarker& mark = it();
                        const int nmarks = mark.getNumberOfMarkers();
                        std::vector<double> X(X_mark_new.begin()+NDIM*marker_offset,X_mark_new.begin()+NDIM*(marker_offset+nmarks));
                        mark.setPositions(X);
                        marker_offset += nmarks;
                    }
                }
            }
        }

        // Set X(n+1/2) = 0.5*(X(n) + X(n+1)).
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                int ierr;
                Vec X_vec = X_data[ln]->getGlobalVec();
                Vec X_new_vec = X_new_data[ln]->getGlobalVec();
                Vec X_half_vec = X_half_data[ln]->getGlobalVec();
                ierr = VecCopy(X_vec, X_half_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecAXPBY(X_half_vec, 0.5, 0.5, X_new_vec);  IBTK_CHKERRQ(ierr);
            }
        }

        // Compute the pressure at the updated locations of any distributed
        // internal fluid sources or sinks.
        if (!d_source_strategy.isNull())
        {
            computeSourcePressures(coarsest_ln, finest_ln, current_time+0.5*dt, X_half_data);
        }
    }
    d_ins_hier_integrator->integrateHierarchy_finalize(current_time, new_time);

    // Update the instrumentation data.
    updateIBInstrumentationData(d_integrator_step+1,new_time);
    if (d_instrument_panel->isInstrumented())
    {
        const std::vector<std::string>& instrument_name = d_instrument_panel->getInstrumentNames();
        const std::vector<double>& flow_data = d_instrument_panel->getFlowValues();
        for (unsigned m = 0; m < flow_data.size(); ++m)
        {
            // NOTE: Flow volume is calculated in default units.
            d_total_flow_volume[m] += flow_data[m]*dt;
            SAMRAI::tbox::plog << "flow volume through " << instrument_name[m] << ":\t " << d_total_flow_volume[m] << "\n";
        }
    }

    // Reset X to equal X_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_new_vec = X_new_data[ln]->getGlobalVec();
            int ierr = VecCopy(X_new_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Reset X_mark to equal X_mark_new.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            mark_data->copy(*mark_scratch_data);
        }
        level->deallocatePatchData(d_mark_scratch_idx);
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_V_idx);
        level->deallocatePatchData(d_F_idx);
        if (!d_source_strategy.isNull())
        {
            level->deallocatePatchData(d_Q_idx);
        }
    }

    // Synchronize all data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): synchronizing data\n";
    synchronizeHierarchy();

    // Reset all time dependent data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): resetting time dependent data\n";
    resetTimeDependentHierData(new_time);

    // Determine the next stable timestep.
    double dt_next = d_ins_hier_integrator->getStableTimestep(getCurrentContext());

    if (d_integrator_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    if (d_integrator_time >= d_dt_max_time_min && d_integrator_time <= d_dt_max_time_max)
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

void
IBStaggeredHierarchyIntegrator::postProcessData()
{
    if (d_post_processor.isNull()) return;

    const double current_time = d_integrator_time;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Initialize X_data, F_data, and U_data on each level of the patch
    // hierarchy.
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data(finest_ln+1);
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::VEL_DATA_NAME,ln);
            F_data[ln] = d_lag_data_manager->createLNodeLevelData("F",ln,NDIM);
        }
    }

    // Interpolate u(n) from the Cartesian grid onto the Lagrangian mesh.
    interp(U_data, d_V_idx, true, X_data, true, coarsest_ln, finest_ln);
    resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Compute F(n) = F(X(n),U(n),n), the Lagrangian force corresponding to
    // configuration X(n) and velocity U(n) at time t_{n}.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec F_vec = F_data[ln]->getGlobalVec();
            int ierr = VecSet(F_vec, 0.0);  IBTK_CHKERRQ(ierr);
            d_force_strategy->computeLagrangianForce(F_data[ln], X_data[ln], U_data[ln], d_hierarchy, ln, current_time, d_lag_data_manager);
        }
    }
    resetAnchorPointValues(F_data, coarsest_ln, finest_ln);

    // Perform the user-defined post-processing.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int U_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getCurrentContext());
    const int P_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getCurrentContext());
    const int F_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getForceVar(), d_ins_hier_integrator->getCurrentContext());
    d_post_processor->postProcessData(U_current_idx, P_current_idx, F_current_idx, F_data, X_data, U_data, d_hierarchy, coarsest_ln, finest_ln, current_time, d_lag_data_manager);
}// postProcessData

bool
IBStaggeredHierarchyIntegrator::atRegridPoint() const
{
    static const int level_number = 0;
    return ((d_integrator_step>0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
IBStaggeredHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
IBStaggeredHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
IBStaggeredHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
IBStaggeredHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
IBStaggeredHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
IBStaggeredHierarchyIntegrator::stepsRemaining() const
{
    return d_integrator_step < d_max_integrator_steps;
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
IBStaggeredHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
IBStaggeredHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

IBTK::LDataManager*
IBStaggeredHierarchyIntegrator::getLDataManager() const
{
    return d_lag_data_manager;
}// getLDataManager

SAMRAI::tbox::Pointer<IBInstrumentPanel>
IBStaggeredHierarchyIntegrator::getIBInstrumentPanel() const
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
///  allow the IBStaggeredHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
IBStaggeredHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    // Update the marker data.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): resetting markers particles.\n";
    const int num_marks = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);

    collectMarkersOnPatchHierarchy();
    const int num_marks_after_collection = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);
    const int num_marks_after_collection_level_0 = countMarkers(0,0,false);
    if (num_marks != num_marks_after_collection || num_marks != num_marks_after_collection_level_0)
    {
        TBOX_ERROR(d_object_name << "::regridHierarchy()\n"
                   << "  number of marker particles changed during collection to coarsest level\n"
                   << "  number of markers in hierarchy before collection to coarsest level = " << num_marks << "\n"
                   << "  number of markers in hierarchy after  collection to coarsest level = " << num_marks_after_collection << "\n"
                   << "  number of markers on level 0   after  collection to coarsest level = " << num_marks_after_collection_level_0 << "\n");
    }

    // Update the workload pre-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Before regriding, begin Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): starting Lagrangian data movement.\n";
    d_lag_data_manager->beginDataRedistribution();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): calling INSStaggeredHierarchyIntegrator::regridHierarchy().\n";
    d_ins_hier_integrator->regridHierarchy();

    // After regridding, finish Lagrangian data movement.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): finishing Lagrangian data movement.\n";
    d_lag_data_manager->endDataRedistribution();

    // Update the workload post-regridding.
    if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::regridHierarchy(): updating workload estimates.\n";
    d_lag_data_manager->updateWorkloadData(0,d_hierarchy->getFinestLevelNumber());

    // Prune duplicate markers following regridding.
    pruneDuplicateMarkers(0,d_hierarchy->getFinestLevelNumber());

    // Ensure that we haven't misplaced any of the markers.
    const int num_marks_after_regrid = countMarkers(0,d_hierarchy->getFinestLevelNumber(),true);
    if (num_marks != num_marks_after_regrid)
    {
        TBOX_ERROR(d_object_name << "::regridHierarchy()\n"
                   << "  number of marker particles changed during regrid\n"
                   << "  number of markers before regrid = " << num_marks << "\n"
                   << "  number of markers after  regrid = " << num_marks_after_regrid << "\n");
    }

    // Indicate that the force and (optional) source strategies and
    // post-processor need to be re-initialized.
    d_force_strategy_needs_init  = true;
    d_source_strategy_needs_init = true;
    d_post_processor_needs_init  = true;

    // Lookup the re-distributed Lagrangian position data.
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data(d_hierarchy->getFinestLevelNumber()+1);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            X_data[ln] = d_lag_data_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
        }
    }

    // Compute the set of local anchor points.
    static const double eps = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        d_anchor_point_local_idxs[ln].clear();
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(lag_node_index_idx);
                for (IBTK::LNodeIndexData2::Iterator it(patch_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                    const IBTK::LNodeIndexSet& node_set = (*idx_data)(i);
                    for (IBTK::LNodeIndexSet::const_iterator n = node_set.begin();
                         n != node_set.end(); ++n)
                    {
                        const IBTK::LNodeIndexSet::value_type& node_idx = *n;
                        const std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> >& stash_data = node_idx->getStashData();
                        for (unsigned l = 0; l < stash_data.size(); ++l)
                        {
                            SAMRAI::tbox::Pointer<IBAnchorPointSpec> anchor_point_spec = stash_data[l];
                            if (!anchor_point_spec.isNull())
                            {
                                d_anchor_point_local_idxs[ln].insert(node_idx->getLocalPETScIndex());
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < X_data[ln]->getLocalNodeCount(); ++i)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    if ((periodic_shift[d] == 0) &&
                        ((*X_data[ln])(i,d) - grid_xLower[d] <= eps ||
                         grid_xUpper[d] - (*X_data[ln])(i,d) <= eps))
                    {
                        d_anchor_point_local_idxs[ln].insert(i);
                        break;
                    }
                }
            }
        }
    }

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

void
IBStaggeredHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeHierarchy();

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
IBStaggeredHierarchyIntegrator::synchronizeNewLevels(
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
    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->synchronizeNewLevels(
        hierarchy, coarsest_level, finest_level,
        sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
IBStaggeredHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_data->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetTimeDependentHierData(new_time);

    // Reset the time of the current data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(new_time, d_mark_current_idx);
    }

    t_reset_time_dependent_data->stop();
    return;
}// resetTimeDependentHierData

void
IBStaggeredHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_data_to_preadvance_state->start();

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
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
IBStaggeredHierarchyIntegrator::initializeLevelData(
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

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_lag_data_manager->setPatchHierarchy(hierarchy);
    d_lag_data_manager->resetLevels(0,hierarchy->getFinestLevelNumber());
    d_lag_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

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

    // On the coarsest level of the patch hierarchy, copy marker data from the
    // old coarse level.  Otherwise, refine marker data from the coarsest level
    // of the patch hierarchy.
    if (!old_level.isNull() && level_number == 0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > copy_mark_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        copy_mark_alg->registerRefine(d_mark_current_idx, d_mark_current_idx, d_mark_scratch_idx, NULL);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > dst_level = level;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > src_level = old_level;
        SAMRAI::xfer::RefinePatchStrategy<NDIM>* refine_mark_op = NULL;
        level->allocatePatchData(d_mark_scratch_idx, init_data_time);
        copy_mark_alg->createSchedule(dst_level, src_level, refine_mark_op)->fillData(init_data_time);
        level->deallocatePatchData(d_mark_scratch_idx);
    }
    else if (level_number > 0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_mark_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        refine_mark_alg->registerRefine(d_mark_current_idx, d_mark_current_idx, d_mark_scratch_idx, new IBTK::LagMarkerRefine());
        SAMRAI::xfer::RefinePatchStrategy<NDIM>* refine_mark_op = NULL;
        for (int ln = 1; ln <= level_number; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > dst_level = hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_mark_scratch_idx, init_data_time);
            refine_mark_alg->createSchedule(dst_level, NULL, ln-1, hierarchy, refine_mark_op)->fillData(init_data_time);
            level->deallocatePatchData(d_mark_scratch_idx);
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

            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
            for (unsigned k = 0; k < d_mark_init_posns.size()/NDIM; ++k)
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
                    const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(X, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
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
IBStaggeredHierarchyIntegrator::resetHierarchyConfiguration(
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

    // We use the INSStaggeredHierarchyIntegrator to handle as much structured
    // data management as possible.
    d_ins_hier_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // We use the LDataManager to handle as much unstructured data management as
    // possible.
    d_lag_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

   // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the schedule vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build generic refine communication schedules.  These are created for
    // all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin(); it != d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            d_rscheds[(*it).first][ln] = (*it).second->createSchedule(level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
IBStaggeredHierarchyIntegrator::applyGradientDetector(
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
    // INSStaggeredHierarchyIntegrator.
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
    if (!d_source_strategy.isNull() && !initial_time && hierarchy->finerLevelExists(level_number))
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

SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > >
IBStaggeredHierarchyIntegrator::getLagMarkerVar() const
{
    return d_mark_var;
}// getLagMarkerVar

///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getScratchContext()
///
///  allow access to the various variable contexts maintained by the integrator.
///

///
/// We simply reuse the SAMRAI::hier::VariableContext objects defined in the
/// INSStaggeredHierarchyIntegrator object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getCurrentContext() const
{
    return d_ins_hier_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getNewContext() const
{
    return d_ins_hier_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
IBStaggeredHierarchyIntegrator::getScratchContext() const
{
    return d_ins_hier_integrator->getScratchContext();
}// getScratchContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
IBStaggeredHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);

    db->putString("d_interp_delta_fcn", d_interp_delta_fcn);
    db->putString("d_spread_delta_fcn", d_spread_delta_fcn);
    db->putInteger("d_total_flow_volume_sz", d_total_flow_volume.size());
    if (!d_total_flow_volume.empty())
    {
        db->putDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_num_cycles", d_num_cycles);
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

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStaggeredHierarchyIntegrator::spread(
    const int f_data_idx,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data,
    const bool F_data_ghost_node_update,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    const bool X_data_ghost_node_update,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (F_data_ghost_node_update) F_data[ln]->beginGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->beginGhostUpdate();
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            if (F_data_ghost_node_update) F_data[ln]->endGhostUpdate();
            if (X_data_ghost_node_update) X_data[ln]->endGhostUpdate();
        }
    }

    d_hier_sc_data_ops->setToScalar(f_data_idx, 0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            TBOX_ASSERT(ln == d_hierarchy->getFinestLevelNumber());
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_data = patch->getPatchData(f_data_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                const SAMRAI::hier::Box<NDIM>& box = idx_data->getGhostBox();
                IBTK::LEInteractor::spread(f_data, F_data[ln], X_data[ln], idx_data, patch, box, periodic_shift, d_spread_delta_fcn);
            }
        }
    }
    return;
}// spread

void
IBStaggeredHierarchyIntegrator::interp(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data,
    const int u_data_idx,
    const bool u_data_ghost_cell_update,
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
    const bool X_data_ghost_node_update,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? d_hierarchy->getFinestLevelNumber() : finest_ln_in);
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    if (u_data_ghost_cell_update)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var;
        var_db->mapIndexToVariable(u_data_idx, var);
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        refine_alg->registerRefine(u_data_idx, u_data_idx, u_data_idx, refine_op);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > refine_sched = d_rscheds["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"][ln];
            refine_alg->resetSchedule(refine_sched);
            refine_sched->fillData(d_integrator_time);
            d_ralgs["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"]->resetSchedule(refine_sched);
        }
    }

    if (X_data_ghost_node_update)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                X_data[ln]->beginGhostUpdate();
            }
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                X_data[ln]->endGhostUpdate();
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            TBOX_ASSERT(ln == d_hierarchy->getFinestLevelNumber());
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(d_lag_data_manager->getLNodeIndexPatchDescriptorIndex());
                const SAMRAI::hier::Box<NDIM>& box = idx_data->getBox();
                IBTK::LEInteractor::interpolate(U_data[ln], X_data[ln], idx_data, u_data, patch, box, periodic_shift, d_interp_delta_fcn);
            }
        }
    }
    return;
}// interp

void
IBStaggeredHierarchyIntegrator::resetLagrangianForceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_force_strategy->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_lag_data_manager);
        }
    }
    return;
}// resetLagrangianForceStrategy

void
IBStaggeredHierarchyIntegrator::resetLagrangianSourceStrategy(
    const double init_data_time,
    const bool initial_time)
{
    if (!d_source_strategy.isNull())
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            if (d_lag_data_manager->levelContainsLagrangianData(ln))
            {
                d_source_strategy->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_lag_data_manager);
            }
        }
    }
    return;
}// resetLagrangianSourceStrategy

void
IBStaggeredHierarchyIntegrator::resetPostProcessor(
    const double init_data_time,
    const bool initial_time)
{
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            d_post_processor->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_lag_data_manager);
        }
    }
    return;
}// resetPostProcessor

void
IBStaggeredHierarchyIntegrator::updateIBInstrumentationData(
    const int timestep_num,
    const double data_time)
{
    if (!d_instrument_panel->isInstrumented()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute the positions of the flow meter nets.
    d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_lag_data_manager, timestep_num, data_time);

    // Compute the flow rates and pressures.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int U_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVar(), d_ins_hier_integrator->getScratchContext());
    const int P_scratch_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVar(), d_ins_hier_integrator->getScratchContext());

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

    d_instrument_panel->readInstrumentData(U_scratch_idx, P_scratch_idx, d_hierarchy, d_lag_data_manager, timestep_num, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (deallocate_U_scratch_data[ln]) level->deallocatePatchData(U_scratch_idx);
        if (deallocate_P_scratch_data[ln]) level->deallocatePatchData(P_scratch_idx);
    }
    return;
}// updateIBInstrumentationData

void
IBStaggeredHierarchyIntegrator::collectMarkersOnPatchHierarchy()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const int num_marks_before_coarsening = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);

    // Collect all marker data on the patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > mark_coarsen_alg = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    mark_coarsen_alg->registerCoarsen(d_mark_scratch_idx, d_mark_current_idx, new IBTK::LagMarkerCoarsen());
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln-1);

        // Allocate scratch data.
        coarser_level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);

        // Coarsen fine data onto coarser level.
        SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* mark_coarsen_op = NULL;
        mark_coarsen_alg->createSchedule(coarser_level, level, mark_coarsen_op)->coarsenData();

        // Merge the coarsened fine data with the coarse data.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(coarser_level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = coarser_level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_current_data = patch->getPatchData(d_mark_current_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_scratch_data = patch->getPatchData(d_mark_scratch_idx);
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*mark_scratch_data); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
                if (!mark_current_data->isElement(i))
                {
                    mark_current_data->appendItem(i,IBTK::LagMarker());
                }
                IBTK::LagMarker& dst_mark = *(mark_current_data->getItem(i));
                const IBTK::LagMarker& src_mark = it();
                dst_mark.addMarker(src_mark);
            }
        }

        // Clear the fine data.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_current_data = patch->getPatchData(d_mark_current_idx);
            mark_current_data->removeAllItems();
        }

        // Deallocate scratch data.
        coarser_level->deallocatePatchData(d_mark_scratch_idx);
    }

    // Ensure that the total number of markers is correct.
    const int num_marks_after_coarsening = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);
    const int num_marks_after_coarsening_level_0 = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);
    if (num_marks_before_coarsening != num_marks_after_coarsening || num_marks_before_coarsening != num_marks_after_coarsening_level_0)
    {
        TBOX_ERROR(d_object_name << "::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during collection to coarsest level\n"
                   << "  number of markers in hierarchy before collection to coarsest level = " << num_marks_before_coarsening << "\n"
                   << "  number of markers in hierarchy after  collection to coarsest level = " << num_marks_after_coarsening << "\n"
                   << "  number of markers on level 0   after  collection to coarsest level = " << num_marks_after_coarsening_level_0 << "\n");
    }

    // Reset the assignment of markers to Cartesian grid cells on the coarsest
    // level of the patch hierarchy.
    //
    // NOTE: It is important to do this only *after* collecting markers on the
    // patch hierarchy, as markers which have left a fine level through the
    // coarse-fine interface would be discarded by this procedure.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > mark_level_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    mark_level_fill_alg->registerRefine(d_mark_current_idx, d_mark_current_idx, d_mark_scratch_idx, NULL);
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(coarsest_ln);
    level->allocatePatchData(d_mark_scratch_idx, d_integrator_time);
    mark_level_fill_alg->createSchedule(level,NULL)->fillData(d_integrator_time);
    level->deallocatePatchData(d_mark_scratch_idx);
    int added_counter = 0;
    int interior_counter = 0;
    int ghost_counter = 0;
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

        SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > current_mark_data = patch->getPatchData(d_mark_current_idx);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > >     new_mark_data =
            new SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >(current_mark_data->getBox(), current_mark_data->getGhostCellWidth());
        for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*current_mark_data); it; it++)
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
                    const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(X_shifted, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
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
                    added_counter++;
                }

                if (current_mark_data->getBox().contains(it.getIndex()))
                {
                    interior_counter++;
                }
                else
                {
                    ghost_counter++;
                }
            }
        }

        // Swap the old and new patch data pointers.
        patch->setPatchData(d_mark_current_idx, new_mark_data);
    }

    // Ensure that the total number of markers is correct.
    SAMRAI::tbox::plog << "added_counter = " << SAMRAI::tbox::SAMRAI_MPI::sumReduction(added_counter) << "\n";
    SAMRAI::tbox::plog << "interior_counter = " << SAMRAI::tbox::SAMRAI_MPI::sumReduction(interior_counter) << "\n";
    SAMRAI::tbox::plog << "ghost_counter = " << SAMRAI::tbox::SAMRAI_MPI::sumReduction(ghost_counter) << "\n";
    const int num_marks_after_posn_reset = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);
    const int num_marks_after_posn_reset_level_0 = countMarkers(0,d_hierarchy->getFinestLevelNumber(),false);
    if (num_marks_before_coarsening != num_marks_after_posn_reset || num_marks_before_coarsening != num_marks_after_posn_reset_level_0)
    {
        TBOX_ERROR(d_object_name << "::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during position reset on coarsest level\n"
                   << "  number of markers in hierarchy before position reset on coarsest level = " << num_marks_before_coarsening << "\n"
                   << "  number of markers in hierarchy after  position reset on coarsest level = " << num_marks_after_posn_reset << "\n"
                   << "  number of markers on level 0   after  position reset on coarsest level = " << num_marks_after_posn_reset_level_0 << "\n");
    }
    return;
}// collectMarkersOnPatchHierarchy

void
IBStaggeredHierarchyIntegrator::pruneDuplicateMarkers(
    const int coarsest_ln,
    const int finest_ln)
{
    const int finest_hier_level_number = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= std::min(finest_ln,finest_hier_level_number-1); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(ln+1);
        SAMRAI::hier::BoxArray<NDIM> refined_region_boxes = finer_level->getBoxes();
        const SAMRAI::hier::IntVector<NDIM>& ratio = finer_level->getRatioToCoarserLevel();
        refined_region_boxes.coarsen(ratio);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
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
IBStaggeredHierarchyIntegrator::countMarkers(
    const int coarsest_ln,
    const int finest_ln,
    const bool log_results)
{
    std::vector<int> num_marks(finest_ln+1,0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_current_idx);
            for (SAMRAI::pdat::IndexData<NDIM,IBTK::LagMarker,SAMRAI::pdat::CellGeometry<NDIM> >::Iterator it(*mark_data); it; it++)
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
    if (log_results)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::plog << "number of markers on level " << ln << ": " << num_marks[ln] << "\n";
        }
    }
    return std::accumulate(num_marks.begin(), num_marks.end(), 0);
}// countMarkers

void
IBStaggeredHierarchyIntegrator::resetAnchorPointValues(
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > V_data,
    const int coarsest_ln,
    const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_lag_data_manager->levelContainsLagrangianData(ln))
        {
            const int depth = V_data[ln]->getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(depth == NDIM);
#endif
            Vec V_vec = V_data[ln]->getGlobalVec();
            double* V_arr;
            int ierr = VecGetArray(V_vec, &V_arr);  IBTK_CHKERRQ(ierr);
            for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin(); cit != d_anchor_point_local_idxs[ln].end(); ++cit)
            {
                const int& i = *cit;
                for (int d = 0; d < depth; ++d)
                {
                    V_arr[depth*i+d] = 0.0;
                }
            }
            ierr = VecRestoreArray(V_vec, &V_arr);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// resetAnchorPointValues

void
IBStaggeredHierarchyIntegrator::computeSourceStrengths(
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
    d_hier_cc_data_ops->setToScalar(d_Q_idx, 0.0);
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        if (d_n_src[ln] > 0)
        {
            TBOX_ASSERT(ln == d_hierarchy->getFinestLevelNumber());
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << "::advanceHierarchy(): spreading fluid source strengths to the Cartesian grid on level number " << ln << "\n";
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
    double Q_max = 0.0;
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        Q_sum = std::accumulate(d_Q_src[ln].begin(), d_Q_src[ln].end(), Q_sum);
        for (unsigned k = 0; k < d_Q_src[ln].size(); ++k)
        {
            Q_max = std::max(Q_max,std::abs(d_Q_src[ln][k]));
        }
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

    if (std::abs(q_total-Q_sum)                     > 1.0e-12 &&
        std::abs(q_total-Q_sum)/std::max(Q_max,1.0) > 1.0e-12)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Lagrangian and Eulerian source/sink strengths are inconsistent.");
    }

    // When necessary, balance the net inflow/outflow with outflow/inflow along
    // the upper/lower boundaries of the computational domain.
    if (std::abs(q_total) > 1.0e-12)
    {
        SAMRAI::tbox::plog << "    adding ``external'' source/sink to offset net inflow/outflow into domain.\n";
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

        if (std::abs(d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx) > 1.0e-12))
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "``external'' source/sink does not correctly offset net inflow/outflow into domain.\n"
                       << "integral{q} = " << d_hier_cc_data_ops->integral(d_Q_idx, wgt_idx) << " != 0.\n");
        }
    }
    return;
}// computeSourceStrengths

void
IBStaggeredHierarchyIntegrator::computeSourcePressures(
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
IBStaggeredHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault("max_integrator_steps", d_max_integrator_steps);

    d_num_cycles = db->getIntegerWithDefault("num_cycles", d_num_cycles);

    d_regrid_interval = db->getIntegerWithDefault("regrid_interval", d_regrid_interval);

    d_dt_max = db->getDoubleWithDefault("dt_max",d_dt_max);
    d_dt_max_time_max = db->getDoubleWithDefault("dt_max_time_max", d_dt_max_time_max);
    d_dt_max_time_min = db->getDoubleWithDefault("dt_max_time_min", d_dt_max_time_min);

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    if (!is_from_restart)
    {
        if (db->isString("interp_delta_fcn") && db->isString("spread_delta_fcn"))
        {
            d_interp_delta_fcn = db->getStringWithDefault("interp_delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("spread_delta_fcn", d_spread_delta_fcn);
        }
        else
        {
            d_interp_delta_fcn = db->getStringWithDefault("delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("delta_fcn", d_spread_delta_fcn);
        }
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
        d_mark_input_file_name = db->getStringWithDefault("marker_input_file_name", d_mark_input_file_name);
    }
    return;
}// getFromInput

void
IBStaggeredHierarchyIntegrator::getFromRestart()
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

    int ver = db->getInteger("IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_interp_delta_fcn = db->getString("d_interp_delta_fcn");
    d_spread_delta_fcn = db->getString("d_spread_delta_fcn");
    const int total_flow_volume_sz = db->getInteger("d_total_flow_volume_sz");
    d_total_flow_volume.resize(total_flow_volume_sz, std::numeric_limits<double>::quiet_NaN());
    if (!d_total_flow_volume.empty())
    {
        db->getDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_num_cycles = db->getInteger("d_num_cycles");
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
template class SAMRAI::tbox::Pointer<IBAMR::IBStaggeredHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
