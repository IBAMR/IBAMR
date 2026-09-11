// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "CoupledApplication.h"
// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/CarmanKozenyDragForce.h>
#include <ibamr/EnthalpyHierarchyIntegrator.h>
#include <ibamr/HeavisideForcingFunction.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/PhaseChangeDivUSourceFunction.h>
#include <ibamr/PhaseChangeUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
CoupledApplication::CoupledApplication(Pointer<AppInitializer> app_initializer)
    : d_app_initializer(app_initializer),
      d_input_db(app_initializer->getInputDatabase()),
      d_velocity_bc_coefs(NDIM),
      d_output(app_initializer)
{
    d_time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
        "INSVCStaggeredConservativeHierarchyIntegrator",
        d_app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
}

CoupledApplication::~CoupledApplication() = default;

void
CoupledApplication::run()
{
    d_output.registerDataWriter(d_app_initializer, d_time_integrator);

    // Initialize hierarchy configuration and data on all patches.
    d_time_integrator->initializePatchHierarchy(d_patch_hierarchy, d_gridding_algorithm);

    // Remove the AppInitializer
    d_app_initializer.setNull();

    // Print the input database contents to the log file.
    plog << "Input database:\n";
    d_input_db->printClassData(plog);

    const int iteration_num = d_time_integrator->getIntegratorStep();
    const double loop_time = d_time_integrator->getIntegratorTime();
    d_output.writeInitial(d_time_integrator, d_patch_hierarchy, iteration_num, loop_time);

    initializeDiagnostics(loop_time);

    run_time_loop(d_time_integrator, d_patch_hierarchy, d_output, [this](const double time) { postprocessStep(time); });

    finalizeDiagnostics();
}

void
CoupledApplication::useEnthalpy()
{
    d_enthalpy_integrator = new EnthalpyHierarchyIntegrator(
        "EnthalpyHierarchyIntegrator", d_app_initializer->getComponentDatabase("EnthalpyHierarchyIntegrator"));
    d_phase_change_integrator = d_enthalpy_integrator;
    initializeCoupling();
}

void
CoupledApplication::useAllenCahn()
{
    d_allen_cahn_integrator = new AllenCahnHierarchyIntegrator(
        "AllenCahnHierarchyIntegrator", d_app_initializer->getComponentDatabase("AllenCahnHierarchyIntegrator"));
    d_phase_change_integrator = d_allen_cahn_integrator;
    initializeCoupling();
}

void
CoupledApplication::registerLevelSet()
{
    d_level_set_var = new CellVariable<NDIM, double>("ls_var");
    d_phase_change_integrator->registerTransportedQuantity(d_level_set_var, true);
    d_phase_change_integrator->setDiffusionCoefficient(d_level_set_var, 0.0);

    d_level_set_ops =
        new RelaxationLSMethod("RelaxationLSMethod", d_app_initializer->getComponentDatabase("RelaxationLSMethod"));
    d_level_set_initial = new muParserCartGridFunction(
        "ls_init", d_app_initializer->getComponentDatabase("LevelSetInitialConditions"), d_grid_geometry);
    d_locator = std::make_unique<MultiphaseExamples::LSLocateInterface>(
        d_phase_change_integrator, d_level_set_var, d_level_set_initial);
    d_level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&MultiphaseExamples::call_locate_interface,
                                                              static_cast<void*>(d_locator.get()));
    d_level_set_properties =
        std::make_unique<IBAMR::LevelSetUtilities::SetLSProperties>("SetLSProperties", d_level_set_ops);
    d_phase_change_integrator->registerResetFunction(d_level_set_var,
                                                     &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy,
                                                     static_cast<void*>(d_level_set_properties.get()));
}

void
CoupledApplication::registerLiquidFraction()
{
    d_liquid_fraction_var = new CellVariable<NDIM, double>("lf_var");
    d_phase_change_integrator->registerLiquidFractionVariable(d_liquid_fraction_var, true);
    d_extrapolated_liquid_fraction_var = d_liquid_fraction_var;
}

void
CoupledApplication::registerLiquidFractionGradient()
{
    d_liquid_fraction_gradient_var = new CellVariable<NDIM, double>("lf_gradient_var", NDIM);
    d_phase_change_integrator->registerLiquidFractionGradientVariable(d_liquid_fraction_gradient_var, true);
}

void
CoupledApplication::registerExtrapolatedLiquidFraction()
{
    d_extrapolated_liquid_fraction_var = new CellVariable<NDIM, double>("lf_extrap_var");
    d_enthalpy_integrator->registerLiquidFractionVariableForExtrapolation(d_extrapolated_liquid_fraction_var);
}

void
CoupledApplication::registerLevelSetForExtrapolation()
{
    d_enthalpy_integrator->registerLevelSetVariable(d_level_set_var);
}

void
CoupledApplication::registerSpecificEnthalpy()
{
    d_enthalpy_var = new CellVariable<NDIM, double>("h_var");
    d_enthalpy_integrator->registerSpecificEnthalpyVariable(d_enthalpy_var, true);
}

void
CoupledApplication::registerHeaviside()
{
    d_heaviside_var = new CellVariable<NDIM, double>("heaviside_var");
    d_phase_change_integrator->registerTransportedQuantity(d_heaviside_var, true);
    d_phase_change_integrator->setDiffusionCoefficient(d_heaviside_var, 0.0);
}

void
CoupledApplication::registerTransportFields()
{
    d_phase_change_integrator->registerHeavisideVariable(d_heaviside_var);

    d_temperature_var = new CellVariable<NDIM, double>("Temperature");
    d_phase_change_integrator->registerTemperatureVariable(d_temperature_var, true);

    d_phase_change_integrator->AdvDiffHierarchyIntegrator::setAdvectionVelocity(
        d_level_set_var, d_time_integrator->getAdvectionVelocityVariable());
    d_phase_change_integrator->AdvDiffHierarchyIntegrator::setAdvectionVelocity(
        d_heaviside_var, d_time_integrator->getAdvectionVelocityVariable());
    d_phase_change_integrator->setAdvectionVelocity(d_time_integrator->getAdvectionVelocityVariable());

    const ConvectiveDifferencingType ls_difference_form =
        IBAMR::string_to_enum<ConvectiveDifferencingType>(d_input_db->getString("LS_CONVECTIVE_FORM"));
    d_phase_change_integrator->setConvectiveDifferencingType(d_level_set_var, ls_difference_form);

    const ConvectiveDifferencingType H_difference_form =
        IBAMR::string_to_enum<ConvectiveDifferencingType>(d_input_db->getString("H_CONVECTIVE_FORM"));
    d_phase_change_integrator->setConvectiveDifferencingType(d_heaviside_var, H_difference_form);

    d_phase_change_integrator->setResetPriority(d_level_set_var, 0);
    d_phase_change_integrator->setResetPriority(d_heaviside_var, 1);

    d_phase_change_integrator->setInitialConditions(d_level_set_var, d_level_set_initial);

    // H is initialized by the level-set synchronization callback.

    setThermalInitialConditions();

    if (d_input_db->keyExists("VelocityInitialConditions"))
    {
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", d_app_initializer->getComponentDatabase("VelocityInitialConditions"), d_grid_geometry);
        d_time_integrator->registerVelocityInitialConditions(u_init);
    }

    if (d_input_db->keyExists("PressureInitialConditions"))
    {
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", d_app_initializer->getComponentDatabase("PressureInitialConditions"), d_grid_geometry);
        d_time_integrator->registerPressureInitialConditions(p_init);
    }

    d_heaviside_sync = std::make_unique<HeavisideFromLevelSet>(
        d_phase_change_integrator, d_level_set_var, d_input_db->getDouble("NUMBER_OF_INTERFACE_CELLS"));

    d_phase_change_integrator->registerResetFunction(
        d_heaviside_var,
        &PhaseChangeExamples::HeavisideFromLevelSet::synchronize_levelset_with_heaviside_fcn,
        static_cast<void*>(d_heaviside_sync.get()));
}

void
CoupledApplication::setThermalInitialConditions()
{
    setTemperatureInitialCondition();
    setLiquidFractionInitialCondition();
}

void
CoupledApplication::setTemperatureInitialCondition()
{
    Pointer<CartGridFunction> T_init = new muParserCartGridFunction(
        "T_init", d_app_initializer->getComponentDatabase("TemperatureInitialConditions"), d_grid_geometry);
    d_phase_change_integrator->setTemperatureInitialCondition(d_temperature_var, T_init);
}

void
CoupledApplication::setLiquidFractionInitialCondition()
{
    Pointer<CartGridFunction> lf_init = new muParserCartGridFunction(
        "lf_init", d_app_initializer->getComponentDatabase("LiquidFractionInitialConditions"), d_grid_geometry);
    d_phase_change_integrator->setLiquidFractionInitialCondition(d_liquid_fraction_var, lf_init);
}

void
CoupledApplication::registerMaterialFields()
{
    // Setup the INS maintained material properties.
    d_side_density_var = new SideVariable<NDIM, double>("rho_sc_var");
    d_time_integrator->registerMassDensityVariable(d_side_density_var);

    d_viscosity_var = new CellVariable<NDIM, double>("mu");
    d_time_integrator->registerViscosityVariable(d_viscosity_var);

    d_density_var = new CellVariable<NDIM, double>("rho_cc_var");
    d_phase_change_integrator->registerDensityVariable(d_density_var, true);

    d_specific_heat_var = new CellVariable<NDIM, double>("Cp");
    d_phase_change_integrator->registerSpecificHeatVariable(d_specific_heat_var, true);
}

void
CoupledApplication::registerLevelSetTagging()
{
    const double tag_thresh = d_input_db->getDouble("LS_TAG_ABS_THRESH");
    d_level_set_tagger = std::make_unique<LevelSetUtilities::TagLSRefinementCells>(
        d_phase_change_integrator, d_level_set_var, -tag_thresh, tag_thresh);
    d_time_integrator->registerApplyGradientDetectorCallback(&LevelSetUtilities::tagLSCells, d_level_set_tagger.get());
}

void
CoupledApplication::registerLiquidFractionTagging()
{
    const double min_tag_val = d_input_db->getDouble("MIN_TAG_VAL");
    const double max_tag_val = d_input_db->getDouble("MAX_TAG_VAL");
    d_liquid_fraction_tagger = std::make_unique<IBAMR::PhaseChangeUtilities::TagLiquidFractionRefinementCells>(
        d_phase_change_integrator, d_liquid_fraction_var, d_liquid_fraction_gradient_var, min_tag_val, max_tag_val);
    d_phase_change_integrator->registerApplyGradientDetectorCallback(
        &IBAMR::PhaseChangeUtilities::call_tag_liquid_fraction_cells_callback,
        static_cast<void*>(d_liquid_fraction_tagger.get()));
}

void
CoupledApplication::setHeavisideBoundary()
{
    d_heaviside_bc_coef = createScalarBoundary("H_bc_coef", "HeavisideBcCoefs");
    if (d_heaviside_bc_coef)
    {
        d_phase_change_integrator->setPhysicalBcCoef(d_heaviside_var, d_heaviside_bc_coef.get());
    }
}

void
CoupledApplication::setTemperatureBoundary()
{
    d_temperature_bc_coef = createScalarBoundary("T_bc_coef", "TemperatureBcCoefs");
    if (d_temperature_bc_coef)
    {
        d_phase_change_integrator->setTemperaturePhysicalBcCoef(d_temperature_var, d_temperature_bc_coef.get());
    }
}

void
CoupledApplication::setEnthalpyBoundary()
{
    d_enthalpy_bc_coef = createScalarBoundary("h_bc_coef", "EnthalpyBcCoefs");
    if (d_enthalpy_bc_coef)
    {
        d_enthalpy_integrator->setEnthalpyBcCoef(d_enthalpy_bc_coef.get());
    }
}

void
CoupledApplication::setLiquidFractionBoundary()
{
    d_liquid_fraction_bc_coef = createScalarBoundary("lf_bc_coef", "LiquidFractionBcCoefs");
    if (d_liquid_fraction_bc_coef && d_allen_cahn_integrator)
    {
        d_allen_cahn_integrator->setLiquidFractionPhysicalBcCoef(d_liquid_fraction_var,
                                                                 d_liquid_fraction_bc_coef.get());
    }
}

void
CoupledApplication::setSpecificHeatBoundary()
{
    d_specific_heat_bc_coef = createScalarBoundary("Cp_bc_coef", "SpecificHeatBcCoefs");
    if (d_specific_heat_bc_coef)
    {
        d_allen_cahn_integrator->registerSpecificHeatBoundaryConditions(d_specific_heat_bc_coef.get());
    }
}

void
CoupledApplication::setConductivityBoundary()
{
    d_conductivity_bc_coef = createScalarBoundary("k_bc_coef", "ThermalConductivityBcCoefs");
    if (d_conductivity_bc_coef)
    {
        d_phase_change_integrator->registerThermalConductivityBoundaryConditions(d_conductivity_bc_coef.get());
    }
}

void
CoupledApplication::setVelocityBoundary()
{
    const IntVector<NDIM>& periodic_shift = d_grid_geometry->getPeriodicShift();
    if (periodic_shift.min() == 0)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream bc_coefs_name_stream;
            bc_coefs_name_stream << "u_bc_coefs_" << d;
            const string bc_coefs_name = bc_coefs_name_stream.str();

            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

            d_velocity_bc_coefs[d] = std::make_unique<muParserRobinBcCoefs>(
                bc_coefs_name, d_app_initializer->getComponentDatabase(bc_coefs_db_name), d_grid_geometry);
        }
        d_time_integrator->registerPhysicalBoundaryConditions({
            d_velocity_bc_coefs[0].get(), d_velocity_bc_coefs[1].get()
#if (NDIM == 3)
                                              ,
                d_velocity_bc_coefs[2].get()
#endif
        });
    }
}

void
CoupledApplication::setDensityBoundary()
{
    d_density_bc_coef = createScalarBoundary("rho_bc_coef", "DensityBcCoefs");
    if (d_density_bc_coef)
    {
        d_time_integrator->registerMassDensityBoundaryConditions(d_density_bc_coef.get());
        d_phase_change_integrator->registerMassDensityBoundaryConditions(d_density_bc_coef.get());
    }
}

void
CoupledApplication::setViscosityBoundary()
{
    d_viscosity_bc_coef = createScalarBoundary("mu_bc_coef", "ViscosityBcCoefs");
    if (d_viscosity_bc_coef)
    {
        d_time_integrator->registerViscosityBoundaryConditions(d_viscosity_bc_coef.get());
    }
}

void
CoupledApplication::setLevelSetBoundary()
{
    d_level_set_bc_coef = createScalarBoundary("ls_bc_coef", "LevelSetBcCoefs");
    if (d_level_set_bc_coef)
    {
        d_phase_change_integrator->setPhysicalBcCoef(d_level_set_var, d_level_set_bc_coef.get());
        d_level_set_ops->registerPhysicalBoundaryCondition(d_level_set_bc_coef.get());
    }
}

void
CoupledApplication::registerMaterialProperties()
{
    registerMaterialProperties(d_input_db->getDouble("KAPPA_S"),
                               d_input_db->getDouble("CP_S"),
                               d_input_db->getDouble("RHO_S"),
                               d_input_db->getDouble("MU_S"));
}

void
CoupledApplication::registerMaterialProperties(const double kappa_solid,
                                               const double specific_heat_solid,
                                               const double rho_solid,
                                               const double mu_solid)
{
    d_kappa_liquid = d_input_db->getDouble("KAPPA_L");
    const double kappa_gas = d_input_db->getDouble("KAPPA_G");
    const double Cp_liquid = d_input_db->getDouble("CP_L");
    const double Cp_gas = d_input_db->getDouble("CP_G");
    d_rho_liquid = d_input_db->getDouble("RHO_L");
    d_rho_solid = rho_solid;
    d_rho_gas = d_input_db->getDouble("RHO_G");
    d_mu_liquid = d_input_db->getDouble("MU_L");
    const double mu_gas = d_input_db->getDouble("MU_G");

    // Callback functions can either be registered with the NS integrator, or
    // the advection-diffusion integrator
    d_fluid_properties =
        std::make_unique<IBAMR::PhaseChangeUtilities::SetFluidProperties>("SetFluidProperties",
                                                                          d_phase_change_integrator,
                                                                          d_heaviside_var,
                                                                          d_heaviside_bc_coef.get(),
                                                                          d_liquid_fraction_var,
                                                                          d_liquid_fraction_bc_coef.get(),
                                                                          d_rho_liquid,
                                                                          d_rho_solid,
                                                                          d_rho_gas,
                                                                          d_kappa_liquid,
                                                                          kappa_solid,
                                                                          kappa_gas,
                                                                          Cp_liquid,
                                                                          specific_heat_solid,
                                                                          Cp_gas,
                                                                          d_mu_liquid,
                                                                          mu_solid,
                                                                          mu_gas);

    d_time_integrator->registerResetFluidDensityFcn(&IBAMR::PhaseChangeUtilities::call_set_density_callback,
                                                    static_cast<void*>(d_fluid_properties.get()));
    d_time_integrator->registerResetFluidViscosityFcn(&IBAMR::PhaseChangeUtilities::call_set_viscosity_callback,
                                                      static_cast<void*>(d_fluid_properties.get()));

    d_phase_change_integrator->registerResetDiffusionCoefficientFcn(
        &IBAMR::PhaseChangeUtilities::call_set_thermal_conductivity_callback,
        static_cast<void*>(d_fluid_properties.get()));

    d_phase_change_integrator->registerResetSpecificHeatFcn(
        &IBAMR::PhaseChangeUtilities::call_set_specific_heat_callback, static_cast<void*>(d_fluid_properties.get()));

    d_phase_change_integrator->registerResetDensityFcn(&IBAMR::PhaseChangeUtilities::call_set_density_callback,
                                                       static_cast<void*>(d_fluid_properties.get()));
}

void
CoupledApplication::registerPhaseChangeSources(const std::string& source_name)
{
    // Register H Div U term in the Heaviside equation.
    Pointer<CellVariable<NDIM, double>> F_var = new CellVariable<NDIM, double>(source_name);
    d_phase_change_integrator->registerSourceTerm(F_var, true);
    Pointer<CartGridFunction> H_forcing_fcn = new HeavisideForcingFunction(
        "H_forcing_fcn", d_phase_change_integrator, d_heaviside_var, d_time_integrator->getAdvectionVelocityVariable());
    d_phase_change_integrator->setSourceTermFunction(F_var, H_forcing_fcn);
    d_phase_change_integrator->setSourceTerm(d_heaviside_var, F_var);

    // Register source term for Div U equation.
    Pointer<CartGridFunction> Div_U_forcing_fcn =
        new PhaseChangeDivUSourceFunction("Div_U_forcing_fcn", d_phase_change_integrator);
    d_time_integrator->registerVelocityDivergenceFunction(Div_U_forcing_fcn);
}

void
CoupledApplication::registerSolidDrag()
{
    // Configure the drag force object to enforce solid velocity to be zero.
    Pointer<CarmanKozenyDragForce> drag_force =
        new CarmanKozenyDragForce("drag_force",
                                  d_heaviside_var,
                                  d_liquid_fraction_var,
                                  d_phase_change_integrator,
                                  d_time_integrator,
                                  d_app_initializer->getComponentDatabase("CarmanKozenyDragForce"),
                                  /*register_for_restart*/ true);
    d_time_integrator->registerBrinkmanPenalizationStrategy(drag_force);
}

void
CoupledApplication::registerGravity()
{
    d_time_integrator->registerBodyForceFunction(createGravityForce());
}

void
CoupledApplication::registerGravityAndSurfaceTension(Pointer<CellVariable<NDIM, double>> mask_var)
{
    Pointer<SurfaceTensionForceFunction> surface_tension_force =
        new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                        d_app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                        d_phase_change_integrator,
                                        d_level_set_var);
    registerForceMask(surface_tension_force, mask_var);
    Pointer<CartGridFunction> grav_force = createGravityForce();
    Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
    eul_forces->addFunction(grav_force);
    eul_forces->addFunction(surface_tension_force);
    d_time_integrator->registerBodyForceFunction(eul_forces);
}

void
CoupledApplication::registerForceMask(Pointer<SurfaceTensionForceFunction> force,
                                      Pointer<CellVariable<NDIM, double>> mask_var)
{
    d_force_mask = std::make_unique<LiquidFractionForceMask>(mask_var,
                                                             d_liquid_fraction_bc_coef.get(),
                                                             d_phase_change_integrator,
                                                             d_time_integrator,
                                                             d_rho_liquid,
                                                             d_rho_gas);
    force->registerSurfaceTensionForceMasking(&LiquidFractionForceMask::mask_surface_tension_force, d_force_mask.get());
}

void
CoupledApplication::initializeDiagnostics(double /*time*/)
{
}

void
CoupledApplication::postprocessStep(double /*time*/)
{
}

void
CoupledApplication::finalizeDiagnostics()
{
}

void
CoupledApplication::initializeCoupling()
{
    d_time_integrator->registerAdvDiffHierarchyIntegrator(d_phase_change_integrator);
    d_grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                      d_app_initializer->getComponentDatabase("CartesianGeometry"));
    d_patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", d_grid_geometry);

    d_error_detector =
        new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                           d_time_integrator,
                                           d_app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    d_box_generator = new BergerRigoutsos<NDIM>();
    d_load_balancer = new LoadBalancer<NDIM>("LoadBalancer", d_app_initializer->getComponentDatabase("LoadBalancer"));
    d_gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                                       d_app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                       d_error_detector,
                                                       d_box_generator,
                                                       d_load_balancer);
}

std::unique_ptr<RobinBcCoefStrategy<NDIM>>
CoupledApplication::createScalarBoundary(const std::string& object_name, const std::string& database_name)
{
    if (!(d_grid_geometry->getPeriodicShift().min() > 0) && d_input_db->keyExists(database_name))
    {
        return std::make_unique<muParserRobinBcCoefs>(
            object_name, d_app_initializer->getComponentDatabase(database_name), d_grid_geometry);
    }
    return nullptr;
}

Pointer<CartGridFunction>
CoupledApplication::createGravityForce()
{
    // Register gravity force.
    std::vector<double> grav_const(NDIM);
    d_input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
    Pointer<CartGridFunction> grav_force =
        new IBAMR::VCINSUtilities::GravityForcing("GravityForcing", d_time_integrator, grav_const);
    return grav_force;
}

int
run_coupled(int argc, char* argv[], const ApplicationFactory factory)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);
    {
        auto application = factory(new AppInitializer(argc, argv, "INS.log"));
        application->run();
    }
    return 0;
}
} // namespace PhaseChangeExamples
