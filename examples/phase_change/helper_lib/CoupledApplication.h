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

#ifndef included_PhaseChangeExamples_CoupledApplication
#define included_PhaseChangeExamples_CoupledApplication

#include <ibamr/config.h>

#include <ibamr/AllenCahnHierarchyIntegrator.h>
#include <ibamr/EnthalpyHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/PhaseChangeUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>

#include <ibtk/AppInitializer.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LSLocateInterface.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <memory>
#include <string>
#include <vector>

#include "ExampleOutput.h"
#include "HeavisideFromLevelSet.h"
#include "LiquidFractionForceMask.h"

namespace PhaseChangeExamples
{
// Own registered callbacks and boundary data through integration. Each example
// explicitly selects its physical setup and registration order before run().
class CoupledApplication
{
public:
    explicit CoupledApplication(SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer);
    virtual ~CoupledApplication();
    void run();

protected:
    void useEnthalpy();
    void useAllenCahn();
    void registerLevelSet();
    void registerLiquidFraction();
    void registerLiquidFractionGradient();
    void registerExtrapolatedLiquidFraction();
    void registerLevelSetForExtrapolation();
    void registerSpecificEnthalpy();
    void registerHeaviside();
    void registerTransportFields();
    virtual void setThermalInitialConditions();
    void setTemperatureInitialCondition();
    void setLiquidFractionInitialCondition();
    void registerMaterialFields();
    void registerLevelSetTagging();
    void registerLiquidFractionTagging();
    void setHeavisideBoundary();
    void setTemperatureBoundary();
    void setEnthalpyBoundary();
    void setLiquidFractionBoundary();
    void setSpecificHeatBoundary();
    void setConductivityBoundary();
    void setVelocityBoundary();
    void setDensityBoundary();
    void setViscosityBoundary();
    void setLevelSetBoundary();
    void registerMaterialProperties();
    void registerMaterialProperties(double kappa_solid, double specific_heat_solid, double rho_solid, double mu_solid);
    void registerPhaseChangeSources(const std::string& source_name);
    void registerSolidDrag();
    void registerGravity();
    void registerGravityAndSurfaceTension(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> mask_var);
    void registerForceMask(SAMRAI::tbox::Pointer<IBAMR::SurfaceTensionForceFunction> force,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> mask_var);
    virtual void initializeDiagnostics(double time);
    virtual void postprocessStep(double time);
    virtual void finalizeDiagnostics();

    SAMRAI::tbox::Pointer<IBTK::AppInitializer> d_app_initializer;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_input_db;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_time_integrator;
    SAMRAI::tbox::Pointer<IBAMR::PhaseChangeHierarchyIntegrator> d_phase_change_integrator;
    SAMRAI::tbox::Pointer<IBAMR::EnthalpyHierarchyIntegrator> d_enthalpy_integrator;
    SAMRAI::tbox::Pointer<IBAMR::AllenCahnHierarchyIntegrator> d_allen_cahn_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> d_grid_geometry;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_patch_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_level_set_var, d_liquid_fraction_var,
        d_liquid_fraction_gradient_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_enthalpy_var, d_heaviside_var, d_temperature_var,
        d_extrapolated_liquid_fraction_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_density_var, d_specific_heat_var, d_viscosity_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_side_density_var;
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> d_liquid_fraction_bc_coef, d_temperature_bc_coef;
    double d_rho_liquid = 0.0, d_rho_solid = 0.0, d_rho_gas = 0.0, d_kappa_liquid = 0.0, d_mu_liquid = 0.0;

private:
    void initializeCoupling();
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> createScalarBoundary(const std::string& object_name,
                                                                                  const std::string& database_name);
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> createGravityForce();
    SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitialize<NDIM>> d_error_detector;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::BergerRigoutsos<NDIM>> d_box_generator;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM>> d_load_balancer;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> d_gridding_algorithm;
    SAMRAI::tbox::Pointer<IBAMR::RelaxationLSMethod> d_level_set_ops;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_level_set_initial;
    std::unique_ptr<MultiphaseExamples::LSLocateInterface> d_locator;
    std::unique_ptr<IBAMR::LevelSetUtilities::SetLSProperties> d_level_set_properties;
    std::unique_ptr<HeavisideFromLevelSet> d_heaviside_sync;
    std::unique_ptr<IBAMR::LevelSetUtilities::TagLSRefinementCells> d_level_set_tagger;
    std::unique_ptr<LiquidFractionForceMask> d_force_mask;
    std::unique_ptr<IBAMR::PhaseChangeUtilities::TagLiquidFractionRefinementCells> d_liquid_fraction_tagger;
    std::unique_ptr<IBAMR::PhaseChangeUtilities::SetFluidProperties> d_fluid_properties;
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> d_heaviside_bc_coef, d_enthalpy_bc_coef;
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> d_specific_heat_bc_coef, d_conductivity_bc_coef,
        d_density_bc_coef, d_viscosity_bc_coef, d_level_set_bc_coef;
    std::vector<std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>>> d_velocity_bc_coefs;
    ExampleOutput d_output;
};

using ApplicationFactory = std::unique_ptr<CoupledApplication> (*)(SAMRAI::tbox::Pointer<IBTK::AppInitializer>);
// One initialization/finalization owner for the coupled example entry points.
int run_coupled(int argc, char* argv[], ApplicationFactory factory);
} // namespace PhaseChangeExamples

#endif
