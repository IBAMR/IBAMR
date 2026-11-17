// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/CarreauYasudaViscosity.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
// Wraps a viscosity model to check each viscosity update. The integrator should compute the viscosity of the initial
// velocity once. Then, in each time step, it should update the viscosity at the current time, run the preprocess
// callbacks, and update the viscosity at the new time at the start of each cycle and after the last one.
class ObservedViscosity : public ViscosityModel
{
public:
    ObservedViscosity(std::unique_ptr<ViscosityModel> model,
                      const INSVCStaggeredHierarchyIntegrator& ins_integrator,
                      const double mu_min,
                      const double mu_max)
        : d_model(std::move(model)), d_ins_integrator(ins_integrator), d_mu_min(mu_min), d_mu_max(mu_max)
    {
    }

    int getVelocityGhostWidth() const override
    {
        return d_model->getVelocityGhostWidth();
    }

    void setDataOnPatchHierarchy(const int mu_idx,
                                 const int U_idx,
                                 IBTK::HierarchyMathOps& hier_math_ops,
                                 const double time) override
    {
        d_model->setDataOnPatchHierarchy(mu_idx, U_idx, hier_math_ops, time);
        ++num_updates;

        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> mu_var = d_ins_integrator.getViscosityVariable();
        if (mu_idx == var_db->mapVariableAndContextToIndex(mu_var, d_ins_integrator.getCurrentContext()))
        {
            d_events.push_back({ Event::CURRENT, time, time, 0 });
        }
        else if (mu_idx == var_db->mapVariableAndContextToIndex(mu_var, d_ins_integrator.getNewContext()))
        {
            d_events.push_back({ Event::NEW, time, time, 0 });
        }
        else
        {
            error("viscosity update of unknown patch data");
        }

        // Every updated value must lie between the infinite and zero shear viscosities.
        int num_bad_values = 0;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy = hier_math_ops.getPatchHierarchy();
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
                for (SAMRAI::pdat::CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
                {
                    const double mu = (*mu_data)(ic());
                    if (!(std::isfinite(mu) && mu >= d_mu_min && mu <= d_mu_max)) ++num_bad_values;
                }
            }
        }
        if (IBTK::IBTK_MPI::sumReduction(num_bad_values) > 0) error("viscosity out of range");
    }

    // Record the start of a time step; ctx is this object.
    static void start_step(const double current_time, const double new_time, const int num_cycles, void* ctx)
    {
        static_cast<ObservedViscosity*>(ctx)->d_events.push_back({ Event::STEP, current_time, new_time, num_cycles });
    }

    // Check the recorded sequence of updates and time steps.
    void finish()
    {
        const std::size_t num_events = d_events.size();
        if (num_events == 0 || d_events[0].type != Event::CURRENT)
        {
            error("the viscosity of the initial velocity was not computed");
            return;
        }
        for (std::size_t i = 1; i < num_events; ++num_steps)
        {
            const std::string step = " in time step " + std::to_string(num_steps);
            if (i + 1 >= num_events || d_events[i].type != Event::CURRENT || d_events[i + 1].type != Event::STEP ||
                d_events[i].time != d_events[i + 1].time)
            {
                error("no update of the viscosity at the current time" + step);
                return;
            }
            const Event& start = d_events[i + 1];
            if (num_cycles >= 0 && start.num_cycles != num_cycles) error("number of cycles changed" + step);
            num_cycles = start.num_cycles;
            i += 2;
            for (int k = 0; k <= num_cycles; ++k, ++i)
            {
                if (i >= num_events || d_events[i].type != Event::NEW || d_events[i].time != start.new_time)
                {
                    error("wrong updates of the viscosity at the new time" + step);
                    return;
                }
            }
        }
    }

    int num_updates = 0, num_steps = 0, num_cycles = -1, num_errors = 0;

private:
    struct Event
    {
        enum Type
        {
            CURRENT,
            NEW,
            STEP
        } type;
        double time, new_time;
        int num_cycles;
    };

    void error(const std::string& message)
    {
        SAMRAI::tbox::pout << "error: " << message << '\n';
        ++num_errors;
    }

    std::unique_ptr<ViscosityModel> d_model;
    const INSVCStaggeredHierarchyIntegrator& d_ins_integrator;
    double d_mu_min, d_mu_max;
    std::vector<Event> d_events;
};

// Place one passive marker point at the center of the domain on the finest level.
int finest_ln = 0;

void
generate_structure(const unsigned int& /*struct_num*/,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* /*ctx*/)
{
    num_vertices = ln == finest_ln ? 1 : 0;
    vertex_posn.assign(num_vertices, IBTK::Point::Constant(0.5));
    return;
}
} // namespace

// Advance a Carreau-Yasuda fluid with the integrator whose input database is present, check every viscosity update,
// and report the final velocity norms and viscosity range. The conservative integrator also transports a
// side-centered density. If an IBHierarchyIntegrator database is present, the fluid integrator is a child of an IB
// integrator with one passive marker point, which regrids the hierarchy.
int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender(new TestAppender());

    SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv, "output");
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = app_initializer->getInputDatabase();

    // The integrator borrows the boundary condition objects.
    std::vector<std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>>> u_bc_coefs_owner(NDIM);
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM>> rho_bc_coef;
    const bool conservative = input_db->keyExists("INSVCStaggeredConservativeHierarchyIntegrator");
    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
    if (conservative)
    {
        time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
            "INSVCStaggeredConservativeHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
    }
    else
    {
        time_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
            "INSVCStaggeredNonConservativeHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
    }
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> viscosity_db =
        app_initializer->getComponentDatabase("CarreauYasudaViscosity");
    auto viscosity = std::make_unique<ObservedViscosity>(
        std::make_unique<CarreauYasudaViscosity>("CarreauYasudaViscosity", viscosity_db),
        *time_integrator,
        viscosity_db->getDouble("infinite_shear_viscosity"),
        viscosity_db->getDouble("zero_shear_viscosity"));
    ObservedViscosity& observer = *viscosity;
    time_integrator->registerPreprocessIntegrateHierarchyCallback(&ObservedViscosity::start_step, &observer);
    time_integrator->setViscosityModel(std::move(viscosity));

    // The integrator that advances the solution, which is either the fluid integrator or its IB parent.
    SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> driver = time_integrator;
    SAMRAI::tbox::Pointer<IBMethod> ib_method_ops;
    if (input_db->keyExists("IBHierarchyIntegrator"))
    {
        ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        driver = new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                                   app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                   ib_method_ops,
                                                   time_integrator);
    }
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> grid_geometry =
        new SAMRAI::geom::CartesianGridGeometry<NDIM>("CartesianGeometry",
                                                      app_initializer->getComponentDatabase("CartesianGeometry"));
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy =
        new SAMRAI::hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitialize<NDIM>> error_detector =
        new SAMRAI::mesh::StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", driver, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    SAMRAI::tbox::Pointer<SAMRAI::mesh::BergerRigoutsos<NDIM>> box_generator =
        new SAMRAI::mesh::BergerRigoutsos<NDIM>();
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM>> load_balancer =
        new SAMRAI::mesh::LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_algorithm =
        new SAMRAI::mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                                  app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                  error_detector,
                                                  box_generator,
                                                  load_balancer);

    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_init = new IBTK::muParserCartGridFunction(
        "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
    time_integrator->registerVelocityInitialConditions(u_init);
    if (conservative)
    {
        time_integrator->registerMassDensityVariable(new SAMRAI::pdat::SideVariable<NDIM, double>("rho"));
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_init = new IBTK::muParserCartGridFunction(
            "rho_init", app_initializer->getComponentDatabase("MassDensityInitialConditions"), grid_geometry);
        time_integrator->registerMassDensityInitialConditions(rho_init);
    }

    // Directions without periodicity use the boundary conditions in the input file.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    if (grid_geometry->getPeriodicShift().min() == 0)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs_owner[d] = std::make_unique<IBTK::muParserRobinBcCoefs>(
                "u_bc_coefs_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry);
            u_bc_coefs[d] = u_bc_coefs_owner[d].get();
        }
        time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        if (conservative)
        {
            rho_bc_coef = std::make_unique<IBTK::muParserRobinBcCoefs>(
                "rho_bc_coef", app_initializer->getComponentDatabase("MassDensityBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef.get());
        }
    }

    if (ib_method_ops)
    {
        SAMRAI::tbox::Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        finest_ln = gridding_algorithm->getMaxLevels() - 1;
        ib_initializer->setStructureNamesOnLevel(finest_ln, { "marker" });
        ib_initializer->registerInitStructureFunction(generate_structure);
        ib_method_ops->registerLInitStrategy(ib_initializer);
    }

    driver->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
    app_initializer.setNull();

    while (!IBTK::rel_equal_eps(driver->getIntegratorTime(), driver->getEndTime()) && driver->stepsRemaining())
    {
        driver->advanceHierarchy(driver->getMaximumTimeStepSize());
    }
    observer.finish();

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int u_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                           time_integrator->getCurrentContext());
    const int mu_idx = var_db->mapVariableAndContextToIndex(time_integrator->getViscosityVariable(),
                                                            time_integrator->getCurrentContext());
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, 0, finest_ln);
    SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, 0, finest_ln);
    const int wgt_sc_idx = time_integrator->getHierarchyMathOps()->getSideWeightPatchDescriptorIndex();

    SAMRAI::tbox::plog << std::setprecision(10) << "time = " << driver->getIntegratorTime() << '\n'
                       << "u L2-norm = " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << '\n'
                       << "u max-norm = " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << '\n'
                       << "mu min = " << hier_cc_data_ops.min(mu_idx) << '\n'
                       << "mu max = " << hier_cc_data_ops.max(mu_idx) << '\n'
                       << "time steps = " << observer.num_steps << '\n'
                       << "cycles per time step = " << observer.num_cycles << '\n'
                       << "viscosity updates = " << observer.num_updates << '\n'
                       << "viscosity update errors = " << observer.num_errors << '\n';
    return 0;
}
