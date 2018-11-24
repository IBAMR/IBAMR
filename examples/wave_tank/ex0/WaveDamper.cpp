// Filename: WaveDamper.cpp
// Created on Sep 04, 2018 by Amneet Bhalla and Nishant Nangia

// APPLICATION INCLUDES
#include "WaveDamper.h"
#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callWaveDamperCallbackFunction(double current_time,
                               double new_time,
                               bool skip_synchronize_new_state_data,
                               int num_cycles,
                               void* ctx)
{
    static WaveDamper* ptr_WaveDamper = static_cast<WaveDamper*>(ctx);
    ptr_WaveDamper->DampOutletWaves(current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;

} // callWaveDamperCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

WaveDamper::WaveDamper(const std::string& object_name,
                       Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator,
                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                       Pointer<Variable<NDIM> > phi_var,
                       Pointer<Database> input_db)
    : d_object_name(object_name),
      d_ins_hier_integrator(ins_hier_integrator),
      d_adv_diff_hier_integrator(adv_diff_hier_integrator),
      d_phi_var(phi_var)
{
    d_x_zone_start = input_db->getDouble("x_zone_start");
    d_x_zone_end = input_db->getDouble("x_zone_end");
    d_depth = input_db->getDouble("depth");
    d_alpha = input_db->getDoubleWithDefault("alpha", 2.0);
    return;
} // WaveDamper

WaveDamper::~WaveDamper()
{
    // intentionally left blank
    return;

} //~WaveDamper

void
WaveDamper::DampOutletWaves(double /*current_time*/,
                            double /*new_time*/,
                            bool /*skip_synchronize_new_state_data*/,
                            int /*num_cycles*/)
{
    pout << "Damping waves at outlet..." << std::endl;
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_ins_hier_integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                         d_ins_hier_integrator->getNewContext());

    Pointer<CellVariable<NDIM, double> > phi_cc_var = d_phi_var;
    if (!phi_cc_var) TBOX_ERROR("WaveDamper::DampOutletWaves: Level set variable must be cell centered!");
    int phi_new_idx = var_db->mapVariableAndContextToIndex(phi_cc_var, d_adv_diff_hier_integrator->getNewContext());
    static const int dir = NDIM == 2 ? 1 : 2;

    // Compute alpha
    pout << "Damping factor alpha = " << d_alpha << std::endl;

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_new_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, SideIndex<NDIM>::Lower);
                    const double x_posn =
                        patch_x_lower[0] +
                        patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + axis == 0 ? 0.0 : 0.5);
                    if (x_posn >= d_x_zone_start)
                    {
                        double xtilde = (x_posn - d_x_zone_start) / (d_x_zone_end - d_x_zone_start);
                        double gamma = 1.0 - (exp(std::pow(xtilde, d_alpha)) - 1.0) / (exp(1.0) - 1.0);
                        (*u_data)(i_side, 0) *= gamma;
                    }
                }
            }
        }
    }

    // Modify the level set in the absorption zone.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_new_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                Index<NDIM> i = it();

                const double x_posn =
                    patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const double z_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);

                if (x_posn >= d_x_zone_start && x_posn <= d_x_zone_end)
                {
                    double xtilde = (x_posn - d_x_zone_start) / (d_x_zone_end - d_x_zone_start);
                    double gamma = 1.0 - (exp(std::pow(xtilde, d_alpha)) - 1.0) / (exp(1.0) - 1.0);
                    double analytical = z_posn - d_depth;

                    (*phi_data)(i, 0) = (1.0 - gamma) * analytical + gamma * (*phi_data)(i, 0);
                }
            }
        }
    }

    return;
} // WaveDamperCoefPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
