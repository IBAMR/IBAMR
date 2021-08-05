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

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>

#include "LSLocateStructureInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

LSLocateStructureInterface::LocateStructureMethod LSLocateStructureInterface::s_locate_method =
    LSLocateStructureInterface::GEOMETRY_METHOD;

void
callLSLocateStructureInterfaceCallbackFunction(int D_idx,
                                               Pointer<HierarchyMathOps> hier_math_ops,
                                               double time,
                                               bool initial_time,
                                               void* ctx)
{
    // Set the level set information
    static LSLocateStructureInterface* ptr_LSLocateStructureInterface = static_cast<LSLocateStructureInterface*>(ctx);
    ptr_LSLocateStructureInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateStructureInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateStructureInterface::LSLocateStructureInterface(const std::string& object_name,
                                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                       Pointer<CellVariable<NDIM, double> > ls_var,
                                                       LDataManager* lag_data_manager,
                                                       double vol_elem,
                                                       WedgeInterface* wedge)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_lag_data_manager(lag_data_manager),
      d_vol_elem(vol_elem),
      d_wedge(wedge)
{
    if (d_wedge->wedge_locate_method == "GEOMETRY_METHOD")
        s_locate_method = GEOMETRY_METHOD;
    else if (d_wedge->wedge_locate_method == "IB_SPREADING_METHOD")
        s_locate_method = IB_SPREADING_METHOD;
    else
        TBOX_ERROR("Unknown method to locate wedge specified.\n");

    // intentionally left blank
    return;
} // LSLocateStructureInterface

LSLocateStructureInterface::~LSLocateStructureInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateStructureInterface::setLevelSetPatchData(int D_idx,
                                                 Pointer<HierarchyMathOps> hier_math_ops,
                                                 double time,
                                                 bool initial_time)
{
    if (s_locate_method == GEOMETRY_METHOD)
    {
        setLevelSetPatchDataByGeometry(D_idx, hier_math_ops, time, initial_time);
    }
    else if (s_locate_method == IB_SPREADING_METHOD)
    {
        setLevelSetPatchDataBySpreading(D_idx, hier_math_ops, time, initial_time);
    }
    else
    {
        TBOX_ERROR(
            "LSLocateStructureInterface::setLevelSetPatchData(): Unknown method for locating structure interface. \n");
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSLocateStructureInterface::setLevelSetPatchDataBySpreading(int D_idx,
                                                            Pointer<HierarchyMathOps> hier_math_ops,
                                                            double /*time*/,
                                                            bool /*initial_time*/)
{
    // In this version of this class, the initial level set location is set by spreading
    // from the Lagrangian structure

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    // Spread
    std::vector<Pointer<LData> > X_data(finest_ln + 1, Pointer<LData>(NULL));
    X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    std::vector<Pointer<LData> > lag_phi(finest_ln + 1, Pointer<LData>(NULL));
    lag_phi[finest_ln] = d_lag_data_manager->getLData("PHI", finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        // Set D_idx to zero
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
            D_data->fill(0.0);
        }

        // Set lag data to -1.0
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;

        Vec petsc_vec = lag_phi[ln]->getVec();
        VecSet(petsc_vec, 1.0 * d_vol_elem);
    }

    std::vector<Pointer<LData> > PHI_data(finest_ln + 1, Pointer<LData>(NULL));
    PHI_data[finest_ln] = lag_phi[finest_ln];
    d_lag_data_manager->spread(D_idx, PHI_data, X_data, (RobinPhysBdryPatchStrategy*)NULL);
    return;
} // setLevelSetPatchDataBySpreading

void
LSLocateStructureInterface::setLevelSetPatchDataByGeometry(int D_idx,
                                                           Pointer<HierarchyMathOps> hier_math_ops,
                                                           double /*time*/,
                                                           bool /*initial_time*/)
{
    // In this version of this class, the initial level set location is set to be
    // exact since we always know the wedge geometry

    static const double m = tan(d_wedge->wedge_angle);
    static const double wedge_height = 0.5 * (d_wedge->wedge_length) * m;

#if (NDIM == 3)
    // Normal of right plane
    static const double ar = -sin(d_wedge->wedge_angle);
    static const double br = 0.0;
    static const double cr = cos(d_wedge->wedge_angle);

    // Normal of left plane
    static const double al = sin(d_wedge->wedge_angle);
    static const double bl = 0.0;
    static const double cl = cos(d_wedge->wedge_angle);
#endif

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the initial condition for locating the interface
    IBTK::Vector& X0 = d_wedge->X0;

    // Get the coordinate of the bottom point of the wedge
    X0(NDIM - 1) = getMinimumWedgeCoord(hier_math_ops);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Get physical coordinates
                IBTK::Vector X = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const SAMRAI::hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

#if (NDIM == 2)

                double distance[3];

                // a) Top plane
                const double H = wedge_height + X0(1);
                distance[0] = X(1) - H;

                // b) Right line of the form y = slope*x + c
                const double cr = X0(1) - m * X0(0);
                distance[1] = (m * X(0) - X(1) + cr) / sqrt(1.0 + m * m);

                // c) Left line of the form y = slope*x + c
                const double cl = X0(1) + m * X0(0);
                distance[2] = (-m * X(0) - X(1) + cl) / sqrt(1.0 + m * m);

                (*D_data)(ci) = *std::max_element(distance, distance + 3);

#endif

#if (NDIM == 3)

                double distance[5];

                // a) Top plane
                const double H = wedge_height + X0(2);
                distance[0] = X(2) - H;

                // b) Right plane of the form: ax + by + cz + d = 0
                const double dr = 0.0 - ar * X0(0) - br * X0(1) - cr * X0(2);
                distance[1] = -dr - ar * X(0) - br * X(1) - cr * X(2);

                // c) Left plane  of the form ax + by + cz + d = 0
                const double dl = 0.0 - al * X0(0) - bl * X0(1) - cl * X0(2);
                distance[2] = -dl - al * X(0) - bl * X(1) - cl * X(2);

                // d) Front plane
                const double Yf = X0(1) - (d_wedge->wedge_length) / 2.0;
                distance[3] = Yf - X(1);

                // d) Back plane
                const double Yb = X0(1) + (d_wedge->wedge_length) / 2.0;
                distance[4] = X(1) - Yb;

                (*D_data)(ci) = *std::max_element(distance, distance + 5);
#endif
            }
        }
    }
    return;
} // setLevelSetPatchDataByGeometry

double
LSLocateStructureInterface::getMinimumWedgeCoord(Pointer<HierarchyMathOps> hier_math_ops)
{
    double wedge_min = 1.0e12;

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> > X_data(finest_ln + 1, Pointer<LData>(NULL));
    X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get pointer to LData
        boost::multi_array_ref<double, 2>& X_boost_data = *X_data[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_lag_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            double* const X = &X_boost_data[local_idx][0];

            // Find the minimum coordinate point of the wedge
            wedge_min = std::min(X[NDIM - 1], wedge_min);
        }

        X_data[ln]->restoreArrays();
    } // all levels
    wedge_min = IBTK_MPI::minReduction(wedge_min);

    return wedge_min;
} // getMinimumWedgeCoord
