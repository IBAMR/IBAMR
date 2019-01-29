// Filename LSLocateStructureInterface.cpp
// Created by Nishant Nangia

#include "LSLocateStructureInterface.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// STATIC ///////////////////////////////////////
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
                                                       BargeInterface* barge)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_lag_data_manager(lag_data_manager),
      d_vol_elem(vol_elem),
      d_barge(barge)
{
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
    setLevelSetPatchDataByGeometry(D_idx, hier_math_ops, time, initial_time);
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSLocateStructureInterface::setLevelSetPatchDataByGeometry(int D_idx,
                                                           Pointer<HierarchyMathOps> hier_math_ops,
                                                           double /*time*/,
                                                           bool /*initial_time*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the initial condition for locating the interface
    double Xcom = d_barge->COM(0);
    double Ycom = d_barge->COM(1);
    double L = d_barge->length;
    double W = d_barge->width;

    // Note that these must correspond to one of the corners of the barge
    std::vector<IBTK::Vector> corners;
    getExtremeCoords(corners, hier_math_ops);

    // Store the coordinates
    // Note that I think this assumes that the corners aren't crossing the COM threshold
    for (std::vector<IBTK::Vector>::iterator it = corners.begin(); it != corners.end(); ++it)
    {
        IBTK::Vector c = *it;
        if (c(0) < Xcom)
        {
            if (c(1) > Ycom)
                d_barge->TL = c;
            else
                d_barge->BL = c;
        }
        else
        {
            if (c(1) > Ycom)
                d_barge->TR = c;
            else
                d_barge->BR = c;
        }
    }

    // Store the angle with the horizontal
    IBTK::Vector a = d_barge->BR - d_barge->BL;
    d_barge->theta = std::atan(a(1) / a(0));

    // Analytical distance away from the rectangle
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
                const double relx = X[0] - Xcom;
                const double rely = X[1] - Ycom;
                const double rotx = relx * std::cos(-d_barge->theta) - rely * std::sin(-d_barge->theta);
                const double roty = relx * std::sin(-d_barge->theta) + rely * std::cos(-d_barge->theta);

                bool inside = (-0.5 * L < rotx && rotx < 0.5 * L && -0.5 * W < roty && roty < 0.5 * W);
                if (!inside)
                {
                    const double dx = std::max(abs(rotx) - 0.5 * L, 0.0);
                    const double dy = std::max(abs(roty) - 0.5 * W, 0.0);
                    (*D_data)(ci) = std::sqrt(dx * dx + dy * dy);
                }
                else
                {
                    double dx_min = std::min(abs(rotx + 0.5 * L), abs(0.5 * L - rotx));
                    double dy_min = std::min(abs(roty + 0.5 * W), abs(0.5 * W - roty));
                    (*D_data)(ci) = -std::min(dx_min, dy_min);
                }

#endif

#if (NDIM == 3)
                TBOX_ERROR("Presently not implemented");
#endif
            }
        }
    }
    return;
} // setLevelSetPatchDataByGeometry

void
LSLocateStructureInterface::getExtremeCoords(std::vector<IBTK::Vector>& corners,
                                             Pointer<HierarchyMathOps> hier_math_ops)
{
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    double y_xmin = std::numeric_limits<double>::quiet_NaN();
    double y_xmax = std::numeric_limits<double>::quiet_NaN();
    double x_ymin = std::numeric_limits<double>::quiet_NaN();
    double x_ymax = std::numeric_limits<double>::quiet_NaN();

    IBTK::Vector c_xmin, c_xmax, c_ymin, c_ymax;

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

            // Find the corners of the barge
            if (X[0] < xmin)
            {
                xmin = X[0];
                y_xmin = X[1];
            }
            if (xmax < X[0])
            {
                xmax = X[0];
                y_xmax = X[1];
            }
            if (X[1] < ymin)
            {
                ymin = X[1];
                x_ymin = X[0];
            }
            if (ymax < X[1])
            {
                ymax = X[1];
                x_ymax = X[0];
            }
        }

        X_data[ln]->restoreArrays();
    } // all levels

    // For each of the coordinates, carry out reduction but keep track of which processor has the extreme value
    int rank_xmin, rank_xmax, rank_ymin, rank_ymax;
    xmin = SAMRAI_MPI::minReduction(xmin, &rank_xmin);
    xmax = SAMRAI_MPI::maxReduction(xmax, &rank_xmax);
    ymin = SAMRAI_MPI::minReduction(ymin, &rank_ymin);
    ymax = SAMRAI_MPI::maxReduction(ymax, &rank_ymax);

    // Broadcast via minReduction the missing coordinate from the appropriate rank.
    const int num_corners = 4;
    double other_coords[num_corners];
    other_coords[0] = y_xmin;
    other_coords[1] = y_xmax;
    other_coords[2] = x_ymin;
    other_coords[3] = x_ymax;
    if (SAMRAI_MPI::getRank() != rank_xmin) other_coords[0] = std::numeric_limits<double>::max();
    if (SAMRAI_MPI::getRank() != rank_xmax) other_coords[1] = std::numeric_limits<double>::max();
    if (SAMRAI_MPI::getRank() != rank_ymin) other_coords[2] = std::numeric_limits<double>::max();
    if (SAMRAI_MPI::getRank() != rank_ymax) other_coords[3] = std::numeric_limits<double>::max();
    SAMRAI_MPI::minReduction(&other_coords[0], num_corners);
    y_xmin = other_coords[0];
    y_xmax = other_coords[1];
    x_ymin = other_coords[2];
    x_ymax = other_coords[3];

    // Store the coordinates
    c_xmin(0) = xmin;
    c_xmin(1) = y_xmin;
    c_xmax(0) = xmax;
    c_xmax(1) = y_xmax;
    c_ymin(1) = ymin;
    c_ymin(0) = x_ymin;
    c_ymax(1) = ymax;
    c_ymax(0) = x_ymax;
    corners.resize(num_corners);
    corners[0] = c_xmin;
    corners[1] = c_xmax;
    corners[2] = c_ymin;
    corners[3] = c_ymax;

    return;
} // getExtremeCoords
