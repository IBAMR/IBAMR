// Filename LSLocateBargeInterface.h
//
// Copyright (c) 2002-2019, Amneet Bhalla and Nishant Nangia
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <CartesianGridGeometry.h>
#include <tbox/Pointer.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/ConstraintIBMethod.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_utilities.h>

#ifndef included_IBAMR_multiphase_flow_LSLocateBargeInterface
#define included_IBAMR_multiphase_flow_LSLocateBargeInterface

// Struct to maintain the properties of the barge interface
struct BargeInterface
{
    double length;
    double width;
    IBTK::Vector COM;
    IBTK::Vector TR, TL, BR, BL;
    double theta;
};

/*!
 * \brief class LSLocateBargeInterface is a utility class which is used to
 * identify the barge interface for level set computations.
 */
class LSLocateBargeInterface
{
public:
    LSLocateBargeInterface(const std::string& object_name,
                           SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                           IBTK::LDataManager* lag_data_manager,
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

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              double time,
                              bool initial_time)
    {
        setLevelSetPatchDataByGeometry(D_idx, hier_math_ops, time, initial_time);
        return;
    } // setLevelSetPatchData

    /////////////////////////////// PRIVATE //////////////////////////////////////
private:
    LSLocateBargeInterface();
    LSLocateBargeInterface& operator=(const LSLocateBargeInterface& that);
    LSLocateBargeInterface(const LSLocateBargeInterface& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to the advection-diffusion solver
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;

    /*!
     * IB information
     */
    IBTK::LDataManager* d_lag_data_manager;

    /*!
     * Volume element
     */
    double d_vol_elem;

    /*!
     * Barge structure location.
     */
    BargeInterface* d_barge;

    /*!
     * Reinitialize the level set information by geometry.
     */
    void setLevelSetPatchDataByGeometry(int D_idx,
                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        double /*time*/,
                                        bool /*initial_time*/)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
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
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
                for (SAMRAI::hier::Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    SAMRAI::pdat::CellIndex<NDIM> ci(it());

                    // Get physical coordinates
                    IBTK::Vector X = IBTK::Vector::Zero();
                    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();
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

    /*!
     * Get the extreme coordinate points of the barge.
     */
    void getExtremeCoords(std::vector<IBTK::Vector>& corners,
                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops)
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

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data(finest_ln + 1,
                                                                SAMRAI::tbox::Pointer<IBTK::LData>(NULL));
        X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;

            // Get pointer to LData
            boost::multi_array_ref<double, 2>& X_boost_data = *X_data[ln]->getLocalFormVecArray();
            const SAMRAI::tbox::Pointer<IBTK::LMesh> mesh = d_lag_data_manager->getLMesh(ln);
            const std::vector<IBTK::LNode*>& local_nodes = mesh->getLocalNodes();

            for (std::vector<IBTK::LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
            {
                const IBTK::LNode* const node_idx = *cit;
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
        xmin = SAMRAI::tbox::SAMRAI_MPI::minReduction(xmin, &rank_xmin);
        xmax = SAMRAI::tbox::SAMRAI_MPI::maxReduction(xmax, &rank_xmax);
        ymin = SAMRAI::tbox::SAMRAI_MPI::minReduction(ymin, &rank_ymin);
        ymax = SAMRAI::tbox::SAMRAI_MPI::maxReduction(ymax, &rank_ymax);

        // Broadcast via minReduction the missing coordinate from the appropriate rank.
        const int num_corners = 4;
        double other_coords[num_corners];
        other_coords[0] = y_xmin;
        other_coords[1] = y_xmax;
        other_coords[2] = x_ymin;
        other_coords[3] = x_ymax;
        if (SAMRAI::tbox::SAMRAI_MPI::getRank() != rank_xmin) other_coords[0] = std::numeric_limits<double>::max();
        if (SAMRAI::tbox::SAMRAI_MPI::getRank() != rank_xmax) other_coords[1] = std::numeric_limits<double>::max();
        if (SAMRAI::tbox::SAMRAI_MPI::getRank() != rank_ymin) other_coords[2] = std::numeric_limits<double>::max();
        if (SAMRAI::tbox::SAMRAI_MPI::getRank() != rank_ymax) other_coords[3] = std::numeric_limits<double>::max();
        SAMRAI::tbox::SAMRAI_MPI::minReduction(&other_coords[0], num_corners);
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
};

inline void
callLSLocateBargeInterfaceCallbackFunction(int D_idx,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           double time,
                                           bool initial_time,
                                           void* ctx)
{
    // Set the level set information
    static LSLocateBargeInterface* ptr_LSLocateBargeInterface = static_cast<LSLocateBargeInterface*>(ctx);
    ptr_LSLocateBargeInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateBargeInterfaceCallbackFunction

#endif
