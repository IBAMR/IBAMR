// Filename: LMarkerUtilities.cpp
// Created on 28 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <string>
#include <vector>

#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LMarker.h"
#include "ibtk/LMarkerCoarsen.h"
#include "ibtk/LMarkerRefine.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerUtilities.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSetDataIterator.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace xfer
{

class CoarsenPatchStrategy;

class RefinePatchStrategy;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline std::string discard_comments(const std::string& input_string)
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
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

unsigned int LMarkerUtilities::readMarkerPositions(std::vector<Point>& mark_init_posns,
                                                   const std::string& mark_input_file_name,
                                                   const boost::shared_ptr<CartesianGridGeometry>& grid_geom)
{
    if (mark_input_file_name.empty()) return 0;

    // Read in the initial marker positions.
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int mpi_rank = comm.getRank();
    const int mpi_size = comm.getSize();

    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();

    int num_mark = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        if (rank == mpi_rank)
        {
            std::string line_string;
            std::ifstream file_stream(mark_input_file_name.c_str(), std::ios::in);

            // The first entry in the file is the number of markers.
            if (!std::getline(file_stream, line_string))
            {
                TBOX_ERROR(
                    "LMarkerUtilities::readMarkerPositions():\n  Premature end to input file "
                    "encountered before line 1 of file "
                    << mark_input_file_name << "\n");
            }
            else
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
                if (!(line_stream >> num_mark))
                {
                    TBOX_ERROR(
                        "LMarkerUtilities::readMarkerPositions():\n  Invalid entry in input "
                        "file "
                        "encountered on line 1 of file "
                        << mark_input_file_name << "\n");
                }
            }

            if (num_mark <= 0)
            {
                TBOX_ERROR(
                    "LMarkerUtilities::readMarkerPositions():\n  Invalid entry in input file "
                    "encountered on line 1 of file "
                    << mark_input_file_name << "\n");
            }

            // Each successive line provides the initial position of each
            // marker in the input file.
            mark_init_posns.resize(num_mark);
            for (int k = 0; k < num_mark; ++k)
            {
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(
                        "LMarkerUtilities::readMarkerPositions():\n  Premature end to input "
                        "file "
                        "encountered before line "
                        << k + 2 << " of file " << mark_input_file_name << "\n");
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (!(line_stream >> mark_init_posns[k][d]))
                        {
                            TBOX_ERROR(
                                "LMarkerUtilities::readMarkerPositions():\n  Invalid entry in "
                                "input file encountered on line "
                                << k + 2 << " of file " << mark_input_file_name << "\n");
                        }
                    }

                    // Ensure the initial marker position lies within the
                    // physical domain.
                    const double* const X = mark_init_posns[k].data();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (MathUtilities<double>::equalEps(X[d], grid_xLower[d]))
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker intersecting lower physical "
                                          "boundary.\n"
                                       << "  please ensure that all markers are within the "
                                          "computational "
                                          "domain." << std::endl);
                        }
                        else if (X[d] <= grid_xLower[d])
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker below lower physical boundary\n"
                                       << "  please ensure that all markers are within the "
                                          "computational domain." << std::endl);
                        }

                        if (MathUtilities<double>::equalEps(X[d], grid_xUpper[d]))
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker intersecting upper physical "
                                          "boundary.\n"
                                       << "  please ensure that all markers are within the "
                                          "computational "
                                          "domain." << std::endl);
                        }
                        else if (X[d] >= grid_xUpper[d])
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker above upper physical boundary\n"
                                       << "  please ensure that all markers are within the "
                                          "computational domain." << std::endl);
                        }
                    }
                }
            }
        }
    }
    return num_mark;
}

void LMarkerUtilities::eulerStep(const int mark_current_idx,
                                 const int mark_new_idx,
                                 const int u_current_idx,
                                 const double dt,
                                 const std::string& weighting_fcn,
                                 const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                 const int coarsest_ln_in,
                                 const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            auto u_current_data = patch->getPatchData(u_current_idx);
            auto u_cc_current_data = boost::dynamic_pointer_cast<CellData<double> >(u_current_data);
            auto u_sc_current_data = boost::dynamic_pointer_cast<SideData<double> >(u_current_data);
            const bool is_cc_data = u_cc_current_data != NULL;
            const bool is_sc_data = u_sc_current_data != NULL;
            auto mark_current_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_current_idx));
            auto mark_new_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_new_idx));

            const unsigned int num_patch_marks = countMarkersOnPatch(mark_current_data);
            TBOX_ASSERT(num_patch_marks == countMarkersOnPatch(mark_new_data));

            // Collect the local marker positions at time n.
            std::vector<double> X_mark_current;
            collectMarkerPositionsOnPatch(X_mark_current, mark_current_data);

            // Compute U_mark(n) = u(X_mark(n),n).
            std::vector<double> U_mark_current(X_mark_current.size());
            if (is_cc_data)
                LEInteractor::interpolate(U_mark_current, NDIM, X_mark_current, NDIM, u_cc_current_data, patch,
                                          patch_box, weighting_fcn);
            if (is_sc_data)
                LEInteractor::interpolate(U_mark_current, NDIM, X_mark_current, NDIM, u_sc_current_data, patch,
                                          patch_box, weighting_fcn);

            // Compute X_mark(n+1) = X_mark(n) + dt*U_mark(n).
            std::vector<double> X_mark_new(X_mark_current.size());
            for (unsigned int k = 0; k < NDIM * num_patch_marks; ++k)
            {
                X_mark_new[k] = X_mark_current[k] + dt * U_mark_current[k];
            }

            // Prevent markers from leaving the computational domain through
            // physical boundaries (but *not* through periodic boundaries).
            auto grid_geom = BOOST_CAST<CartesianGridGeometry>(hierarchy->getGridGeometry());
            preventMarkerEscape(X_mark_new, grid_geom);

            // Store the local marker velocities at at time n, and the marker
            // positions at time n+1.
            resetMarkerVelocitiesOnPatch(U_mark_current, mark_current_data);
            resetMarkerPositionsOnPatch(X_mark_new, mark_new_data);
        }
    }
    return;
}

void LMarkerUtilities::midpointStep(const int mark_current_idx,
                                    const int mark_new_idx,
                                    const int u_half_idx,
                                    const double dt,
                                    const std::string& weighting_fcn,
                                    const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                    const int coarsest_ln_in,
                                    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            auto u_half_data = patch->getPatchData(u_half_idx);
            auto u_cc_half_data = boost::dynamic_pointer_cast<CellData<double> >(u_half_data);
            auto u_sc_half_data = boost::dynamic_pointer_cast<SideData<double> >(u_half_data);
            const bool is_cc_data = u_cc_half_data != NULL;
            const bool is_sc_data = u_sc_half_data != NULL;
            auto mark_current_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_current_idx));
            auto mark_new_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_new_idx));

            const unsigned int num_patch_marks = countMarkersOnPatch(mark_current_data);
            TBOX_ASSERT(num_patch_marks == countMarkersOnPatch(mark_new_data));

            // Collect the local marker positions at time n and predicted marker
            // positions at time n+1.
            std::vector<double> X_mark_current;
            collectMarkerPositionsOnPatch(X_mark_current, mark_current_data);
            std::vector<double> X_mark_new;
            collectMarkerPositionsOnPatch(X_mark_new, mark_new_data);

            // Set X(n+1/2) = 0.5*(X(n)+X(n+1)).
            std::vector<double> X_mark_half(NDIM * num_patch_marks);
            for (unsigned int k = 0; k < NDIM * num_patch_marks; ++k)
            {
                X_mark_half[k] = 0.5 * (X_mark_current[k] + X_mark_new[k]);
            }

            // Compute U_mark(n+1/) = u(X_mark(n+1/2),n+1/2).
            std::vector<double> U_mark_half(X_mark_half.size());
            if (is_cc_data)
                LEInteractor::interpolate(U_mark_half, NDIM, X_mark_half, NDIM, u_cc_half_data, patch, patch_box,
                                          weighting_fcn);
            if (is_sc_data)
                LEInteractor::interpolate(U_mark_half, NDIM, X_mark_half, NDIM, u_sc_half_data, patch, patch_box,
                                          weighting_fcn);

            // Compute X_mark(n+1) = X_mark(n) + dt*U_mark(n+1/2).
            for (unsigned int k = 0; k < NDIM * num_patch_marks; ++k)
            {
                X_mark_new[k] = X_mark_current[k] + dt * U_mark_half[k];
            }

            // Prevent markers from leaving the computational domain through
            // physical boundaries (but *not* through periodic boundaries).
            auto grid_geom = BOOST_CAST<CartesianGridGeometry>(hierarchy->getGridGeometry());
            preventMarkerEscape(X_mark_new, grid_geom);

            // Store the local marker positions at time n+1.
            resetMarkerPositionsOnPatch(X_mark_new, mark_new_data);
        }
    }
    return;
}

void LMarkerUtilities::trapezoidalStep(const int mark_current_idx,
                                       const int mark_new_idx,
                                       const int u_new_idx,
                                       const double dt,
                                       const std::string& weighting_fcn,
                                       const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                       const int coarsest_ln_in,
                                       const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            auto u_new_data = patch->getPatchData(u_new_idx);
            auto u_cc_new_data = boost::dynamic_pointer_cast<CellData<double> >(u_new_data);
            auto u_sc_new_data = boost::dynamic_pointer_cast<SideData<double> >(u_new_data);
            const bool is_cc_data = u_cc_new_data != NULL;
            const bool is_sc_data = u_sc_new_data != NULL;
            auto mark_current_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_current_idx));
            auto mark_new_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_new_idx));

            const unsigned int num_patch_marks = countMarkersOnPatch(mark_current_data);
            TBOX_ASSERT(num_patch_marks == countMarkersOnPatch(mark_new_data));

            // Collect the local marker positions at time n and predicted marker
            // positions at time n+1.
            std::vector<double> X_mark_current;
            collectMarkerPositionsOnPatch(X_mark_current, mark_current_data);
            std::vector<double> X_mark_new;
            collectMarkerPositionsOnPatch(X_mark_new, mark_new_data);

            // Collect the local marker velocities at time n.
            std::vector<double> U_mark_current;
            collectMarkerVelocitiesOnPatch(U_mark_current, mark_current_data);

            // Compute U_mark(n+1/) = u(X_mark(n+1/2),n+1/2).
            std::vector<double> U_mark_new(X_mark_new.size());
            if (is_cc_data)
                LEInteractor::interpolate(U_mark_new, NDIM, X_mark_new, NDIM, u_cc_new_data, patch, patch_box,
                                          weighting_fcn);
            if (is_sc_data)
                LEInteractor::interpolate(U_mark_new, NDIM, X_mark_new, NDIM, u_sc_new_data, patch, patch_box,
                                          weighting_fcn);

            // Set U(n+1/2) = 0.5*(U(n)+U(n+1)).
            std::vector<double> U_mark_half(NDIM * num_patch_marks);
            for (unsigned int k = 0; k < NDIM * num_patch_marks; ++k)
            {
                U_mark_half[k] = 0.5 * (U_mark_current[k] + U_mark_new[k]);
            }

            // Compute X_mark(n+1) = X_mark(n) + dt*U_mark(n+1/2).
            for (unsigned int k = 0; k < NDIM * num_patch_marks; ++k)
            {
                X_mark_new[k] = X_mark_current[k] + dt * U_mark_half[k];
            }

            // Prevent markers from leaving the computational domain through
            // physical boundaries (but *not* through periodic boundaries).
            auto grid_geom = BOOST_CAST<CartesianGridGeometry>(hierarchy->getGridGeometry());
            preventMarkerEscape(X_mark_new, grid_geom);

            // Store the local marker velocities at at time n, and the marker
            // positions at time n+1.
            resetMarkerVelocitiesOnPatch(U_mark_new, mark_new_data);
            resetMarkerPositionsOnPatch(X_mark_new, mark_new_data);
        }
    }
    return;
}

void LMarkerUtilities::collectMarkersOnPatchHierarchy(const int mark_idx, const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    const unsigned int num_marks_before_coarsening = countMarkers(mark_idx, hierarchy);
    const unsigned int num_marks_before_coarsening_level_0 = countMarkers(mark_idx, hierarchy, 0, 0);

    // Collect all marker data on the patch hierarchy.
    auto var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> var;
    var_db->mapIndexToVariable(mark_idx, var);
    int mark_scratch_idx = var_db->registerClonedPatchDataIndex(var, mark_idx);
    auto mark_coarsen_alg = boost::make_shared<CoarsenAlgorithm>(DIM);
    auto mark_coarsen_op = boost::make_shared<LMarkerCoarsen>();
    mark_coarsen_alg->registerCoarsen(mark_scratch_idx, mark_idx, mark_coarsen_op);
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        auto coarser_level = hierarchy->getPatchLevel(ln - 1);

        // Allocate scratch data.
        coarser_level->allocatePatchData(mark_scratch_idx);

        // Coarsen fine data onto coarser level.
        CoarsenPatchStrategy* mark_coarsen_op = NULL;
        mark_coarsen_alg->createSchedule(coarser_level, level, mark_coarsen_op)->coarsenData();

        // Merge the coarsened fine data with the coarse data.
        for (auto p = coarser_level->begin(), e = coarser_level->end(); p != e; ++p)
        {
            auto patch = *p;
            auto mark_current_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
            auto mark_scratch_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_scratch_idx));
            for (auto it = LMarkerSetData::SetIterator(*mark_scratch_data, /*begin*/ true), e = LMarkerSetData::SetIterator(*mark_scratch_data, /*begin*/ false); it != e; ++it)
            {
                const Index& i = it.getIndex();
                if (!mark_current_data->isElement(i))
                {
                    mark_current_data->appendItemPointer(i, new LMarkerSet());
                }
                LMarkerSet& dst_mark_set = *(mark_current_data->getItem(i));
                const LMarkerSet& src_mark_set = *it;
                dst_mark_set.insert(dst_mark_set.end(), src_mark_set.begin(), src_mark_set.end());
            }
        }

        // Clear the fine data.
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            auto mark_current_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
            mark_current_data->removeAllItems();
        }

        // Deallocate scratch data.
        coarser_level->deallocatePatchData(mark_scratch_idx);
    }
    var_db->removePatchDataIndex(mark_scratch_idx);
    mark_scratch_idx = -1;

    // Ensure that the total number of markers is correct.
    const unsigned int num_marks_after_coarsening = countMarkers(mark_idx, hierarchy);
    const unsigned int num_marks_after_coarsening_level_0 = countMarkers(mark_idx, hierarchy, 0, 0);
    if (num_marks_before_coarsening != num_marks_after_coarsening ||
        num_marks_before_coarsening != num_marks_after_coarsening_level_0)
    {
        TBOX_ERROR("LMarkerUtilities::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during collection to coarsest level\n"
                   << "  number of markers in hierarchy before collection to coarsest level = "
                   << num_marks_before_coarsening << "\n"
                   << "  number of markers on level 0   before collection to coarsest level = "
                   << num_marks_before_coarsening_level_0 << "\n"
                   << "  number of markers in hierarchy after  collection to coarsest level = "
                   << num_marks_after_coarsening << "\n"
                   << "  number of markers on level 0   after  collection to coarsest level = "
                   << num_marks_after_coarsening_level_0 << "\n");
    }

    // Reset the assignment of markers to Cartesian grid cells on the coarsest
    // level of the patch hierarchy.
    //
    // NOTE: It is important to do this only *after* collecting markers on the
    // patch hierarchy.  Otherwise, markers that have left a fine level through
    // the coarse-fine interface would be discarded by this procedure.
    auto mark_level_fill_alg = boost::make_shared<RefineAlgorithm>();
    const boost::shared_ptr<RefineOperator> no_refine_op;
    mark_level_fill_alg->registerRefine(mark_idx, mark_idx, mark_idx, no_refine_op);
    auto level = hierarchy->getPatchLevel(coarsest_ln);
    mark_level_fill_alg->createSchedule(level, NULL)->fillData(0.0);
    for (auto p = level->begin(); p != level->end(); ++p)
    {
        auto patch = *p;
        const Box& patch_box = patch->getBox();

        auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
        const Index& patch_lower = patch_box.lower();
        const Index& patch_upper = patch_box.upper();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const double* const patch_dx = pgeom->getDx();

        auto mark_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
        auto mark_data_new = boost::make_shared<LMarkerSetData>(mark_data->getBox(), mark_data->getGhostCellWidth());
        const Box& it_box = mark_data->getGhostBox();
        for (auto it = mark_data->begin(it_box), e = mark_data->end(it_box); it != e; ++it)
        {
            const LMarkerSet::value_type& mark = *it;
            const Point& X = mark->getPosition();
            const IntVector& offset = mark->getPeriodicOffset();
            Point X_shifted;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_shifted[d] = X[d] + static_cast<double>(offset(d)) * patch_dx[d];
            }
            const bool patch_owns_mark_at_new_loc =
                ((patch_x_lower[0] <= X_shifted[0]) && (X_shifted[0] < patch_x_upper[0]))
#if (NDIM > 1)
                && ((patch_x_lower[1] <= X_shifted[1]) && (X_shifted[1] < patch_x_upper[1]))
#if (NDIM > 2)
                && ((patch_x_lower[2] <= X_shifted[2]) && (X_shifted[2] < patch_x_upper[2]))
#endif
#endif
                ;
            if (patch_owns_mark_at_new_loc)
            {
                const Index i = IndexUtilities::getCellIndex(X_shifted, patch_x_lower, patch_x_upper, patch_dx,
                                                             patch_lower, patch_upper);
                if (!mark_data_new->isElement(i))
                {
                    mark_data_new->appendItemPointer(i, new LMarkerSet());
                }
                LMarkerSet& new_mark_set = *(mark_data_new->getItem(i));
                new_mark_set.push_back(mark);
            }
        }

        // Swap the old and new patch data pointers.
        patch->setPatchData(mark_idx, mark_data_new);
    }

    // Ensure that the total number of markers is correct.
    const unsigned int num_marks_after_posn_reset = countMarkers(mark_idx, hierarchy);
    const unsigned int num_marks_after_posn_reset_level_0 = countMarkers(mark_idx, hierarchy, 0, 0);
    if (num_marks_before_coarsening != num_marks_after_posn_reset ||
        num_marks_before_coarsening != num_marks_after_posn_reset_level_0)
    {
        TBOX_ERROR("LMarkerUtilities::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during position reset on coarsest level\n"
                   << "  number of markers in hierarchy before position reset on coarsest level = "
                   << num_marks_before_coarsening << "\n"
                   << "  number of markers in hierarchy after  position reset on coarsest level = "
                   << num_marks_after_posn_reset << "\n"
                   << "  number of markers on level 0   after  position reset on coarsest level = "
                   << num_marks_after_posn_reset_level_0 << "\n");
    }
    return;
}

void LMarkerUtilities::initializeMarkersOnLevel(const int mark_idx,
                                                const std::vector<Point>& mark_init_posns,
                                                const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                const int level_number,
                                                const bool initial_time,
                                                const boost::shared_ptr<PatchLevel>& old_level)
{
    auto level = hierarchy->getPatchLevel(level_number);

    // On the coarsest level of the patch hierarchy, copy marker data from the
    // old coarse level.  Otherwise, refine marker data from the coarsest level
    // of the patch hierarchy.
    if (initial_time && level_number == 0)
    {
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const Index& patch_lower = patch_box.lower();
            const Index& patch_upper = patch_box.upper();
            const double* const patch_x_lower = pgeom->getXLower();
            const double* const patch_x_upper = pgeom->getXUpper();
            const double* const patch_dx = pgeom->getDx();

            auto mark_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
            for (unsigned int k = 0; k < mark_init_posns.size(); ++k)
            {
                const Point& X = mark_init_posns[k];
                static const Vector U(Vector::Zero());
                const bool patch_owns_mark_at_loc = ((patch_x_lower[0] <= X[0]) && (X[0] < patch_x_upper[0]))
#if (NDIM > 1)
                                                    && ((patch_x_lower[1] <= X[1]) && (X[1] < patch_x_upper[1]))
#if (NDIM > 2)
                                                    && ((patch_x_lower[2] <= X[2]) && (X[2] < patch_x_upper[2]))
#endif
#endif
                    ;
                if (patch_owns_mark_at_loc)
                {
                    const Index i = IndexUtilities::getCellIndex(X, patch_x_lower, patch_x_upper, patch_dx, patch_lower,
                                                                 patch_upper);
                    if (!mark_data->isElement(i))
                    {
                        mark_data->appendItemPointer(i, new LMarkerSet());
                    }
                    LMarkerSet& new_mark_set = *(mark_data->getItem(i));
                    new_mark_set.push_back(boost::make_shared<LMarker>(k, X, U));
                }
            }
        }
    }
    else
    {
        if (old_level && level_number == 0)
        {
            auto copy_mark_alg = boost::make_shared<RefineAlgorithm>();
            boost::shared_ptr<RefineOperator> no_fill_op;
            copy_mark_alg->registerRefine(mark_idx, mark_idx, mark_idx, no_fill_op);
            auto dst_level = level;
            auto src_level = old_level;
            RefinePatchStrategy* refine_mark_op = NULL;
            copy_mark_alg->createSchedule(dst_level, src_level, refine_mark_op)->fillData(0.0);
        }
        else if (level_number > 0)
        {
            auto refine_mark_alg = boost::make_shared<RefineAlgorithm>();
            auto fill_op = boost::make_shared<LMarkerRefine>();
            refine_mark_alg->registerRefine(mark_idx, mark_idx, mark_idx, fill_op);
            RefinePatchStrategy* refine_mark_op = NULL;
            for (int ln = 1; ln <= level_number; ++ln)
            {
                auto dst_level = hierarchy->getPatchLevel(ln);
                refine_mark_alg->createSchedule(dst_level, NULL, ln - 1, hierarchy, refine_mark_op)->fillData(0.0);
            }
        }
    }
    return;
}

void LMarkerUtilities::pruneInvalidMarkers(const int mark_idx,
                                           const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                           const int coarsest_ln_in,
                                           const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    const int finest_hier_level_number = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= std::min(finest_ln, finest_hier_level_number - 1); ++ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        auto finer_level = hierarchy->getPatchLevel(ln + 1);
        BoxContainer refined_region_boxes = finer_level->getGlobalizedBoxLevel().getBoxes();
        const IntVector& ratio = finer_level->getRatioToCoarserLevel();
        refined_region_boxes.coarsen(ratio);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            auto mark_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
            const Box& ghost_box = mark_data->getGhostBox();
            for (auto b = refined_region_boxes.begin(), e = refined_region_boxes.end(); b != e; ++b)
            {
                const Box& refined_box = *b;
                const Box intersection = ghost_box * refined_box;
                if (!intersection.empty())
                {
                    mark_data->removeInsideBox(intersection);
                }
            }
        }
    }
    return;
}

unsigned int LMarkerUtilities::countMarkers(const int mark_idx,
                                            const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                            const int coarsest_ln_in,
                                            const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    int num_marks = 0;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            auto mark_data = BOOST_CAST<LMarkerSetData>(patch->getPatchData(mark_idx));
            num_marks += countMarkersOnPatch(mark_data);
        }
    }
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.AllReduce(&num_marks, 1, MPI_SUM);
    return num_marks;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

unsigned int LMarkerUtilities::countMarkersOnPatch(const boost::shared_ptr<LMarkerSetData>& mark_data)
{
    const Box patch_box = mark_data->getBox();
    unsigned int num_marks = 0;
    for (LMarkerSetData::SetIterator it(*mark_data, /*begin*/ true), e(*mark_data, /*begin*/ false); it != e; ++it)
    {
        if (patch_box.contains(it.getIndex()))
        {
            num_marks += it->size();
        }
    }
    return num_marks;
}

void LMarkerUtilities::collectMarkerPositionsOnPatch(std::vector<double>& X_mark,
                                                     const boost::shared_ptr<LMarkerSetData>& mark_data)
{
    X_mark.resize(NDIM * countMarkersOnPatch(mark_data));
    unsigned int k = 0;
    const Box& it_box = mark_data->getBox();
    for (auto it = mark_data->begin(it_box), e = mark_data->end(it_box); it != e; ++it, ++k)
    {
        const LMarkerSet::value_type& mark = *it;
        const Point& X = mark->getPosition();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_mark[NDIM * k + d] = X[d];
        }
    }
    return;
}

void LMarkerUtilities::resetMarkerPositionsOnPatch(const std::vector<double>& X_mark,
                                                   const boost::shared_ptr<LMarkerSetData>& mark_data)
{
    unsigned int k = 0;
    const Box& it_box = mark_data->getBox();
    for (auto it = mark_data->begin(it_box), e = mark_data->end(it_box); it != e; ++it, ++k)
    {
        const LMarkerSet::value_type& mark = *it;
        Point& X = mark->getPosition();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X[d] = X_mark[NDIM * k + d];
        }
    }
    return;
}

void LMarkerUtilities::collectMarkerVelocitiesOnPatch(std::vector<double>& U_mark,
                                                      const boost::shared_ptr<LMarkerSetData>& mark_data)
{
    U_mark.resize(NDIM * countMarkersOnPatch(mark_data));
    unsigned int k = 0;
    const Box& it_box = mark_data->getBox();
    for (auto it = mark_data->begin(it_box), e = mark_data->end(it_box); it != e; ++it, ++k)
    {
        const LMarkerSet::value_type& mark = *it;
        const Vector& U = mark->getVelocity();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            U_mark[NDIM * k + d] = U[d];
        }
    }
    return;
}

void LMarkerUtilities::resetMarkerVelocitiesOnPatch(const std::vector<double>& U_mark,
                                                    const boost::shared_ptr<LMarkerSetData>& mark_data)
{
    unsigned int k = 0;
    const Box& it_box = mark_data->getBox();
    for (auto it = mark_data->begin(it_box), e = mark_data->end(it_box); it != e; ++it, ++k)
    {
        const LMarkerSet::value_type& mark = *it;
        Vector& U = mark->getVelocity();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            U[d] = U_mark[NDIM * k + d];
        }
    }
    return;
}

void LMarkerUtilities::preventMarkerEscape(std::vector<double>& X_mark,
                                           const boost::shared_ptr<CartesianGridGeometry>& grid_geom)
{
    const IntVector& periodic_shift = grid_geom->getPeriodicShift(IntVector::getOne(DIM));
    if (periodic_shift.min() > 0) return;
    static const double edge_tol = sqrt(std::numeric_limits<double>::epsilon());
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    for (unsigned int k = 0; k < X_mark.size() / NDIM; ++k)
    {
        double* const X = &X_mark[NDIM * k];
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift[d] != 0) continue;
            X[d] = std::max(X[d], x_lower[d] + edge_tol);
            X[d] = std::min(X[d], x_upper[d] - edge_tol);
        }
    }
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
