// Filename: LMarkerUtilities.C
// Created on 28 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE.

#include "LMarkerUtilities.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// C++ STDLIB INCLUDES
#include <limits>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LMarkerCoarsen.h>
#include <ibtk/LMarkerRefine.h>
#include <ibtk/LMarkerSet.h>
#include <ibtk/LMarkerSetData.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CoarsenAlgorithm.h>
#include <RefineAlgorithm.h>
#include <SideData.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

int
LMarkerUtilities::readMarkerPositions(
    std::vector<blitz::TinyVector<double,NDIM> >& mark_init_posns,
    const std::string& mark_input_file_name,
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (mark_input_file_name.empty()) return 0;

    // Read in the initial marker positions.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int mpi_size = SAMRAI_MPI::getNodes();

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
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
                TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n  Premature end to input file encountered before line 1 of file " << mark_input_file_name << "\n");
            }
            else
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
                if (!(line_stream >> num_mark))
                {
                    TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n  Invalid entry in input file encountered on line 1 of file " << mark_input_file_name << "\n");
                }
            }

            if (num_mark <= 0)
            {
                TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n  Invalid entry in input file encountered on line 1 of file " << mark_input_file_name << "\n");
            }

            // Each successive line provides the initial position of each
            // marker in the input file.
            mark_init_posns.resize(num_mark);
            for (int k = 0; k < num_mark; ++k)
            {
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n  Premature end to input file encountered before line " << k+2 << " of file " << mark_input_file_name << "\n");
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (!(line_stream >> mark_init_posns[k][d]))
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n  Invalid entry in input file encountered on line " << k+2 << " of file " << mark_input_file_name << "\n");
                        }
                    }

                    // Ensure the initial marker position lies within the
                    // physical domain.
                    const double* const X = mark_init_posns[k].data();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (MathUtilities<double>::equalEps(X[d],grid_xLower[d]))
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker intersecting lower physical boundary.\n"
                                       << "  please ensure that all markers are within the computational domain."<< std::endl);
                        }
                        else if (X[d] <= grid_xLower[d])
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker below lower physical boundary\n"
                                       << "  please ensure that all markers are within the computational domain."<< std::endl);
                        }

                        if (MathUtilities<double>::equalEps(X[d],grid_xUpper[d]))
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker intersecting upper physical boundary.\n"
                                       << "  please ensure that all markers are within the computational domain."<< std::endl);
                        }
                        else if (X[d] >= grid_xUpper[d])
                        {
                            TBOX_ERROR("LMarkerUtilities::readMarkerPositions():\n"
                                       << "  encountered marker above upper physical boundary\n"
                                       << "  please ensure that all markers are within the computational domain."<< std::endl);
                        }
                    }
                }
            }
        }
    }
    return num_mark;
}// readMarkerPositions

void
LMarkerUtilities::advectMarkers(
    const int mark_current_idx,
    const int mark_new_idx,
    const int u_idx,
    const double dt,
    const std::string& weighting_fcn,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<PatchData<NDIM> > v_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM,double> > v_cc_data = v_data;
            Pointer<SideData<NDIM,double> > v_sc_data = v_data;
            const bool is_cc_data = !v_cc_data.isNull();
            const bool is_sc_data = !v_sc_data.isNull();
            Pointer<LMarkerSetData> mark_current_data = patch->getPatchData(mark_current_idx);
            Pointer<LMarkerSetData> mark_new_data     = patch->getPatchData(mark_new_idx);

            // Collect the current positions of all markers in the patch into a
            // vector.
            std::vector<double> X_mark;
#ifdef DEBUG_CHECK_ASSERTIONS
            std::vector<int> idx_mark;
#endif
            for (LMarkerSetData::DataIterator it = mark_current_data->data_begin(patch_box);
                 it != mark_current_data->data_end(); ++it)
            {
                const LMarker& mark = *it;
                const blitz::TinyVector<double,NDIM>& X = mark.getPosition();
                X_mark.insert(X_mark.end(), X.begin(), X.end());
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& idx = mark.getIndex();
                idx_mark.push_back(idx);
#endif
            }

            // Collect the predicted new positions of all markers in the patch
            // into a vector.
            std::vector<double> X_mark_new;
#ifdef DEBUG_CHECK_ASSERTIONS
            std::vector<int> idx_mark_new;
#endif
            for (LMarkerSetData::DataIterator it = mark_new_data->data_begin(patch_box);
                 it != mark_new_data->data_end(); ++it)
            {
                const LMarker& mark = *it;
                const blitz::TinyVector<double,NDIM>& X = mark.getPosition();
                X_mark_new.insert(X_mark_new.end(), X.begin(), X.end());
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& idx = mark.getIndex();
                idx_mark_new.push_back(idx);
#endif
            }

            // Compute X_mark(n+1/2) = 0.5*(X_mark(n+1) + X_mark(n)).
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
            if (is_cc_data) LEInteractor::interpolate(U_mark_half, NDIM, X_mark_half, NDIM, v_cc_data, patch, patch_box, weighting_fcn);
            if (is_sc_data) LEInteractor::interpolate(U_mark_half, NDIM, X_mark_half, NDIM, v_sc_data, patch, patch_box, weighting_fcn);

            // Compute X_mark(n+1) = X_mark(n) + dt*U_mark(n+1/2).
            for (unsigned k = 0; k < X_mark.size(); ++k)
            {
                X_mark_new[k] = X_mark[k] + dt*U_mark_half[k];
            }

            // Prevent markers from leaving the computational domain through
            // physical boundaries (but *not* through periodic boundaries).
            static const double edge_tol = sqrt(std::numeric_limits<double>::epsilon());
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
            const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
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
            for (LMarkerSetData::DataIterator it = mark_new_data->data_begin(patch_box);
                 it != mark_new_data->data_end(); ++it, ++marker_offset)
            {
                LMarker& mark = *it;
                blitz::TinyVector<double,NDIM>& X = mark.getPosition();
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = X_mark_new[NDIM*marker_offset+d];
                }
            }
        }
    }
    return;
}// advectMarkers

void
LMarkerUtilities::collectMarkersOnPatchHierarchy(
    const int mark_idx,
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    const int num_marks_before_coarsening = countMarkers(mark_idx,hierarchy);
    const int num_marks_before_coarsening_level_0 = countMarkers(mark_idx,hierarchy,0,0);

    // Collect all marker data on the patch hierarchy.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(mark_idx, var);
    int mark_scratch_idx = var_db->registerClonedPatchDataIndex(var, mark_idx);
    Pointer<CoarsenAlgorithm<NDIM> > mark_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    mark_coarsen_alg->registerCoarsen(mark_scratch_idx, mark_idx, new LMarkerCoarsen());
    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);

        // Allocate scratch data.
        coarser_level->allocatePatchData(mark_scratch_idx);

        // Coarsen fine data onto coarser level.
        CoarsenPatchStrategy<NDIM>* mark_coarsen_op = NULL;
        mark_coarsen_alg->createSchedule(coarser_level, level, mark_coarsen_op)->coarsenData();

        // Merge the coarsened fine data with the coarse data.
        for (PatchLevel<NDIM>::Iterator p(coarser_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = coarser_level->getPatch(p());
            Pointer<LMarkerSetData> mark_current_data = patch->getPatchData(mark_idx);
            Pointer<LMarkerSetData> mark_scratch_data = patch->getPatchData(mark_scratch_idx);
            for (LMarkerSetData::Iterator it(*mark_scratch_data); it; it++)
            {
                const Index<NDIM>& i = it.getIndex();
                if (!mark_current_data->isElement(i))
                {
                    mark_current_data->appendItemPointer(i, new LMarkerSet());
                }
                LMarkerSet& dst_mark_set = *(mark_current_data->getItem(i));
                const LMarkerSet& src_mark_set = it();
                dst_mark_set.insert(dst_mark_set.end(), src_mark_set.begin(), src_mark_set.end());
            }
        }

        // Clear the fine data.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LMarkerSetData> mark_current_data = patch->getPatchData(mark_idx);
            mark_current_data->removeAllItems();
        }

        // Deallocate scratch data.
        coarser_level->deallocatePatchData(mark_scratch_idx);
    }
    var_db->removePatchDataIndex(mark_scratch_idx);
    mark_scratch_idx = -1;

    // Ensure that the total number of markers is correct.
    const int num_marks_after_coarsening = countMarkers(mark_idx,hierarchy);
    const int num_marks_after_coarsening_level_0 = countMarkers(mark_idx,hierarchy,0,0);
    if (num_marks_before_coarsening != num_marks_after_coarsening || num_marks_before_coarsening != num_marks_after_coarsening_level_0)
    {
        TBOX_ERROR("LMarkerUtilities::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during collection to coarsest level\n"
                   << "  number of markers in hierarchy before collection to coarsest level = " << num_marks_before_coarsening << "\n"
                   << "  number of markers on level 0   before collection to coarsest level = " << num_marks_before_coarsening_level_0 << "\n"
                   << "  number of markers in hierarchy after  collection to coarsest level = " << num_marks_after_coarsening << "\n"
                   << "  number of markers on level 0   after  collection to coarsest level = " << num_marks_after_coarsening_level_0 << "\n");
    }

    // Reset the assignment of markers to Cartesian grid cells on the coarsest
    // level of the patch hierarchy.
    //
    // NOTE: It is important to do this only *after* collecting markers on the
    // patch hierarchy, because markers which have left a fine level through the
    // coarse-fine interface would otherwise be discarded by this procedure.
    Pointer<RefineAlgorithm<NDIM> > mark_level_fill_alg = new RefineAlgorithm<NDIM>();
    mark_level_fill_alg->registerRefine(mark_idx, mark_idx, mark_idx, NULL);
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(coarsest_ln);
    mark_level_fill_alg->createSchedule(level,NULL)->fillData(0.0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const Index<NDIM>& patch_lower = patch_box.lower();
        const Index<NDIM>& patch_upper = patch_box.upper();
        const double* const patchXLower = patch_geom->getXLower();
        const double* const patchXUpper = patch_geom->getXUpper();
        const double* const patchDx = patch_geom->getDx();

        Pointer<LMarkerSetData> mark_data = patch->getPatchData(mark_idx);
        Pointer<LMarkerSetData> mark_data_new = new LMarkerSetData(mark_data->getBox(), mark_data->getGhostCellWidth());
        for (LMarkerSetData::DataIterator it = mark_data->data_begin(mark_data->getGhostBox());
             it != mark_data->data_end(); ++it)
        {
            const LMarker& mark = *it;
            const blitz::TinyVector<double,NDIM>& X = mark.getPosition();
            const IntVector<NDIM>& offset = mark.getPeriodicOffset();
            double X_shifted[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X_shifted[d] = X[d] + double(offset(d))*patchDx[d];
            }
            const bool patch_owns_mark_at_new_loc =
                ((  patchXLower[0] <= X_shifted[0])&&(X_shifted[0] < patchXUpper[0]))
#if (NDIM > 1)
                &&((patchXLower[1] <= X_shifted[1])&&(X_shifted[1] < patchXUpper[1]))
#if (NDIM > 2)
                &&((patchXLower[2] <= X_shifted[2])&&(X_shifted[2] < patchXUpper[2]))
#endif
#endif
                ;
            if (patch_owns_mark_at_new_loc)
            {
                const Index<NDIM> i = IndexUtilities::getCellIndex(X_shifted, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
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
    const int num_marks_after_posn_reset = countMarkers(mark_idx,hierarchy);
    const int num_marks_after_posn_reset_level_0 = countMarkers(mark_idx,hierarchy,0,0);
    if (num_marks_before_coarsening != num_marks_after_posn_reset || num_marks_before_coarsening != num_marks_after_posn_reset_level_0)
    {
        TBOX_ERROR("LMarkerUtilities::collectMarkersOnPatchHierarchy()\n"
                   << "  number of marker particles changed during position reset on coarsest level\n"
                   << "  number of markers in hierarchy before position reset on coarsest level = " << num_marks_before_coarsening << "\n"
                   << "  number of markers in hierarchy after  position reset on coarsest level = " << num_marks_after_posn_reset << "\n"
                   << "  number of markers on level 0   after  position reset on coarsest level = " << num_marks_after_posn_reset_level_0 << "\n");
    }
    return;
}// collectMarkersOnPatchHierarchy

void
LMarkerUtilities::initializeMarkersOnLevel(
    const int mark_idx,
    const std::vector<blitz::TinyVector<double,NDIM> >& mark_init_posns,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level)
{
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // On the coarsest level of the patch hierarchy, copy marker data from the
    // old coarse level.  Otherwise, refine marker data from the coarsest level
    // of the patch hierarchy.
    if (initial_time && level_number == 0)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            const double* const patchDx = patch_geom->getDx();

            Pointer<LMarkerSetData> mark_data = patch->getPatchData(mark_idx);
            for (unsigned k = 0; k < mark_init_posns.size(); ++k)
            {
                const blitz::TinyVector<double,NDIM>& X = mark_init_posns[k];
                static const blitz::TinyVector<double,NDIM> U(0.0);
                const bool patch_owns_mark_at_loc =
                    ((  patchXLower[0] <= X[0])&&(X[0] < patchXUpper[0]))
#if (NDIM > 1)
                    &&((patchXLower[1] <= X[1])&&(X[1] < patchXUpper[1]))
#if (NDIM > 2)
                    &&((patchXLower[2] <= X[2])&&(X[2] < patchXUpper[2]))
#endif
#endif
                    ;
                if (patch_owns_mark_at_loc)
                {
                    const Index<NDIM> i = IndexUtilities::getCellIndex(X, patchXLower, patchXUpper, patchDx, patch_lower, patch_upper);
                    if (!mark_data->isElement(i))
                    {
                        mark_data->appendItemPointer(i, new LMarkerSet());
                    }
                    LMarkerSet& new_mark_set = *(mark_data->getItem(i));
                    new_mark_set.push_back(LMarker(k, X, U));
                }
            }
        }
    }
    else if (!old_level.isNull() && level_number == 0)
    {
        Pointer<RefineAlgorithm<NDIM> > copy_mark_alg = new RefineAlgorithm<NDIM>();
        copy_mark_alg->registerRefine(mark_idx, mark_idx, mark_idx, NULL);
        Pointer<PatchLevel<NDIM> > dst_level = level;
        Pointer<PatchLevel<NDIM> > src_level = old_level;
        RefinePatchStrategy<NDIM>* refine_mark_op = NULL;
        copy_mark_alg->createSchedule(dst_level, src_level, refine_mark_op)->fillData(0.0);
    }
    else if (level_number > 0)
    {
        Pointer<RefineAlgorithm<NDIM> > refine_mark_alg = new RefineAlgorithm<NDIM>();
        refine_mark_alg->registerRefine(mark_idx, mark_idx, mark_idx, new LMarkerRefine());
        RefinePatchStrategy<NDIM>* refine_mark_op = NULL;
        for (int ln = 1; ln <= level_number; ++ln)
        {
            Pointer<PatchLevel<NDIM> > dst_level = hierarchy->getPatchLevel(ln);
            refine_mark_alg->createSchedule(dst_level, NULL, ln-1, hierarchy, refine_mark_op)->fillData(0.0);
        }
    }
    return;
}// initializeMarkersOnLevel

void
LMarkerUtilities::pruneDuplicateMarkers(
    const int mark_idx,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    const int finest_hier_level_number = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= std::min(finest_ln,finest_hier_level_number-1); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > finer_level = hierarchy->getPatchLevel(ln+1);
        BoxArray<NDIM> refined_region_boxes = finer_level->getBoxes();
        const IntVector<NDIM>& ratio = finer_level->getRatioToCoarserLevel();
        refined_region_boxes.coarsen(ratio);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LMarkerSetData> mark_data = patch->getPatchData(mark_idx);
            const Box<NDIM>& ghost_box = mark_data->getGhostBox();
            for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
            {
                const Box<NDIM>& refined_box = refined_region_boxes[i];
                const Box<NDIM> intersection = ghost_box * refined_box;
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
LMarkerUtilities::countMarkers(
    const int mark_idx,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    int num_marks = 0;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<LMarkerSetData> mark_data = patch->getPatchData(mark_idx);
            for (LMarkerSetData::SetIterator it(*mark_data); it; it++)
            {
                if (patch_box.contains(it.getIndex()))
                {
                    num_marks += it.getItem().size();
                }
            }
        }
    }
    return SAMRAI_MPI::sumReduction(num_marks);
}// countMarkers

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
