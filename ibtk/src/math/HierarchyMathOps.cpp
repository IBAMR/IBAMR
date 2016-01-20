// Filename: HierarchyMathOps.cpp
// Created on 11 Jun 2003 by Boyce Griffith
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

#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArrayDataBasicOps.h"
#include "BasePatchLevel.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "CoarseFineBoundary.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "EdgeData.h"
#include "EdgeVariable.h"
#include "FaceData.h"
#include "FaceDataFactory.h"
#include "FaceGeometry.h"
#include "FaceVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "OuterfaceData.h"
#include "OuterfaceDataFactory.h"
#include "OuterfaceVariable.h"
#include "OutersideData.h"
#include "OutersideDataFactory.h"
#include "OutersideVariable.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SideData.h"
#include "SideDataFactory.h"
#include "SideGeometry.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/PatchMathOps.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define S_TO_C_INTERP_SPECIAL_FC IBTK_FC_FUNC(stocinterp2ndspecial2d, STOCINTERP2NDSPECIAL2D)
#endif
#if (NDIM == 3)
#define S_TO_C_INTERP_SPECIAL_FC IBTK_FC_FUNC(stocinterp2ndspecial3d, STOCINTERP2NDSPECIAL3D)
#endif

// Function interfaces
extern "C" {
void S_TO_C_INTERP_SPECIAL_FC(const int& direction,
                              double* U,
                              const int& U_gcw,
                              const double& alpha,
                              const double* v0,
                              const double* v1,
#if (NDIM == 3)
                              const double* v2,
#endif
                              const int& v_gcw,
                              const double& beta,
                              const double* W,
                              const int& W_gcw,
                              const int& ilower0,
                              const int& iupper0,
                              const int& ilower1,
                              const int& iupper1
#if (NDIM == 3)
                              ,
                              const int& ilower2,
                              const int& iupper2
#endif
                              );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyMathOps::HierarchyMathOps(const std::string& name,
                                   Pointer<PatchHierarchy<NDIM> > hierarchy,
                                   const int coarsest_ln,
                                   const int finest_ln,
                                   const std::string& coarsen_op_name)
    : d_object_name(name), d_hierarchy(), d_grid_geom(), d_coarsest_ln(coarsest_ln), d_finest_ln(finest_ln),
      d_fc_var(new FaceVariable<NDIM, double>(d_object_name + "::scratch_fc")),
      d_sc_var(new SideVariable<NDIM, double>(d_object_name + "::scratch_sc")),
      d_of_var(new OuterfaceVariable<NDIM, double>(d_object_name + "::scratch_of")),
      d_os_var(new OutersideVariable<NDIM, double>(d_object_name + "::scratch_os")), d_fc_idx(-1), d_sc_idx(-1),
      d_of_idx(-1), d_os_idx(-1), d_coarsen_op_name(coarsen_op_name), d_of_coarsen_op(), d_os_coarsen_op(),
      d_of_coarsen_alg(), d_os_coarsen_alg(), d_of_coarsen_scheds(), d_os_coarsen_scheds(), d_hier_cc_data_ops(),
      d_hier_fc_data_ops(), d_hier_sc_data_ops(), d_patch_math_ops(), d_context(),
      d_wgt_cc_var(new CellVariable<NDIM, double>(d_object_name + "::wgt_cc", 1)),
      d_wgt_fc_var(new FaceVariable<NDIM, double>(d_object_name + "::wgt_fc", 1)),
      d_wgt_sc_var(new SideVariable<NDIM, double>(d_object_name + "::wgt_sc", 1)), d_wgt_cc_idx(-1), d_wgt_fc_idx(-1),
      d_wgt_sc_idx(-1), d_volume(0.0)
{
    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");

    static const bool fine_boundary_represents_var = true;
    static const IntVector<NDIM> no_ghosts = 0;
    d_fc_var->setPatchDataFactory(new FaceDataFactory<NDIM, double>(1, no_ghosts, fine_boundary_represents_var));
    d_sc_var->setPatchDataFactory(new SideDataFactory<NDIM, double>(1, no_ghosts, fine_boundary_represents_var));

    d_of_var->setPatchDataFactory(new OuterfaceDataFactory<NDIM, double>(1));
    d_os_var->setPatchDataFactory(new OutersideDataFactory<NDIM, double>(1));

    static const IntVector<NDIM> ghosts = 1;

    if (var_db->checkVariableExists(d_fc_var->getName()))
    {
        d_fc_var = var_db->getVariable(d_fc_var->getName());
        d_fc_idx = var_db->mapVariableAndContextToIndex(d_fc_var, d_context);
    }
    else
    {
        d_fc_idx = var_db->registerVariableAndContext(d_fc_var, d_context, ghosts);
    }

    if (var_db->checkVariableExists(d_sc_var->getName()))
    {
        d_sc_var = var_db->getVariable(d_sc_var->getName());
        d_sc_idx = var_db->mapVariableAndContextToIndex(d_sc_var, d_context);
    }
    else
    {
        d_sc_idx = var_db->registerVariableAndContext(d_sc_var, d_context, ghosts);
    }

    if (var_db->checkVariableExists(d_of_var->getName()))
    {
        d_of_var = var_db->getVariable(d_of_var->getName());
        d_of_idx = var_db->mapVariableAndContextToIndex(d_of_var, d_context);
    }
    else
    {
        d_of_idx = var_db->registerVariableAndContext(d_of_var, d_context);
    }

    if (var_db->checkVariableExists(d_os_var->getName()))
    {
        d_os_var = var_db->getVariable(d_os_var->getName());
        d_os_idx = var_db->mapVariableAndContextToIndex(d_os_var, d_context);
    }
    else
    {
        d_os_idx = var_db->registerVariableAndContext(d_os_var, d_context);
    }

    if (var_db->checkVariableExists(d_wgt_cc_var->getName()))
    {
        d_wgt_cc_var = var_db->getVariable(d_wgt_cc_var->getName());
        d_wgt_cc_idx = var_db->mapVariableAndContextToIndex(d_wgt_cc_var, d_context);
    }
    else
    {
        d_wgt_cc_idx = var_db->registerVariableAndContext(d_wgt_cc_var, d_context);
    }

    if (var_db->checkVariableExists(d_wgt_fc_var->getName()))
    {
        d_wgt_fc_var = var_db->getVariable(d_wgt_fc_var->getName());
        d_wgt_fc_idx = var_db->mapVariableAndContextToIndex(d_wgt_fc_var, d_context);
    }
    else
    {
        d_wgt_fc_idx = var_db->registerVariableAndContext(d_wgt_fc_var, d_context);
    }

    if (var_db->checkVariableExists(d_wgt_sc_var->getName()))
    {
        d_wgt_sc_var = var_db->getVariable(d_wgt_sc_var->getName());
        d_wgt_sc_idx = var_db->mapVariableAndContextToIndex(d_wgt_sc_var, d_context);
    }
    else
    {
        d_wgt_sc_idx = var_db->registerVariableAndContext(d_wgt_sc_var, d_context);
    }

    // Set the patch hierarchy.
    setPatchHierarchy(hierarchy);
    if ((coarsest_ln < 0) || (finest_ln < 0))
    {
        if (d_hierarchy->getNumberOfLevels() == 0)
        {
            d_coarsest_ln = coarsest_ln;
            d_finest_ln = finest_ln;
        }
        else
        {
            resetLevels(0, d_hierarchy->getFinestLevelNumber());
        }
    }
    else
    {
        resetLevels(coarsest_ln, finest_ln);
    }
    return;
} // HierarchyMathOps

HierarchyMathOps::~HierarchyMathOps()
{
    // intentionally blank
    return;
} // ~HierarchyMathOps

void HierarchyMathOps::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    d_grid_geom = hierarchy->getGridGeometry();

    // Obtain the hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();

    Pointer<CellVariable<NDIM, double> > cc_var = new CellVariable<NDIM, double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);

    Pointer<FaceVariable<NDIM, double> > fc_var = new FaceVariable<NDIM, double>("fc_var");
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(fc_var, d_hierarchy, true);

    Pointer<SideVariable<NDIM, double> > sc_var = new SideVariable<NDIM, double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, d_hierarchy, true);

    // Reset the communications operators.
    resetCoarsenOperators();
    resetRefineOperators();
    return;
} // setPatchHierarchy

void HierarchyMathOps::resetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
    TBOX_ASSERT((coarsest_ln >= 0) && (finest_ln >= coarsest_ln) && (finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;

    // Reset the CoarsenSchedule vectors.
    d_of_coarsen_scheds.resize(d_finest_ln);
    d_os_coarsen_scheds.resize(d_finest_ln);
    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM> > src_level = d_hierarchy->getPatchLevel(dst_ln + 1);
        Pointer<PatchLevel<NDIM> > dst_level = d_hierarchy->getPatchLevel(dst_ln);
        d_of_coarsen_scheds[dst_ln] = d_of_coarsen_alg->createSchedule(dst_level, src_level);
        d_os_coarsen_scheds[dst_ln] = d_os_coarsen_alg->createSchedule(dst_level, src_level);
    }

    // Reset the hierarchy data ops.
    d_hier_cc_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_hier_fc_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_hier_sc_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);

    // Reset the cell weights.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_wgt_cc_idx))
        {
            level->allocatePatchData(d_wgt_cc_idx);
        }
        if (!level->checkAllocated(d_wgt_fc_idx))
        {
            level->allocatePatchData(d_wgt_fc_idx);
        }
        if (!level->checkAllocated(d_wgt_sc_idx))
        {
            level->allocatePatchData(d_wgt_sc_idx);
        }
    }

    // Each cell's weight is set to its cell volume, unless the cell is refined
    // on a finer level, in which case the weight is set to zero.  This insures
    // that no part of the physical domain is counted twice when discrete norms
    // and integrals are calculated on the entire hierarchy.
    //
    // Away from coarse-fine interfaces and boundaries of the computational
    // domain, each cell face's weight is set to the cell volume associated with
    // the level of the patch hierarchy.  Along coarse-fine interfaces or
    // physical boundaries, the weights associated with the cell faces are
    // modified so that the sum of the weights equals to volume of the
    // computational domain.
    ArrayDataBasicOps<NDIM, double> array_ops;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        BoxArray<NDIM> refined_region_boxes;

        if (ln < d_finest_ln)
        {
            Pointer<PatchLevel<NDIM> > next_finer_level = d_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
        }

        const IntVector<NDIM> max_gcw(1);
        const CoarseFineBoundary<NDIM> cf_bdry(*d_hierarchy, ln, max_gcw);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            const double* const dx = pgeom->getDx();
            const double cell_vol = dx[0] * dx[1]
#if (NDIM > 2)
                                    * dx[2]
#endif
                ;

            Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(d_wgt_cc_idx);
            wgt_cc_data->fillAll(cell_vol);

            Pointer<FaceData<NDIM, double> > wgt_fc_data = patch->getPatchData(d_wgt_fc_idx);
            wgt_fc_data->fillAll(cell_vol);

            Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(d_wgt_sc_idx);
            wgt_sc_data->fillAll(cell_vol);

            // Rescale values along the edges of the patches.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> face_lower_box = FaceGeometry<NDIM>::toFaceBox(patch_box, axis);
                face_lower_box.upper()(0) = face_lower_box.lower()(0);
                array_ops.scale(wgt_fc_data->getArrayData(axis), 0.5, wgt_fc_data->getArrayData(axis), face_lower_box);
                Box<NDIM> side_lower_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                side_lower_box.upper()(axis) = side_lower_box.lower()(axis);
                array_ops.scale(wgt_sc_data->getArrayData(axis), 0.5, wgt_sc_data->getArrayData(axis), side_lower_box);
            }
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> face_upper_box = FaceGeometry<NDIM>::toFaceBox(patch_box, axis);
                face_upper_box.lower()(0) = face_upper_box.upper()(0);
                array_ops.scale(wgt_fc_data->getArrayData(axis), 0.5, wgt_fc_data->getArrayData(axis), face_upper_box);
                Box<NDIM> side_upper_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                side_upper_box.lower()(axis) = side_upper_box.upper()(axis);
                array_ops.scale(wgt_sc_data->getArrayData(axis), 0.5, wgt_sc_data->getArrayData(axis), side_upper_box);
            }

            // Correct the values along coarse-fine interfaces.
            if (ln > d_coarsest_ln)
            {
                const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
                const int bdry_type = 1;
                const Array<BoundaryBox<NDIM> >& cf_bdry_boxes = cf_bdry.getBoundaries(p(), bdry_type);
                for (int k = 0; k < cf_bdry_boxes.getSize(); ++k)
                {
                    const Box<NDIM>& bdry_box = cf_bdry_boxes[k].getBox();
                    const unsigned int axis = cf_bdry_boxes[k].getLocationIndex() / 2;
                    const int lower_upper = cf_bdry_boxes[k].getLocationIndex() % 2;
                    if (!pgeom->getTouchesRegularBoundary(axis, lower_upper))
                    {
                        const double extra_vol = 0.5 * static_cast<double>(ratio(axis)) * cell_vol;

                        Box<NDIM> face_bdry_box = FaceGeometry<NDIM>::toFaceBox(bdry_box, axis);
                        array_ops.addScalar(wgt_fc_data->getArrayData(axis), wgt_fc_data->getArrayData(axis), extra_vol,
                                            face_bdry_box);

                        Box<NDIM> side_bdry_box = SideGeometry<NDIM>::toSideBox(bdry_box, axis);
                        array_ops.addScalar(wgt_sc_data->getArrayData(axis), wgt_sc_data->getArrayData(axis), extra_vol,
                                            side_bdry_box);
                    }
                }
            }

            // Zero-out weights within the refined region.
            if (ln < d_finest_ln)
            {
                const IntVector<NDIM>& periodic_shift = d_grid_geom->getPeriodicShift(level->getRatio());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        if (periodic_shift(axis) != 0)
                        {
                            for (int sgn = -1; sgn <= 1; sgn += 2)
                            {
                                IntVector<NDIM> periodic_offset = 0;
                                periodic_offset(axis) = sgn * periodic_shift(axis);
                                const Box<NDIM> refined_box =
                                    Box<NDIM>::shift(refined_region_boxes[i], periodic_offset);
                                const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                                if (!intersection.empty())
                                {
                                    wgt_cc_data->fillAll(0.0, intersection);
                                    wgt_fc_data->fillAll(0.0, intersection);
                                    wgt_sc_data->fillAll(0.0, intersection);
                                }
                            }
                        }
                    }
                    const Box<NDIM>& refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        wgt_cc_data->fillAll(0.0, intersection);
                        wgt_fc_data->fillAll(0.0, intersection);
                        wgt_sc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
    }

#if !defined(NDEBUG)
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        // Check for overlapping boxes on this level.
        //
        // (This is potentially fairly expensive and hence is only done when
        // assertion checking is active.)
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        BoxList<NDIM> boxes(level->getBoxes());

        std::vector<bool> patch_overlaps(boxes.getNumberOfItems());
        std::vector<bool>::size_type j, k;

        for (k = 0; k < patch_overlaps.size(); ++k)
        {
            patch_overlaps[k] = false;
        }
        k = 0;
        while (!boxes.isEmpty())
        {
            j = k + 1;
            Box<NDIM> tryme = boxes.getFirstItem();
            boxes.removeFirstItem();

            for (BoxList<NDIM>::Iterator ib(boxes); ib; ib++)
            {
                if (tryme.intersects(ib()))
                {
                    patch_overlaps[k] = true;
                    patch_overlaps[j] = true;
                }
                ++j;
            }
            ++k;
        }

        for (k = 0; k < patch_overlaps.size(); ++k)
        {
            if (patch_overlaps[k])
            {
                TBOX_ERROR(d_object_name << "::initializeLevelData()\n"
                                         << "  patch " << k << " overlaps another patch!\n");
            }
        }
    }
#endif

    // Compute the volume of the physical domain.
    const double volume_cc = d_hier_cc_data_ops->sumControlVolumes(d_wgt_cc_idx, d_wgt_cc_idx);
    const double volume_fc =
        d_hier_fc_data_ops->sumControlVolumes(d_wgt_fc_idx, d_wgt_fc_idx) / static_cast<double>(NDIM);
    const double volume_sc =
        d_hier_sc_data_ops->sumControlVolumes(d_wgt_sc_idx, d_wgt_sc_idx) / static_cast<double>(NDIM);
    if (!MathUtilities<double>::equalEps(volume_cc, volume_fc) ||
        !MathUtilities<double>::equalEps(volume_cc, volume_sc))
    {
        pout << "WARNING:\n"
             << "HierarchyMathOps::resetLevels():\n"
             << "  cell-centered, face-centered, and side-centered weights give different "
                "total "
                "volumes of the computational domain:\n"
             << "      volume (cell-centered): " << volume_cc << "\n"
             << "      volume (face-centered): " << volume_fc << "\n"
             << "      volume (side-centered): " << volume_sc << std::endl;
    }
    d_volume = volume_cc;
    return;
} // resetLevels

Pointer<CellVariable<NDIM, double> > HierarchyMathOps::getCellWeightVariable() const
{
    return d_wgt_cc_var;
} // getCellWeightVariable

int HierarchyMathOps::getCellWeightPatchDescriptorIndex() const
{
    return d_wgt_cc_idx;
} // getCellWeightPatchDescriptorIndex

Pointer<FaceVariable<NDIM, double> > HierarchyMathOps::getFaceWeightVariable() const
{
    return d_wgt_fc_var;
} // getFaceWeightVariable

int HierarchyMathOps::getFaceWeightPatchDescriptorIndex() const
{
    return d_wgt_fc_idx;
} // getFaceWeightPatchDescriptorIndex

Pointer<SideVariable<NDIM, double> > HierarchyMathOps::getSideWeightVariable() const
{
    return d_wgt_sc_var;
} // getSideWeightVariable

int HierarchyMathOps::getSideWeightPatchDescriptorIndex() const
{
    return d_wgt_sc_idx;
} // getSideWeightPatchDescriptorIndex

double HierarchyMathOps::getVolumeOfPhysicalDomain() const
{
    return d_volume;
} // getVolumeOfPhysicalDomain

void HierarchyMathOps::setCoarsenOperatorName(const std::string& coarsen_op_name)
{
    d_coarsen_op_name = coarsen_op_name;
    resetCoarsenOperators();
    return;
} // setCoarsenOperatorName

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<CellVariable<NDIM, double> > src_var,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    if ((d_coarsest_ln == d_finest_ln) && (d_finest_ln == 0))
    {
        const int ln = d_coarsest_ln;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    else
    {
        // Compute the side centered gradient and interpolate.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sc_idx);
        }

        d_hier_cc_data_ops->setToScalar(dst_idx, 0.0, false);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            grad(d_sc_idx, d_sc_var,
                 true, // synch coarse-fine boundary
                 1.0, src_idx, src_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
                 Pointer<SideVariable<NDIM, double> >(NULL), d);

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());

                    Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
                    Pointer<SideData<NDIM, double> > sc_data = patch->getPatchData(d_sc_idx);
#if (NDIM == 2)
                    double* const W = dst_data->getPointer(0);
                    const int W_ghosts = (dst_data->getGhostCellWidth()).max();

                    const double* const g0 = sc_data->getPointer(0);
                    const double* const g1 = sc_data->getPointer(1);
                    const int g_ghosts = (sc_data->getGhostCellWidth()).max();

                    const Box<NDIM>& patch_box = patch->getBox();

                    const int direction = (d == 0) ? 1 : 0;
                    const double alpha = (d == 0) ? -1.0 : 1.0;
                    const double beta = (d == 0) ? 0.0 : 1.0;

                    S_TO_C_INTERP_SPECIAL_FC(direction, W, W_ghosts, alpha, g0, g1, g_ghosts, beta, W, W_ghosts,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1));
#endif
#if (NDIM == 3)
                    const int W_ghosts = (dst_data->getGhostCellWidth()).max();

                    const double* const g0 = sc_data->getPointer(0);
                    const double* const g1 = sc_data->getPointer(1);
                    const double* const g2 = sc_data->getPointer(2);
                    const int g_ghosts = (sc_data->getGhostCellWidth()).max();

                    const Box<NDIM>& patch_box = patch->getBox();

                    int dst_depth, direction;
                    static const double alpha0 = 1.0;
                    static const double beta = 1.0;

                    if (d == 0)
                    {
                        dst_depth = 1;
                        direction = 2;
                    }
                    else if (d == 1)
                    {
                        dst_depth = 2;
                        direction = 0;
                    }
                    else // (d == 2)
                    {
                        dst_depth = 0;
                        direction = 1;
                    }

                    S_TO_C_INTERP_SPECIAL_FC(direction, dst_data->getPointer(dst_depth), W_ghosts, alpha0, g0, g1, g2,
                                             g_ghosts, beta, dst_data->getPointer(dst_depth), W_ghosts,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1), patch_box.lower(2), patch_box.upper(2));

                    static const double alpha1 = -1.0;

                    if (d == 0)
                    {
                        dst_depth = 2;
                        direction = 1;
                    }
                    else if (d == 1)
                    {
                        dst_depth = 0;
                        direction = 2;
                    }
                    else // (d == 2)
                    {
                        dst_depth = 1;
                        direction = 0;
                    }

                    S_TO_C_INTERP_SPECIAL_FC(direction, dst_data->getPointer(dst_depth), W_ghosts, alpha1, g0, g1, g2,
                                             g_ghosts, beta, dst_data->getPointer(dst_depth), W_ghosts,
                                             patch_box.lower(0), patch_box.upper(0), patch_box.lower(1),
                                             patch_box.upper(1), patch_box.lower(2), patch_box.upper(2));
#endif
                }
            }
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sc_idx);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
#if (NDIM != 2)
    TBOX_ERROR("HierarchyMathOps::curl():\n"
               << "  not implemented for NDIM != 2" << std::endl);
#endif
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::curl(const int dst_idx,
                            const Pointer<EdgeVariable<NDIM, double> > /*dst_var*/,
                            const int src_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                            const double src_ghost_fill_time)
{
#if (NDIM != 3)
    TBOX_ERROR("HierarchyMathOps::curl():\n"
               << "  not implemented for NDIM != 3" << std::endl);
#endif
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete curl.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<EdgeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.curl(dst_data, src_data, patch);
        }
    }
    return;
} // curl

void HierarchyMathOps::rot(int dst_idx,
                           Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                           int src_idx,
                           Pointer<NodeVariable<NDIM, double> > /*src_var*/,
                           Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                           double src_ghost_fill_time,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if (NDIM != 2)
    TBOX_ERROR("HierarchyMathOps::rot():\n"
               << "  not implemented for NDIM != 2" << std::endl);
#endif
    CartSideRobinPhysBdryOp robin_bc_op;
    const bool has_bc_coefs = !bc_coefs.empty();
    if (has_bc_coefs)
    {
        robin_bc_op.setPatchDataIndex(dst_idx);
        robin_bc_op.setPhysicalBcCoefs(bc_coefs);
        robin_bc_op.setHomogeneousBc(true);
    }

    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete rot.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.rot(dst_data, src_data, patch, has_bc_coefs ? &robin_bc_op : NULL, src_ghost_fill_time);
        }
    }
    return;
} // rot

void HierarchyMathOps::rot(int dst_idx,
                           Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                           int src_idx,
                           Pointer<CellVariable<NDIM, double> > /*src_var*/,
                           Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                           double src_ghost_fill_time,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if (NDIM != 2)
    TBOX_ERROR("HierarchyMathOps::rot():\n"
               << "  not implemented for NDIM != 2" << std::endl);
#endif
    CartSideRobinPhysBdryOp robin_bc_op;
    const bool has_bc_coefs = !bc_coefs.empty();
    if (has_bc_coefs)
    {
        robin_bc_op.setPatchDataIndex(dst_idx);
        robin_bc_op.setPhysicalBcCoefs(bc_coefs);
        robin_bc_op.setHomogeneousBc(true);
    }

    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete rot.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.rot(dst_data, src_data, patch, has_bc_coefs ? &robin_bc_op : NULL, src_ghost_fill_time);
        }
    }
    return;
} // rot

void HierarchyMathOps::rot(int dst_idx,
                           Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                           int src_idx,
                           Pointer<EdgeVariable<NDIM, double> > /*src_var*/,
                           Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                           double src_ghost_fill_time,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if (NDIM != 3)
    TBOX_ERROR("HierarchyMathOps::rot():\n"
               << "  not implemented for NDIM != 3" << std::endl);
#endif
    CartSideRobinPhysBdryOp robin_bc_op;
    const bool has_bc_coefs = !bc_coefs.empty();
    if (has_bc_coefs)
    {
        robin_bc_op.setPatchDataIndex(dst_idx);
        robin_bc_op.setPhysicalBcCoefs(bc_coefs);
        robin_bc_op.setHomogeneousBc(true);
    }

    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete rot.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<EdgeData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.rot(dst_data, src_data, patch, has_bc_coefs ? &robin_bc_op : NULL, src_ghost_fill_time);
        }
    }
    return;
} // rot

void HierarchyMathOps::rot(int dst_idx,
                           Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                           int src_idx,
                           Pointer<SideVariable<NDIM, double> > /*src_var*/,
                           Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                           double src_ghost_fill_time,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    CartSideRobinPhysBdryOp robin_bc_op;
    const bool has_bc_coefs = !bc_coefs.empty();
    if (has_bc_coefs)
    {
        robin_bc_op.setPatchDataIndex(dst_idx);
        robin_bc_op.setPhysicalBcCoefs(bc_coefs);
        robin_bc_op.setHomogeneousBc(true);
    }

    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete rot.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.rot(dst_data, src_data, patch, has_bc_coefs ? &robin_bc_op : NULL, src_ghost_fill_time);
        }
    }
    return;
} // rot

void HierarchyMathOps::div(const int dst_idx,
                           const Pointer<CellVariable<NDIM, double> > dst_var,
                           const double alpha,
                           const int src1_idx,
                           const Pointer<CellVariable<NDIM, double> > src1_var,
                           const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                           const double src1_ghost_fill_time,
                           const double beta,
                           const int src2_idx,
                           const Pointer<CellVariable<NDIM, double> > src2_var,
                           const int dst_depth,
                           const int src2_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    if ((d_coarsest_ln == d_finest_ln) && (d_finest_ln == 0))
    {
        const int ln = d_finest_ln;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete divergence.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.div(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth, src2_depth);
        }
    }
    else
    {
        // Interpolate to a side centered variable and compute the divergence of
        // the interpolated data.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sc_idx);
        }

        interp(d_sc_idx, d_sc_var,
               true, // synch coarse-fine boundary
               src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0);

        div(dst_idx, dst_var, alpha, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
            false, // don't re-synch cf boundary
            beta, src2_idx, src2_var, dst_depth, src2_depth);

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sc_idx);
        }
    }
    return;
} // div

void HierarchyMathOps::div(const int dst_idx,
                           const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                           const double alpha,
                           const int src1_idx,
                           const Pointer<FaceVariable<NDIM, double> > /*src1_var*/,
                           const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                           const double src1_ghost_fill_time,
                           const bool src1_cf_bdry_synch,
                           const double beta,
                           const int src2_idx,
                           const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                           const int dst_depth,
                           const int src2_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
        {
            level->allocatePatchData(d_of_idx);
        }

        // Compute the discrete divergence and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.div(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth, src2_depth);

            if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
            {
                Pointer<OuterfaceData<NDIM, double> > of_data = patch->getPatchData(d_of_idx);
                of_data->copy(*src1_data);
            }
        }

        // Synchronize the coarse-fine interface of src1 and deallocate
        // temporary data.
        if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
        {
            xeqScheduleOuterfaceRestriction(src1_idx, d_of_idx, ln - 1);
            level->deallocatePatchData(d_of_idx);
        }
    }
    return;
} // div

void HierarchyMathOps::div(const int dst_idx,
                           const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                           const double alpha,
                           const int src1_idx,
                           const Pointer<SideVariable<NDIM, double> > /*src1_var*/,
                           const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                           const double src1_ghost_fill_time,
                           const bool src1_cf_bdry_synch,
                           const double beta,
                           const int src2_idx,
                           const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                           const int dst_depth,
                           const int src2_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
        {
            level->allocatePatchData(d_os_idx);
        }

        // Compute the discrete divergence and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.div(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth, src2_depth);

            if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
            {
                Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*src1_data);
            }
        }

        // Synchronize the coarse-fine interface of src1 and deallocate
        // temporary data.
        if ((ln > d_coarsest_ln) && src1_cf_bdry_synch)
        {
            xeqScheduleOutersideRestriction(src1_idx, d_os_idx, ln - 1);
            level->deallocatePatchData(d_os_idx);
        }
    }
    return;
} // div

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > dst_var,
                            const double alpha,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > src1_var,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    if ((d_coarsest_ln == d_finest_ln) && (d_finest_ln == 0))
    {
        const int ln = d_finest_ln;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete gradient.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.grad(dst_data, alpha, src1_data, beta, src2_data, patch, src1_depth);
        }
    }
    else
    {
        // Compute the side centered gradient and interpolate.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sc_idx);
        }

        grad(d_sc_idx, d_sc_var,
             true, // synch coarse-fine boundary
             alpha, src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
             Pointer<SideVariable<NDIM, double> >(NULL), src1_depth);

        if (beta != 0.0)
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            int cc_idx = var_db->registerClonedPatchDataIndex(dst_var, dst_idx);
            const Pointer<CellVariable<NDIM, double> > cc_var = dst_var;

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                d_hierarchy->getPatchLevel(ln)->allocatePatchData(cc_idx);
            }

            interp(cc_idx, cc_var, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                   false); // don't re-synch cf boundary

            d_hier_cc_data_ops->linearSum(dst_idx,   // dst
                                          1.0,       // alpha
                                          cc_idx,    // src1
                                          beta,      // beta
                                          src2_idx); // src2

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                d_hierarchy->getPatchLevel(ln)->deallocatePatchData(cc_idx);
            }

            var_db->removePatchDataIndex(cc_idx);
            cc_idx = -1;
        }
        else
        {
            interp(dst_idx, dst_var, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                   false); // don't re-synch cf boundary
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sc_idx);
        }
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                            const bool dst_cf_bdry_synch,
                            const double alpha,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_of_idx);
        }

        // Compute the discrete gradient and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<FaceData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.grad(dst_data, alpha, src1_data, beta, src2_data, patch, src1_depth);

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OuterfaceData<NDIM, double> > of_data = patch->getPatchData(d_of_idx);
                of_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOuterfaceRestriction(dst_idx, d_of_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_of_idx);
        }
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                            const bool dst_cf_bdry_synch,
                            const double alpha,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_os_idx);
        }

        // Compute the discrete gradient and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.grad(dst_data, alpha, src1_data, beta, src2_data, patch, src1_depth);

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOutersideRestriction(dst_idx, d_os_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_os_idx);
        }
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > dst_var,
                            const int alpha_idx,
                            const Pointer<FaceVariable<NDIM, double> > alpha_var,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > src1_var,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    // Compute the face centered gradient and interpolate.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_fc_idx);
    }

    grad(d_fc_idx, d_fc_var,
         true, // synch coarse-fine boundary
         alpha_idx, alpha_var, src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
         Pointer<FaceVariable<NDIM, double> >(NULL), src1_depth);

    if (beta != 0.0)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int cc_idx = var_db->registerClonedPatchDataIndex(dst_var, dst_idx);
        const Pointer<CellVariable<NDIM, double> > cc_var = dst_var;

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(cc_idx);
        }

        interp(cc_idx, cc_var, d_fc_idx, d_fc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
               false); // don't re-synch cf boundary

        d_hier_cc_data_ops->linearSum(dst_idx,   // dst
                                      1.0,       // alpha
                                      cc_idx,    // src1
                                      beta,      // beta
                                      src2_idx); // src2

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(cc_idx);
        }

        var_db->removePatchDataIndex(cc_idx);
        cc_idx = -1;
    }
    else
    {
        interp(dst_idx, dst_var, d_fc_idx, d_fc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
               false); // don't re-synch cf boundary
    }

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_fc_idx);
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<CellVariable<NDIM, double> > dst_var,
                            const int alpha_idx,
                            const Pointer<SideVariable<NDIM, double> > alpha_var,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > src1_var,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    // Compute the side centered gradient and interpolate.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sc_idx);
    }

    grad(d_sc_idx, d_sc_var,
         true, // synch coarse-fine boundary
         alpha_idx, alpha_var, src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
         Pointer<SideVariable<NDIM, double> >(NULL), src1_depth);

    if (beta != 0.0)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int cc_idx = var_db->registerClonedPatchDataIndex(dst_var, dst_idx);
        const Pointer<CellVariable<NDIM, double> > cc_var = dst_var;

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(cc_idx);
        }

        interp(cc_idx, cc_var, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
               false); // don't re-synch cf boundary

        d_hier_cc_data_ops->linearSum(dst_idx,   // dst
                                      1.0,       // alpha
                                      cc_idx,    // src1
                                      beta,      // beta
                                      src2_idx); // src2

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(cc_idx);
        }

        var_db->removePatchDataIndex(cc_idx);
        cc_idx = -1;
    }
    else
    {
        interp(dst_idx, dst_var, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
               false); // don't re-synch cf boundary
    }

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sc_idx);
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                            const bool dst_cf_bdry_synch,
                            const int alpha_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*alpha_var*/,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<FaceVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_of_idx);
        }

        // Compute the discrete gradient and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<FaceData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<FaceData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.grad(dst_data, alpha_data, src1_data, beta, src2_data, patch, src1_depth);

            // Zero-out data on physical boundaries in the case of non-grid
            // aligned anisotropy.  (This is equivalent to enforcing no-flux
            // boundary conditions at the physical boundary.)
            if (alpha_data->getDepth() > 1)
            {
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    static const int gcw = 1;
                    Box<NDIM> boundary_box = Box<NDIM>::grow(patch_box, gcw);
                    const unsigned int axis_lower = patch_box.lower()[axis];
                    const unsigned int axis_upper = patch_box.upper()[axis];
                    for (int upperlower = 0; upperlower <= 1; ++upperlower)
                    {
                        if (pgeom->getTouchesRegularBoundary(axis, upperlower))
                        {
#if !defined(NDEBUG)
                            TBOX_ASSERT(!pgeom->getTouchesPeriodicBoundary(axis, upperlower));
#endif
                            if (upperlower == 0)
                            {
                                boundary_box.lower()[axis] = axis_lower - gcw;
                                boundary_box.upper()[axis] = axis_lower - 1;
                            }
                            else
                            {
                                boundary_box.lower()[axis] = axis_upper + 1;
                                boundary_box.upper()[axis] = axis_upper + gcw;
                            }
                            dst_data->fill(0.0, boundary_box);
                        }
                    }
                }
            }

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OuterfaceData<NDIM, double> > of_data = patch->getPatchData(d_of_idx);
                of_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOuterfaceRestriction(dst_idx, d_of_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_of_idx);
        }
    }
    return;
} // grad

void HierarchyMathOps::grad(const int dst_idx,
                            const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                            const bool dst_cf_bdry_synch,
                            const int alpha_idx,
                            const Pointer<SideVariable<NDIM, double> > /*alpha_var*/,
                            const int src1_idx,
                            const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                            const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                            const double src1_ghost_fill_time,
                            const double beta,
                            const int src2_idx,
                            const Pointer<SideVariable<NDIM, double> > /*src2_var*/,
                            const int src1_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_os_idx);
        }

        // Compute the discrete gradient and extract data on the coarse-fine
        // interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<SideData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.grad(dst_data, alpha_data, src1_data, beta, src2_data, patch, src1_depth);

            // Zero-out data on physical boundaries in the case of non-grid
            // aligned anisotropy.  (This is equivalent to enforcing no-flux
            // boundary conditions at the physical boundary.)
            if (alpha_data->getDepth() > 1)
            {
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    static const int gcw = 1;
                    Box<NDIM> boundary_box = Box<NDIM>::grow(patch_box, gcw);
                    const unsigned int axis_lower = patch_box.lower()[axis];
                    const unsigned int axis_upper = patch_box.upper()[axis];
                    for (int upperlower = 0; upperlower <= 1; ++upperlower)
                    {
                        if (pgeom->getTouchesRegularBoundary(axis, upperlower))
                        {
#if !defined(NDEBUG)
                            TBOX_ASSERT(!pgeom->getTouchesPeriodicBoundary(axis, upperlower));
#endif
                            if (upperlower == 0)
                            {
                                boundary_box.lower()[axis] = axis_lower - gcw;
                                boundary_box.upper()[axis] = axis_lower - 1;
                            }
                            else
                            {
                                boundary_box.lower()[axis] = axis_upper + 1;
                                boundary_box.upper()[axis] = axis_upper + gcw;
                            }
                            dst_data->fill(0.0, boundary_box);
                        }
                    }
                }
            }

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOutersideRestriction(dst_idx, d_os_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_os_idx);
        }
    }
    return;
} // grad

void HierarchyMathOps::interp(const int dst_idx,
                              const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                              const int src_idx,
                              const Pointer<FaceVariable<NDIM, double> > /*src_var*/,
                              const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                              const double src_ghost_fill_time,
                              const bool src_cf_bdry_synch)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
        {
            level->allocatePatchData(d_of_idx);
        }

        // Interpolate and extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.interp(dst_data, src_data, patch);

            if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
            {
                Pointer<OuterfaceData<NDIM, double> > of_data = patch->getPatchData(d_of_idx);
                of_data->copy(*src_data);
            }
        }

        // Synchronize the coarse-fine interface and deallocate temporary data.
        if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
        {
            xeqScheduleOuterfaceRestriction(src_idx, d_of_idx, ln - 1);
            level->deallocatePatchData(d_of_idx);
        }
    }
    return;
} // interp

void HierarchyMathOps::interp(const int dst_idx,
                              const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                              const int src_idx,
                              const Pointer<SideVariable<NDIM, double> > /*src_var*/,
                              const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                              const double src_ghost_fill_time,
                              const bool src_cf_bdry_synch)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
        {
            level->allocatePatchData(d_os_idx);
        }

        // Interpolate and extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.interp(dst_data, src_data, patch);

            if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
            {
                Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*src_data);
            }
        }

        // Synchronize the coarse-fine interface and deallocate temporary data.
        if ((ln > d_coarsest_ln) && src_cf_bdry_synch)
        {
            xeqScheduleOutersideRestriction(src_idx, d_os_idx, ln - 1);
            level->deallocatePatchData(d_os_idx);
        }
    }
    return;
} // interp

void HierarchyMathOps::interp(const int dst_idx,
                              const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                              const bool dst_cf_bdry_synch,
                              const int src_idx,
                              const Pointer<CellVariable<NDIM, double> > /*src_var*/,
                              const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                              const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_of_idx);
        }

        // Interpolate and extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.interp(dst_data, src_data, patch);

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OuterfaceData<NDIM, double> > of_data = patch->getPatchData(d_of_idx);
                of_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOuterfaceRestriction(dst_idx, d_of_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_of_idx);
        }
    }
    return;
} // interp

void HierarchyMathOps::interp(const int dst_idx,
                              const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                              const bool dst_cf_bdry_synch,
                              const int src_idx,
                              const Pointer<CellVariable<NDIM, double> > /*src_var*/,
                              const Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                              const double src_ghost_fill_time)
{
    if (src_ghost_fill) src_ghost_fill->fillData(src_ghost_fill_time);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Allocate temporary data to synchronize the coarse-fine interface.
        if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
        {
            level->allocatePatchData(d_os_idx);
        }

        // Interpolate and extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.interp(dst_data, src_data, patch);

            if ((ln > d_coarsest_ln) && dst_cf_bdry_synch)
            {
                Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
                os_data->copy(*dst_data);
            }
        }
    }

    // Synchronize the coarse-fine interface and deallocate temporary data.
    if (dst_cf_bdry_synch)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln + 1; --ln)
        {
            xeqScheduleOutersideRestriction(dst_idx, d_os_idx, ln - 1);
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_os_idx);
        }
    }
    return;
} // interp

void HierarchyMathOps::laplace(const int dst_idx,
                               const Pointer<CellVariable<NDIM, double> > dst_var,
                               const PoissonSpecifications& poisson_spec,
                               const int src1_idx,
                               const Pointer<CellVariable<NDIM, double> > src1_var,
                               const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                               const double src1_ghost_fill_time,
                               const double gamma,
                               const int src2_idx,
                               const Pointer<CellVariable<NDIM, double> > src2_var,
                               const int dst_depth,
                               const int src1_depth,
                               const int src2_depth)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    const double alpha = poisson_spec.dIsConstant() ? poisson_spec.getDConstant() : 0.0;
    const double beta = poisson_spec.cIsConstant() ? poisson_spec.getCConstant() : 0.0;

    const int alpha_idx = (poisson_spec.dIsConstant()) ? -1 : poisson_spec.getDPatchDataId();
    const int beta_idx = (poisson_spec.cIsConstant() || poisson_spec.cIsZero()) ? -1 : poisson_spec.getCPatchDataId();

    Pointer<SideVariable<NDIM, double> > alpha_var;
    Pointer<CellVariable<NDIM, double> > beta_var;

    bool nonaligned_anisotropy = false;
    if (!poisson_spec.dIsConstant())
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM> > dummy_var;
        var_db->mapIndexToVariable(alpha_idx, dummy_var);
        alpha_var = dummy_var;
#if !defined(NDEBUG)
        TBOX_ASSERT(alpha_var);
#endif
        Pointer<SideDataFactory<NDIM, double> > alpha_fac =
            var_db->getPatchDescriptor()->getPatchDataFactory(alpha_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(alpha_fac);
#endif
        nonaligned_anisotropy = alpha_fac->getDefaultDepth() > 1;
    }

    if (!(poisson_spec.cIsConstant() || poisson_spec.cIsZero()))
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM> > dummy_var;
        var_db->mapIndexToVariable(beta_idx, dummy_var);
        beta_var = dummy_var;
#if !defined(NDEBUG)
        TBOX_ASSERT(beta_var);
#endif
    }

    if ((d_coarsest_ln == d_finest_ln) && (alpha_idx == -1) && (!nonaligned_anisotropy))
    {
        // Compute dst = div alpha grad src1 + beta src1 + gamma src2.
        const int ln = d_finest_ln;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Compute the discrete Laplacian.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.laplace(dst_data, alpha, beta, src1_data, gamma, src2_data, patch, dst_depth, src1_depth,
                                     src2_depth);
        }
    }
    else
    {
        // Allocate temporary data.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_sc_idx);
        }

        // Compute the side centered normal flux of src1[m(i)] and put the
        // result in sc_var.
        if (alpha_idx == -1)
        {
            grad(d_sc_idx, d_sc_var,
                 true, // synch coarse-fine boundary
                 alpha, src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
                 Pointer<SideVariable<NDIM, double> >(NULL), src1_depth);
        }
        else
        {
            grad(d_sc_idx, d_sc_var,
                 true, // synch coarse-fine boundary
                 alpha_idx, alpha_var, src1_idx, src1_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0, 0.0, -1,
                 Pointer<SideVariable<NDIM, double> >(NULL), src1_depth);
        }

        // Take the divergence of the flux.
        if (MathUtilities<double>::equalEps(beta, 0.0) && MathUtilities<double>::equalEps(gamma, 0.0))
        {
            div(dst_idx, dst_var, 1.0, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                false, // don't re-synch coarse-fine boundary
                0.0, -1, Pointer<CellVariable<NDIM, double> >(NULL), dst_depth);
        }
        else if (MathUtilities<double>::equalEps(beta, 0.0))
        {
            div(dst_idx, dst_var, 1.0, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                false, // don't re-synch coarse-fine boundary
                gamma, src2_idx, src2_var, dst_depth, src2_depth);
        }
        else if (MathUtilities<double>::equalEps(gamma, 0.0))
        {
            div(dst_idx, dst_var, 1.0, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                false, // don't re-synch coarse-fine boundary
                beta, src1_idx, src1_var, dst_depth, src1_depth);
        }
        else
        {
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            int cc_idx = var_db->registerClonedPatchDataIndex(dst_var, dst_idx);
            const Pointer<CellVariable<NDIM, double> > cc_var = dst_var;
            const int cc_depth = dst_depth;

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                d_hierarchy->getPatchLevel(ln)->allocatePatchData(cc_idx);
            }

            div(cc_idx, cc_var, 1.0, d_sc_idx, d_sc_var, Pointer<HierarchyGhostCellInterpolation>(NULL), 0.0,
                false, // don't re-synch coarse-fine boundary
                beta, src1_idx, src1_var, cc_depth, src1_depth);

            pointwiseMultiply(dst_idx, dst_var, gamma, src2_idx, src2_var, 1.0, cc_idx, cc_var, dst_depth, src2_depth,
                              cc_depth);

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                d_hierarchy->getPatchLevel(ln)->deallocatePatchData(cc_idx);
            }

            var_db->removePatchDataIndex(cc_idx);
            cc_idx = -1;
        }

        // Deallocate temporary data.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_sc_idx);
        }
    }

    // Take care of the case where beta is spatially varying.
    if (beta_idx != -1)
    {
        pointwiseMultiply(dst_idx, dst_var, beta_idx, beta_var, src1_idx, src1_var, 1.0, dst_idx, dst_var, dst_depth,
                          src1_depth, dst_depth);
    }
    return;
} // laplace

void HierarchyMathOps::laplace(const int dst_idx,
                               const Pointer<SideVariable<NDIM, double> > dst_var,
                               const PoissonSpecifications& poisson_spec,
                               const int src1_idx,
                               const Pointer<SideVariable<NDIM, double> > src1_var,
                               const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                               const double src1_ghost_fill_time,
                               const double gamma,
                               const int src2_idx,
                               const Pointer<SideVariable<NDIM, double> > src2_var)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    const double alpha = poisson_spec.dIsConstant() ? poisson_spec.getDConstant() : 0.0;
    const double beta = poisson_spec.cIsConstant() ? poisson_spec.getCConstant() : 0.0;

    const int alpha_idx = (poisson_spec.dIsConstant()) ? -1 : poisson_spec.getDPatchDataId();
    const int beta_idx = (poisson_spec.cIsConstant() || poisson_spec.cIsZero()) ? -1 : poisson_spec.getCPatchDataId();

    if (alpha_idx != -1)
    {
        TBOX_ERROR("HierarchyMathOps::laplace():\n"
                   << "  side-centered Laplacian requires spatially constant scalar-valued "
                      "diffusivity" << std::endl);
    }

    if (beta_idx != -1)
    {
        TBOX_ERROR("HierarchyMathOps::laplace():\n"
                   << "  side-centered Laplacian requires spatially constant scalar-valued "
                      "damping factor" << std::endl);
    }

    if (!src1_var->fineBoundaryRepresentsVariable())
    {
        TBOX_WARNING("HierarchyMathOps::laplace():\n"
                     << "  recommended usage for side-centered Laplace operator is\n"
                     << "  src1_var->fineBoundaryRepresentsVariable() == true" << std::endl);
    }

    Pointer<SideDataFactory<NDIM, double> > dst_factory = dst_var->getPatchDataFactory();
    Pointer<SideDataFactory<NDIM, double> > src1_factory = src1_var->getPatchDataFactory();
    if (dst_factory->getDefaultDepth() != 1 || src1_factory->getDefaultDepth() != 1)
    {
        TBOX_ERROR("HierarchyMathOps::laplace():\n"
                   << "  side-centered Laplacian requires scalar-valued data" << std::endl);
    }
    if (src2_var)
    {
        Pointer<SideDataFactory<NDIM, double> > src2_factory = src2_var->getPatchDataFactory();
        if (src2_factory->getDefaultDepth() != 1)
        {
            TBOX_ERROR("HierarchyMathOps::laplace():\n"
                       << "  side-centered Laplacian requires scalar-valued data" << std::endl);
        }
    }

    // Compute dst = div grad src1 independently on each level.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.laplace(dst_data, alpha, beta, src1_data, gamma, src2_data, patch);
        }
    }

    // Allocate temporary data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_os_idx);
    }

    // Synchronize data along the coarse-fine interface.
    for (int ln = d_finest_ln; ln > d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
            os_data->copy(*dst_data);
        }

        // Synchronize the coarse-fine interface of dst.
        xeqScheduleOutersideRestriction(dst_idx, d_os_idx, ln - 1);
    }

    // Deallocate temporary data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_os_idx);
    }
    return;
} // laplace

void HierarchyMathOps::vc_laplace(const int dst_idx,
                                  const Pointer<SideVariable<NDIM, double> > dst_var,
                                  const double alpha,
                                  const double beta,
                                  const int coef_idx,
                                  const Pointer<NodeVariable<NDIM, double> > /*coef_var*/,
                                  const int src1_idx,
                                  const Pointer<SideVariable<NDIM, double> > src1_var,
                                  const Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                                  const double src1_ghost_fill_time,
                                  const double gamma,
                                  const int src2_idx,
                                  const Pointer<SideVariable<NDIM, double> > src2_var)
{
    if (src1_ghost_fill) src1_ghost_fill->fillData(src1_ghost_fill_time);

    Pointer<SideDataFactory<NDIM, double> > dst_factory = dst_var->getPatchDataFactory();
    Pointer<SideDataFactory<NDIM, double> > src1_factory = src1_var->getPatchDataFactory();
    if (dst_factory->getDefaultDepth() != 1 || src1_factory->getDefaultDepth() != 1)
    {
        TBOX_ERROR("HierarchyMathOps::vc_laplace():\n"
                   << "  side-centered variable-coefficient Laplacian requires scalar-valued data" << std::endl);
    }
    if (src2_var)
    {
        Pointer<SideDataFactory<NDIM, double> > src2_factory = src2_var->getPatchDataFactory();
        if (src2_factory->getDefaultDepth() != 1)
        {
            TBOX_ERROR("HierarchyMathOps::vc_laplace():\n"
                       << "  side-centered variable-coefficient Laplacian requires scalar-valued data" << std::endl);
        }
    }

    // Compute dst = alpha div coef ((grad src1) + (grad src1)^T) + beta src1 +
    // gamma src2 independently on each level.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > coef_data = patch->getPatchData(coef_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.vc_laplace(dst_data, alpha, beta, coef_data, src1_data, gamma, src2_data, patch);
        }
    }

    // Allocate temporary data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_os_idx);
    }

    // Synchronize data along the coarse-fine interface.
    for (int ln = d_finest_ln; ln > d_coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Extract data on the coarse-fine interface.
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<OutersideData<NDIM, double> > os_data = patch->getPatchData(d_os_idx);
            os_data->copy(*dst_data);
        }

        // Synchronize the coarse-fine interface of dst.
        xeqScheduleOutersideRestriction(dst_idx, d_os_idx, ln - 1);
    }

    // Deallocate temporary data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_os_idx);
    }
    return;
} // vc_laplace

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                         const double alpha,
                                         const int src1_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<CellData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src1_var*/,
                                         const int beta_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*beta_var*/,
                                         const int src2_idx,
                                         const Pointer<CellVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth,
                                         const int beta_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<CellData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<CellData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);
            Pointer<CellData<NDIM, double> > beta_data = patch->getPatchData(beta_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta_data, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth, beta_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                                         const double alpha,
                                         const int src1_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<FaceData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<FaceData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<FaceData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src1_var*/,
                                         const int beta_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*beta_var*/,
                                         const int src2_idx,
                                         const Pointer<FaceVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth,
                                         const int beta_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<FaceData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<FaceData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<FaceData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<FaceData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);
            Pointer<FaceData<NDIM, double> > beta_data = patch->getPatchData(beta_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta_data, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth, beta_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                         const double alpha,
                                         const int src1_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<NodeData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<NodeData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<NodeData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src1_var*/,
                                         const int beta_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*beta_var*/,
                                         const int src2_idx,
                                         const Pointer<NodeVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth,
                                         const int beta_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<NodeData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<NodeData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);
            Pointer<NodeData<NDIM, double> > beta_data = patch->getPatchData(beta_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta_data, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth, beta_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                                         const double alpha,
                                         const int src1_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src1_var*/,
                                         const double beta,
                                         const int src2_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<SideData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseMultiply(const int dst_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*dst_var*/,
                                         const int alpha_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*alpha_var*/,
                                         const int src1_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src1_var*/,
                                         const int beta_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*beta_var*/,
                                         const int src2_idx,
                                         const Pointer<SideVariable<NDIM, double> > /*src2_var*/,
                                         const int dst_depth,
                                         const int src1_depth,
                                         const int src2_depth,
                                         const int alpha_depth,
                                         const int beta_depth)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<SideData<NDIM, double> > src1_data = patch->getPatchData(src1_idx);
            Pointer<SideData<NDIM, double> > src2_data =
                (src2_idx >= 0) ? patch->getPatchData(src2_idx) : Pointer<PatchData<NDIM> >();
            Pointer<SideData<NDIM, double> > alpha_data = patch->getPatchData(alpha_idx);
            Pointer<SideData<NDIM, double> > beta_data = patch->getPatchData(beta_idx);

            d_patch_math_ops.pointwiseMultiply(dst_data, alpha_data, src1_data, beta_data, src2_data, patch, dst_depth,
                                               src1_depth, src2_depth, alpha_depth, beta_depth);
        }
    }
    return;
} // pointwiseMultiply

void HierarchyMathOps::pointwiseL1Norm(const int dst_idx,
                                       const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                       const int src_idx,
                                       const Pointer<CellVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseL1Norm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseL1Norm

void HierarchyMathOps::pointwiseL2Norm(const int dst_idx,
                                       const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                       const int src_idx,
                                       const Pointer<CellVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseL2Norm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseL2Norm

void HierarchyMathOps::pointwiseMaxNorm(const int dst_idx,
                                        const Pointer<CellVariable<NDIM, double> > /*dst_var*/,
                                        const int src_idx,
                                        const Pointer<CellVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<CellData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseMaxNorm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseMaxNorm

void HierarchyMathOps::pointwiseL1Norm(const int dst_idx,
                                       const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                       const int src_idx,
                                       const Pointer<NodeVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseL1Norm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseL1Norm

void HierarchyMathOps::pointwiseL2Norm(const int dst_idx,
                                       const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                       const int src_idx,
                                       const Pointer<NodeVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseL2Norm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseL2Norm

void HierarchyMathOps::pointwiseMaxNorm(const int dst_idx,
                                        const Pointer<NodeVariable<NDIM, double> > /*dst_var*/,
                                        const int src_idx,
                                        const Pointer<NodeVariable<NDIM, double> > /*src_var*/)
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<NodeData<NDIM, double> > dst_data = patch->getPatchData(dst_idx);
            Pointer<NodeData<NDIM, double> > src_data = patch->getPatchData(src_idx);

            d_patch_math_ops.pointwiseMaxNorm(dst_data, src_data, patch);
        }
    }
    return;
} // pointwiseMaxNorm

/////////////////////////////// PRIVATE //////////////////////////////////////

void HierarchyMathOps::resetCoarsenOperators()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_grid_geom);
#endif
    d_of_coarsen_op = d_grid_geom->lookupCoarsenOperator(d_of_var, d_coarsen_op_name);
    d_os_coarsen_op = d_grid_geom->lookupCoarsenOperator(d_os_var, d_coarsen_op_name);

    d_of_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_of_coarsen_alg->registerCoarsen(d_fc_idx, // destination
                                      d_of_idx, // source
                                      d_of_coarsen_op);

    d_os_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_os_coarsen_alg->registerCoarsen(d_sc_idx, // destination
                                      d_os_idx, // source
                                      d_os_coarsen_op);
    return;
} // resetCoarsenOperators

void HierarchyMathOps::resetRefineOperators()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_grid_geom);
#endif
    // intentionally blank
    return;
} // resetRefineOperators

void HierarchyMathOps::xeqScheduleOuterfaceRestriction(const int dst_idx, const int src_idx, const int dst_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(dst_ln >= d_coarsest_ln);
    TBOX_ASSERT(dst_ln + 1 <= d_finest_ln);
#endif
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_alg->registerCoarsen(dst_idx, src_idx, d_of_coarsen_op);
    if (coarsen_alg->checkConsistency(d_of_coarsen_scheds[dst_ln]))
    {
        coarsen_alg->resetSchedule(d_of_coarsen_scheds[dst_ln]);
        d_of_coarsen_scheds[dst_ln]->coarsenData();
        d_of_coarsen_alg->resetSchedule(d_of_coarsen_scheds[dst_ln]);
    }
    else
    {
        Pointer<PatchLevel<NDIM> > src_level = d_hierarchy->getPatchLevel(dst_ln + 1);
        Pointer<PatchLevel<NDIM> > dst_level = d_hierarchy->getPatchLevel(dst_ln);
        coarsen_alg->createSchedule(dst_level, src_level)->coarsenData();
    }
    return;
} // xeqScheduleOuterfaceRestriction

void HierarchyMathOps::xeqScheduleOutersideRestriction(const int dst_idx, const int src_idx, const int dst_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(dst_ln >= d_coarsest_ln);
    TBOX_ASSERT(dst_ln + 1 <= d_finest_ln);
#endif
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_alg->registerCoarsen(dst_idx, src_idx, d_os_coarsen_op);
    if (coarsen_alg->checkConsistency(d_os_coarsen_scheds[dst_ln]))
    {
        coarsen_alg->resetSchedule(d_os_coarsen_scheds[dst_ln]);
        d_os_coarsen_scheds[dst_ln]->coarsenData();
        d_os_coarsen_alg->resetSchedule(d_os_coarsen_scheds[dst_ln]);
    }
    else
    {
        Pointer<PatchLevel<NDIM> > src_level = d_hierarchy->getPatchLevel(dst_ln + 1);
        Pointer<PatchLevel<NDIM> > dst_level = d_hierarchy->getPatchLevel(dst_ln);
        coarsen_alg->createSchedule(dst_level, src_level)->coarsenData();
    }
    return;
} // xeqScheduleOutersideRestriction

////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
