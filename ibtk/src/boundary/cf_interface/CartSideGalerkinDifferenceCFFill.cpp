#include "ibtk/CartSideLinearGalerkinDifferenceCFFill.h"

#include "BoundaryBox.h"
#include "CoarseFineBoundary.h"
#include "ComponentSelector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"             
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>
    
#include "ibtk/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES !!! TODO
#if (NDIM == 2)
#define SC_LIN_GALERKIN_DIFFERENCE_FILL_FC IBTK_FC_FUNC(sclingalerkindifferencefill2d, SCLINGALERKDIFFFILL2D)
#endif
#if (NDIM == 3)
#define SC_LIN_GALERKIN_DIFFERENCE_FILL_FC IBTK_FC_FUNC(sclingalerkindifferencefill3d, SCLINGALERKINDIFFFERENCEILL3D)
#endif
// fortran function interfaces
extern "C"
{
}

namespace IBTK
{

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_TILL = 1;
} // namespace


// public
CartSideLinearGalerkinDifferenceCFFill::CartSideLinearGalerkinDifferenceCFFill()
{
	// set up scratch
	VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
	Pointer<VariableContext> context = var_db->getContext("CartSideLinearGalerkinDifferenceCFFill::CONTEXT");
	if (var_db->checkVariableExists(d_sc_indicator_var->getName());
	{
		d_sc_indicator_var = var_db->getVariable(d_sc_indicator_var->getName());
		d_sc_indicator_idx = var_db->mapVariableAndContextToIndex(d_sc_indicator_var, context);
	}
	else
	{
		d_sc_indicator_indx = var_db_>registerVariableAndContext(d_sc_indicator_var, context, GHOST_WIDTH_TO_FILL);
	}
	return;
} // 
	
CartSideLinearGalerkinDifferenceCFFill::~CartSideLinearGalerkinDifferenceCFFill()
{
    clearPatchHierarchy();
    return;
} // ~CartSideLinearGalerkinDifferenceCFFill

void
CartSideLinearGalerkinDifferenceCFFill::setPhysicalBoundaryConditions(Patch<NDIM>& /*patch*/,
                                                                      const double /*fill_time*/,
                                                                      const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideLinearGalerkinDifferenceCFFill::getRefineOpStencilWidth() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_refine_op->getStencilWidth().max() <= REFINE_OP_STENCIL_WIDTH);
#endif
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartSideLinearGalerkinDifferenceCFFill::preprocessRefine(Patch<NDIM>& /*fine*/,
                                                         const Patch<NDIM>& /*coarse*/,
                                                         const Box<NDIM>& /*fine_box*/,
                                                         const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartSideLinearGalerkinDifferenceCFFill::postprocessRefine(Patch<NDIM>& fine,
                                                          const Patch<NDIM>& coarse,
                                                          const Box<NDIM>& fine_box,
                                                          const IntVector<NDIM>& ratio)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!fine.inHierarchy())
    {
        for (const auto& patch_data_index : d_patch_data_indices)
        {
            // this should not be used hopefully
            d_refine_op->refine(fine, coarse, patch_data_index, patch_data_index, fine_box, ratio);
        }
        return;
    }
#if !defined(NDEBUG)
    else
    {
        // Ensure the fine patch corresponds to the expected patch in the cached
        // patch hierarchy.
        const int patch_num = fine.getPatchNumber();
        const int fine_patch_level_num = fine.getPatchLevelNumber();
        Pointer<PatchLevel<NDIM> > fine_level = d_hierarchy->getPatchLevel(fine_patch_level_num);
        TBOX_ASSERT(&fine == fine_level->getPatch(patch_num).getPointer());
    }
#endif
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = fine.getPatchNumber();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM> >& cf_bdry_codim1_boxes =
        d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 1);
    if (cf_bdry_codim1_boxes.size() == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        Pointer<SideData<NDIM, double> > fdata = fine.getPatchData(patch_data_index);
        Pointer<SideData<NDIM, double> > cdata = coarse.getPatchData(patch_data_index);
        Pointer<SideData<NDIM, int> > indicator_data = fine.getPatchData(d_sc_indicator_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
        TBOX_ASSERT(indicator_data);
#endif
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideLinearGalerkinDifferenceCFFill::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideLinearGalerkinDifferenceCFFill::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
#endif
        const int data_depth = fdata->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM> > pgeom_fine = fine.getPatchGeometry();
        const Box<NDIM>& patch_box_fine = fine.getBox();
        const Box<NDIM>& patch_box_crse = coarse.getBox();
        for (int k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif

            for (int depth = 0; depth < data_depth; ++depth)
            {
                for (int dim = 0; dim < NDIM; dim++){
                    double* const U_fine0 = fdata->getPointer(dim, depth);
                    const double* const U_crse0 = cdata->getPointer(dim, depth);
                    int ibeg = max(bc_fill_box.lower(dim),patch_box_fine.lower(dim));
                    int iend = min(bc_fill_box.upper(dim),patch_box_fine.upper(dim));
                    int jbeg = max(bc_fill_box.lower(1-dim),patch_box_fine.lower(1-dim));
                    int jend = min(bc_fill_box.upper(1-dim),patch_box_fine.upper(1-dim));
                    // let's review what everything is
                    // U_fine012 fine side 012 data
                    // U_fine_ghosts seems like integer ghost cell width
                    // U_crse012 coarse side 012 data
                    // u coarse ghosts seems like integer with coarse ghost cell width
                    // indicator012 pointer to fine patch data 012 using patch indicator index
                    // indicator ghosts uses sc indicator index not patch data
                    // patch box indices for fine and coarse
                    // location index indicates the location of the boundary box in relation to the location of the associated patch
                    // ratio should be refinement ratio
                    // boundary fill box is from patch geometry

                    //TODO need to add a check for corner case -- want it handled by aligned axis refinement
                    // i.e. for u we want it on li = 0/1 and for v we want li = 2/3

                    if (location_index==0 || location_index==1){
                        for (int j=jbeg;j<=jend; j++){
                            if (j < 0)
                                int j_c = (j+1)/ratio(1-dim) - 1;
                            else
                                int j_c = j/ratio(1-dim);
                            for (int i=bc_fill_box.lower(dim),i<=bc_fill_box.upper(dim),i++){
                                if (i<0)
                                    int i_c = (i+1)/ratio(dim)-1;
                                else
                                    int i_c = i/ratio(dim);
                                int j_wgt = (j%2) ? 1 : -1;
                                if (dim==0)
                                    U_fine(i,j) = .75*U_crse(i_c,j_c)+.25*U_crse(i_c,j_c+j_wgt);
                                else
                                    if (i%2)
                                        U_fine(j,i) = .5*U_crse(j_c,i_c)+.5*U_crse(j_c+j_wgt,i_c) +
                                                   .5*U_crse(j_c,i_c+1)+.5*U_crse(j_c+j_wgt,i_c+1);
                                    else
                                        U_fine(j,i) = U_crse(j_c,i_c)+U_crse(j_c+j_wgt,i_c);
                            }
                        }
                        // set upper and lower x side of patch
                    }
                    else if (location_index==2 || location_index==3){
                        // set upper and lower y side of patch 
                        for (int j=bc_fill_box.lower(1-dim),j<=bc_fill_box.upper(1-dim),j++){
                            if (j < 0)
                                int j_c = (j+1)/ratio(1-dim) - 1;
                            else
                                int j_c = j/ratio(1-dim);
                            for (int i=ibeg;i<=iend; i++){
                                if (i<0)
                                    int i_c = (i+1)/ratio(dim)-1;
                                else
                                    int i_c = i/ratio(dim);
                                int j_wgt = (j%2) ? 1 : -1;

                                if (dim==0)
                                    if (i%2)
                                        U_fine(i,j) = .5*U_crse(i_c,j_c)+.5*U_crse(i_c,j_c+j_wgt) +
                                                   .5*U_crse(i_c+1,j_c)+.5*U_crse(i_c+1,j_c+j_wgt);
                                    else
                                        U_fine(i,j) = U_crse(i_c,j_c)+U_crse(i_c,j_c+j_wgt);
                                else
                                    U_fine(j,i) = .75*U_crse(j_c,i_c)+.25*U_crse(j_c+j_wgt,i_c);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
} // postprocessRefine

void
CartSideLinearGalerkinDifferenceCFFill::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartSideLinearGalerkinDifferenceCFFill::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartSideLinearGalerkinDifferenceCFFill::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartSideLinearGalerkinDifferenceCFFill::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
} // setPatchDataIndices

void
CartSideLinearGalerkinDifferenceCFFill::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();

    d_cf_boundary.resize(finest_level_number + 1);
    const IntVector<NDIM>& max_ghost_width = GHOST_WIDTH_TO_FILL;
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = CoarseFineBoundary<NDIM>(*d_hierarchy, ln, max_ghost_width);
    }

    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = nullptr;
    refine_alg->registerRefine(d_sc_indicator_idx, // destination
                               d_sc_indicator_idx, // source
                               d_sc_indicator_idx, // temporary work space
                               refine_op);
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_sc_indicator_idx))
        {
            level->allocatePatchData(d_sc_indicator_idx, 0.0);
        }
        else
        {
            level->setTime(0.0, d_sc_indicator_idx);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, int> > sc_indicator_data = patch->getPatchData(d_sc_indicator_idx);
            sc_indicator_data->fillAll(0, sc_indicator_data->getGhostBox());
            sc_indicator_data->fillAll(1, sc_indicator_data->getBox());
        }
        refine_alg->createSchedule(d_hierarchy->getPatchLevel(ln))->fillData(0.0);
    }
	return;
) // set patch hierarchy

void
CartSideLinearGalerkinDifferenceCFFill::clearPatchHierarchy()
{
    d_hierarchy.setNull();
    d_cf_boundary.clear();
    return;
} // clearPatchHierarchy

oid
CartSideLinearGalerkinDifferenceCFFill::computeNormalExtension(Patch<NDIM>& patch,
                                                               const IntVector<NDIM>& ratio,
                                                               const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy())
    {
        return;
    }
#if !defined(NDEBUG)
    else
    {
        const int patch_num = patch.getPatchNumber();
        const int patch_level_num = patch.getPatchLevelNumber();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(patch_level_num);
        TBOX_ASSERT(&patch == level->getPatch(patch_num).getPointer());
    }
#endif
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = patch.getPatchNumber();
    const int patch_level_num = patch.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM> >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        Pointer<SideData<NDIM, double> > data = patch.getPatchData(patch_data_index);
        SideData<NDIM, double> data_copy(data->getBox(), data->getDepth(), data->getGhostCellWidth());
        data_copy.copyOnBox(*data, data->getGhostBox());
        Pointer<SideData<NDIM, int> > indicator_data = patch.getPatchData(d_sc_indicator_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
        TBOX_ASSERT(indicator_data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
        const int W_ghosts = (data_copy.getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideLinearGalerkinDifferenceCFFill::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (W_ghosts != (data_copy.getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideLinearGalerkinDifferenceCFFill::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
#endif
        const int data_depth = data->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
        const Box<NDIM>& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif
            for (int depth = 0; depth < data_depth; ++depth)
            {
                for (int dim = 0; dim < NDIM; dim++){
                    double* const U = data->getPointer(dim, depth);
                    const double* const W = data_copy->getPointer(dim, depth);
                    int ibeg = max(bc_fill_box.lower(dim),patch_box_fine.lower(dim));
                    int iend = min(bc_fill_box.upper(dim),patch_box_fine.upper(dim));
                    int jbeg = max(bc_fill_box.lower(1-dim),patch_box_fine.lower(1-dim));
                    int jend = min(bc_fill_box.upper(1-dim),patch_box_fine.upper(1-dim));

                    if (dim==1 & (location_index==0 || location_index==1)){
                        for (int j=jbeg;j<=jend; j++){
                            for (int i=bc_fill_box.lower(dim),i<=bc_fill_box.upper(dim),i++){
                                int j_wgt = (j%2) ? 1 : -1;
                                U(j,i) = W(j,i) - W(j+j_wgt,i);
                            }
                        }
                        // set upper and lower x side of patch
                    }
                    else if (dim==0 & (location_index==2 || location_index==3)){
                        // set upper and lower y side of patch 
                        for (int j=bc_fill_box.lower(1-dim),j<=bc_fill_box.upper(1-dim),j++){
                            for (int i=ibeg;i<=iend; i++){
                                int j_wgt = (j%2) ? 1 : -1;
                                U(i,j) = W(i,j) - W(i,j+j_wgt);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
} // computeNormalExtension


} //namespace IBTK
	
