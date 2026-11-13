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

// Set-up and checks of the ghost values on the finest level of a hierarchy, shared by the cf_interface and
// interface_face_ghost_copy tests.

#ifndef included_ibamr_coarse_fine_ghost_values_h
#define included_ibamr_coarse_fine_ghost_values_h

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <GriddingAlgorithm.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideIndex.h>
#include <SideIterator.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace coarse_fine_ghost_values
{
// Value of the ghost cells before the fill, so that cells that nothing sets are equal on all patches.
constexpr double SENTINEL = -12345.678;

// The fine level: the union F of its boxes, and the boxes of the physical domain, whose union need not be a box,
// in the fine index space (domain) and in that of level 0 (coarse_domain).
struct FineLevel
{
    BoxArray<NDIM> boxes;
    BoxArray<NDIM> domain;
    BoxArray<NDIM> coarse_domain;
    IntVector<NDIM> ratio;            // To level 0.
    IntVector<NDIM> ratio_to_coarser; // To the next coarser level.
};

// Cell data have axis -1. Side data have the faces normal to axis (axis >= 0), indexed by the cell above the face.
inline Box<NDIM>
index_box(const Box<NDIM>& cell_box, const int axis)
{
    Box<NDIM> box = cell_box;
    if (axis >= 0)
    {
        box.upper(axis) += 1;
    }
    return box;
}

// Whether one of the cell boxes contains the cell or the face i.
inline bool
in_boxes(const BoxArray<NDIM>& boxes, const int axis, const hier::Index<NDIM>& i)
{
    for (int k = 0; k < boxes.getNumberOfBoxes(); ++k)
    {
        if (index_box(boxes[k], axis).contains(i))
        {
            return true;
        }
    }
    return false;
}

inline bool
in_fine_level(const FineLevel& fine, const int axis, const hier::Index<NDIM>& i)
{
    return in_boxes(fine.boxes, axis, i);
}

inline int
floor_div(const int i, const int r)
{
    return i >= 0 ? i / r : -((-i + r - 1) / r);
}

// Cells (and faces) whose coarse stencil is far from the physical boundary, where the test does not set boundary
// values: the cells of level 0 within two cells of the one that contains them lie in the domain.
inline bool
away_from_boundary(const FineLevel& fine, const hier::Index<NDIM>& i)
{
    hier::Index<NDIM> c;
    for (int d = 0; d < NDIM; ++d)
    {
        c(d) = floor_div(i(d), fine.ratio(d));
    }
    BoxList<NDIM> outside(Box<NDIM>::grow(Box<NDIM>(c, c), IntVector<NDIM>(2)));
    outside.removeIntersections(BoxList<NDIM>(fine.coarse_domain));
    return outside.isEmpty();
}

// The directions (d, s) along which a ghost cell (or face) i gets a coarse-fine value: i is not in F, the first cell
// (face) of F along the direction is at the distance k, at most the depth limit, and the cell (face) after it is also
// in F.
inline std::vector<std::pair<int, int>>
coarse_fine_terms(const FineLevel& fine, const int axis, const hier::Index<NDIM>& i, const int depth_limit)
{
    std::vector<std::pair<int, int>> terms;
    if (!in_boxes(fine.domain, axis, i) || in_fine_level(fine, axis, i))
    {
        return terms;
    }
    for (int d = 0; d < NDIM; ++d)
    {
        for (const int s : { -1, +1 })
        {
            const int limit = std::min(depth_limit, fine.ratio_to_coarser(d));
            for (int k = 1; k <= limit; ++k)
            {
                hier::Index<NDIM> j1 = i, j2 = i;
                j1(d) += k * s;
                j2(d) += (k + 1) * s;
                if (!in_fine_level(fine, axis, j1))
                {
                    continue;
                }
                if (in_fine_level(fine, axis, j2))
                {
                    terms.emplace_back(d, s);
                }
                break;
            }
        }
    }
    return terms;
}

// A ghost value of one patch.
struct GhostValue
{
    int axis;
    int idx[NDIM];
    int first_layer; // In the first ghost layer of this patch, within one cell of its tangential range.
    double value;
};

inline bool
less_key(const GhostValue& a, const GhostValue& b)
{
    if (a.axis != b.axis)
    {
        return a.axis < b.axis;
    }
    for (int d = NDIM - 1; d >= 0; --d)
    {
        if (a.idx[d] != b.idx[d])
        {
            return a.idx[d] < b.idx[d];
        }
    }
    return false;
}

inline bool
same_key(const GhostValue& a, const GhostValue& b)
{
    return !less_key(a, b) && !less_key(b, a);
}

inline bool
same_bits(const double a, const double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

inline std::string
value_string(const double v)
{
    std::ostringstream os;
    os.precision(17);
    os << v;
    return os.str();
}

inline std::string
key_string(const GhostValue& v)
{
    std::string s = "axis " + std::to_string(v.axis) + " index (";
    for (int d = 0; d < NDIM; ++d)
    {
        s += (d ? "," : "") + std::to_string(v.idx[d]);
    }
    return s + ")";
}

// Value of a patch datum at index i, for cell data (axis = -1) or side data.
inline double
get_value(Pointer<Patch<NDIM>> patch, const int data_idx, const int axis, const hier::Index<NDIM>& i)
{
    if (axis < 0)
    {
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
        return (*data)(CellIndex<NDIM>(i));
    }
    Pointer<SideData<NDIM, double>> data = patch->getPatchData(data_idx);
    return (*data)(SideIndex<NDIM>(i, axis, SideIndex<NDIM>::Lower));
}

// Box of the cell indices (axis = -1) or the side indices in the ghost box of the patch data.
inline Box<NDIM>
ghost_index_box(Pointer<Patch<NDIM>> patch, const int data_idx, const int axis)
{
    Pointer<PatchData<NDIM>> data = patch->getPatchData(data_idx);
    return index_box(data->getGhostBox(), axis);
}

template <typename F>
inline void
for_each_index(const Box<NDIM>& box, const int axis, F f)
{
    for (BoxIterator<NDIM> bi(box); bi; bi++)
    {
        f(axis, bi());
    }
}

// The values of all processes, in the order of the ranks.
template <typename T>
std::vector<T>
gather_all(const std::vector<T>& mine)
{
    const int nbytes = static_cast<int>(mine.size() * sizeof(T));
    const int total = IBTK_MPI::sumReduction(nbytes);
    std::vector<T> all(total / sizeof(T));
    IBTK_MPI::allGather(reinterpret_cast<const char*>(mine.data()), nbytes, reinterpret_cast<char*>(all.data()), total);
    return all;
}

// The axes of the data: -1 for cell data, and 0 to NDIM - 1 for side data.
inline std::vector<int>
data_axes(const bool is_cell)
{
    std::vector<int> axes;
    if (is_cell)
    {
        axes.push_back(-1);
    }
    else
    {
        for (int d = 0; d < NDIM; ++d)
        {
            axes.push_back(d);
        }
    }
    return axes;
}

// The patch hierarchy that the input describes: the levels in FINE_BOXES and FINE_BOXES_2, or the levels that the
// gridding algorithm makes. Register the variables first, to set the maximum ghost width.
inline Pointer<PatchHierarchy<NDIM>>
make_hierarchy(Pointer<AppInitializer> app_initializer, Pointer<CartesianGridGeometry<NDIM>> grid_geometry)
{
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);
    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    if (input_db->keyExists("FINE_BOXES"))
    {
        // The boxes of the fine levels, spread over the processes. The refinement ratios are one integer for every
        // axis or one integer for each axis.
        auto read_ratio = [&](const std::string& key)
        {
            const tbox::Array<int> values = input_db->getIntegerArray(key);
            IntVector<NDIM> ratio;
            for (int d = 0; d < NDIM; ++d)
            {
                ratio(d) = values[values.size() == 1 ? 0 : d];
            }
            return ratio;
        };
        IntVector<NDIM> ratio_to_level_0(1);
        for (int ln = 1; ln <= 2; ++ln)
        {
            const std::string suffix = ln == 1 ? "" : "_2";
            if (!input_db->keyExists("FINE_BOXES" + suffix))
            {
                break;
            }
            const tbox::Array<tbox::DatabaseBox> db_boxes = input_db->getDatabaseBoxArray("FINE_BOXES" + suffix);
            BoxArray<NDIM> fine_boxes(db_boxes.size());
            ProcessorMapping mapping(db_boxes.size());
            for (int k = 0; k < db_boxes.size(); ++k)
            {
                fine_boxes[k] = Box<NDIM>(db_boxes[k]);
                mapping.setProcessorAssignment(k, k % IBTK_MPI::getNodes());
            }
            ratio_to_level_0 *= read_ratio("REF_RATIO" + suffix);
            patch_hierarchy->makeNewPatchLevel(ln, ratio_to_level_0, fine_boxes, mapping);
        }
    }
    else
    {
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }
    }
    return patch_hierarchy;
}

// The quadratic function of the position. Along a periodic direction, it is a quadratic of the distance in the period,
// so that the faces on the two sides of a periodic boundary have the same value.
inline std::function<double(const VectorNd&)>
make_quadratic(Pointer<CartesianGridGeometry<NDIM>> grid_geometry)
{
    const IntVector<NDIM> periodic_shift = grid_geometry->getPeriodicShift(IntVector<NDIM>(1));
    return [=](const VectorNd& x) -> double
    {
        double val = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift(d) != 0)
            {
                const double s = (x[d] - grid_geometry->getXLower()[d]) /
                                 (grid_geometry->getXUpper()[d] - grid_geometry->getXLower()[d]);
                val += s * (1.0 - s);
            }
            else
            {
                val += x[d] + x[d] * x[d];
            }
        }
        return val;
    };
}

// The exact value of the quadratic function at a cell or a face of a patch.
inline double
exact_value(Pointer<Patch<NDIM>> patch,
            const int axis,
            const hier::Index<NDIM>& i,
            const std::function<double(const VectorNd&)>& Q_fcn)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const xlow = pgeom->getXLower();
    const hier::Index<NDIM>& idx_low = patch->getBox().lower();
    VectorNd x;
    for (int d = 0; d < NDIM; ++d)
    {
        x[d] = xlow[d] + dx[d] * (static_cast<double>(i(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
    }
    return Q_fcn(x);
}

// Set the interiors of the patch data to the quadratic function and the ghost values to the sentinel.
inline void
fill_quadratic(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
               const std::vector<int>& data_idxs,
               const std::function<double(const VectorNd&)>& Q_fcn)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            for (const int data_idx : data_idxs)
            {
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(data_idx);
                Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(data_idx);
                if (Q_cc_data)
                {
                    Q_cc_data->fillAll(SENTINEL);
                    for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        }
                        (*Q_cc_data)(idx) = Q_fcn(x);
                    }
                }
                else if (Q_sc_data)
                {
                    Q_sc_data->fillAll(SENTINEL);
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                x[d] = xlow[d] +
                                       dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            }
                            (*Q_sc_data)(idx) = Q_fcn(x);
                        }
                    }
                }
            }
        }
    }
}

// The finest level of the hierarchy as a FineLevel.
inline FineLevel
make_fine_level(Pointer<PatchHierarchy<NDIM>> patch_hierarchy, Pointer<CartesianGridGeometry<NDIM>> grid_geometry)
{
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
    FineLevel fine;
    fine.boxes = level->getBoxes();
    fine.ratio = level->getRatio();
    fine.ratio_to_coarser = level->getRatioToCoarserLevel();
    fine.domain = grid_geometry->getPhysicalDomain();
    fine.domain.refine(level->getRatio());
    fine.coarse_domain = grid_geometry->getPhysicalDomain();
    return fine;
}

// Check that every patch that contains a cell or a face of the fine level, as an interior value or as a ghost
// value, has the same value there, to the bit. A face on the boundary of a patch along a coarse-fine interface
// belongs to the fine level, so the ghost copies of it in other patches, on this or another process, must equal
// the value of the patch that contains it, also for patches that do not touch that patch. Cells and faces at
// periodic images are the same cell or face. Print the result and add a failure to failures.
inline void
check_fine_level_agreement(Pointer<PatchLevel<NDIM>> level,
                           const FineLevel& fine,
                           Pointer<CartesianGridGeometry<NDIM>> grid_geometry,
                           const int Q_idx,
                           const std::vector<int>& axes,
                           const bool is_cell,
                           std::ostringstream& failures)
{
    const IntVector<NDIM> shift = grid_geometry->getPeriodicShift(level->getRatio());
    std::vector<GhostValue> mine;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        for (const int axis : axes)
        {
            for_each_index(ghost_index_box(patch, Q_idx, axis),
                           axis,
                           [&](const int ax, const hier::Index<NDIM>& i)
                           {
                               hier::Index<NDIM> j = i;
                               for (int d = 0; d < NDIM; ++d)
                               {
                                   if (shift(d) == 0)
                                   {
                                       continue;
                                   }
                                   const int lo = fine.domain.getBoundingBox().lower(d);
                                   j(d) = lo + ((i(d) - lo) % shift(d) + shift(d)) % shift(d);
                               }
                               if (!in_fine_level(fine, ax, j))
                               {
                                   return;
                               }
                               GhostValue v;
                               v.axis = ax;
                               for (int d = 0; d < NDIM; ++d)
                               {
                                   v.idx[d] = j(d);
                               }
                               v.first_layer = index_box(patch->getBox(), ax).contains(i) ? 0 : 1;
                               v.value = get_value(patch, Q_idx, ax, i);
                               mine.push_back(v);
                           });
        }
    }
    std::vector<GhostValue> all = gather_all(mine);
    std::stable_sort(all.begin(), all.end(), less_key);

    int num_shared = 0, num_differ = 0;
    std::string examples;
    for (std::size_t first = 0; first < all.size();)
    {
        std::size_t last = first;
        while (last < all.size() && same_key(all[first], all[last]))
        {
            ++last;
        }
        if (last - first > 1)
        {
            ++num_shared;
            bool differ = false;
            for (std::size_t k = first + 1; k < last; ++k)
            {
                if (!same_bits(all[first].value, all[k].value))
                {
                    differ = true;
                }
            }
            if (differ)
            {
                if (num_differ < 40)
                {
                    examples += "    " + key_string(all[first]) + ":";
                    for (std::size_t k = first; k < last; ++k)
                    {
                        examples += " " + value_string(all[k].value) + (all[k].first_layer ? "*" : "");
                    }
                    examples += "\n";
                }
                ++num_differ;
            }
        }
        first = last;
    }
    pout << "Compared " << num_shared << (is_cell ? " cells" : " faces")
         << " of the fine level that several patches contain: " << num_differ << " differ.\n";
    if (num_differ > 0)
    {
        failures << "  patches disagree on " << num_differ << (is_cell ? " cells" : " faces")
                 << " of the fine level (* marks a ghost copy):\n"
                 << examples;
    }
}

// Check that the fill does not change the values in the interiors of the patches, including the faces of their
// boundaries, to the bit. Print the result and add a failure to failures.
inline void
check_interior_unchanged(Pointer<PatchLevel<NDIM>> level,
                         const std::vector<int>& data_idxs,
                         const std::vector<int>& axes,
                         const std::function<double(const VectorNd&)>& Q_fcn,
                         std::ostringstream& failures)
{
    int num_compared = 0, num_changed = 0;
    std::string first_changed;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        for (const int data_idx : data_idxs)
        {
            for (const int axis : axes)
            {
                for_each_index(index_box(patch->getBox(), axis),
                               axis,
                               [&](const int ax, const hier::Index<NDIM>& i)
                               {
                                   ++num_compared;
                                   const double value = get_value(patch, data_idx, ax, i);
                                   if (!same_bits(value, exact_value(patch, ax, i, Q_fcn)) && num_changed++ == 0)
                                   {
                                       std::ostringstream os;
                                       os << "axis " << ax << " index " << i << " patch data " << data_idx << " value "
                                          << value_string(value);
                                       first_changed = os.str();
                                   }
                               });
            }
        }
    }
    IBTK_MPI::sumReduction(&num_compared, 1);
    int num_changed_all = num_changed;
    IBTK_MPI::sumReduction(&num_changed_all, 1);
    pout << "Compared " << num_compared
         << " interior values of the patches with their values before the fill: " << num_changed_all << " changed.\n";
    if (num_changed_all > 0)
    {
        failures << "  " << num_changed_all
                 << " interior values changed by the fill; first on process with a change: " << first_changed << "\n";
    }
}
} // namespace coarse_fine_ghost_values

#endif // included_ibamr_coarse_fine_ghost_values_h
