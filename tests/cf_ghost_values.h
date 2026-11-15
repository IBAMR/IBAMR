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

// Comparison of the coarse-fine ghost values of two sets of patch data on the same level.

#ifndef included_ibamr_tests_cf_ghost_values_h
#define included_ibamr_tests_cf_ghost_values_h

#include <ibamr/config.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideIterator.h>

#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// A smooth function of position with integer wave numbers, so that it is periodic on the unit square or cube.
inline double
cf_ghost_test_function(const IBTK::VectorNd& x, const double salt)
{
    double arg_a = 0.0, arg_b = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        arg_a += (d == 0 ? 2.0 : 1.0) * x[d];
        arg_b += (d == 1 ? -3.0 : 1.0) * x[d];
    }
    return std::sin(2.0 * M_PI * arg_a + salt) + 0.5 * std::cos(2.0 * M_PI * arg_b - salt);
}

// Set the interior values of the cell-centered or side-centered data with index data_idx on the levels 0 and 1 to the
// test function, and the ghost values to a sentinel.
inline void
set_cf_ghost_test_data(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                       const int data_idx,
                       const double salt)
{
    using namespace SAMRAI;

    for (int ln = 0; ln <= 1; ++ln)
    {
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            tbox::Pointer<geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            tbox::Pointer<pdat::CellData<NDIM, double>> cc_data = patch->getPatchData(data_idx);
            tbox::Pointer<pdat::SideData<NDIM, double>> sc_data = patch->getPatchData(data_idx);
            if (cc_data)
            {
                cc_data->fillAll(-12345.678);
                for (pdat::CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    IBTK::VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(ci()(d) - idx_low(d)) + 0.5);
                    for (int depth = 0; depth < cc_data->getDepth(); ++depth)
                    {
                        (*cc_data)(ci(), depth) = cf_ghost_test_function(x, salt + 0.5 * depth);
                    }
                }
            }
            else
            {
                sc_data->fillAll(-12345.678);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (pdat::SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        IBTK::VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] =
                                xlow[d] + dx[d] * (static_cast<double>(si()(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                        for (int depth = 0; depth < sc_data->getDepth(); ++depth)
                        {
                            (*sc_data)(si(), depth) = cf_ghost_test_function(x, salt + axis + 0.5 * depth);
                        }
                    }
                }
            }
        }
    }
}

// Compare bit for bit the values of the data with indices u_idx and w_idx on the cells (faces) next to the boxes of
// the patches of level ln, without corners, that the level does not cover and that lie in the physical domain, at
// least two cells of the next coarser level away from a physical boundary, where the two fills differ (a known
// defect). Return the number of entries compared on all processes. Stop with an error at the first difference.
inline int
compare_cf_ghost_values(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                        const int ln,
                        const int u_idx,
                        const int w_idx,
                        const std::string& name)
{
    using namespace SAMRAI;

    tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
    const hier::BoxArray<NDIM>& level_boxes = level->getBoxes();
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    tbox::Pointer<hier::PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(ln - 1);
    const hier::Box<NDIM> coarse_domain = coarse_level->getPhysicalDomain()[0];
    const hier::IntVector<NDIM> ratio_to_coarser = level->getRatioToCoarserLevel();
    const hier::IntVector<NDIM> shift = hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio());

    // Cell data have comp -1. The faces normal to comp are indexed by the cell above the face.
    auto index_box = [](hier::Box<NDIM> box, const int comp)
    {
        if (comp >= 0) box.upper(comp) += 1;
        return box;
    };

    // The images of an index along the periodic directions.
    std::vector<hier::IntVector<NDIM>> offsets;
    {
        int num = 1;
        for (int d = 0; d < NDIM; ++d) num *= 3;
        for (int k = 0; k < num; ++k)
        {
            hier::IntVector<NDIM> offset(0);
            int rem = k;
            bool allowed = true;
            for (int d = 0; d < NDIM; ++d)
            {
                const int m = rem % 3 - 1;
                rem /= 3;
                if (m != 0 && shift(d) == 0) allowed = false;
                offset(d) = m * shift(d);
            }
            if (allowed) offsets.push_back(offset);
        }
    }

    int count = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> u_cc = patch->getPatchData(u_idx);
        tbox::Pointer<pdat::CellData<NDIM, double>> w_cc = patch->getPatchData(w_idx);
        tbox::Pointer<pdat::SideData<NDIM, double>> u_sc = patch->getPatchData(u_idx);
        tbox::Pointer<pdat::SideData<NDIM, double>> w_sc = patch->getPatchData(w_idx);
        const int first_comp = u_cc ? -1 : 0;
        const int last_comp = u_cc ? -1 : NDIM - 1;
        for (int comp = first_comp; comp <= last_comp; ++comp)
        {
            const hier::Box<NDIM> box = index_box(patch->getBox(), comp);
            const hier::Box<NDIM> domain_box = index_box(domain, comp);
            for (hier::Box<NDIM>::Iterator b(hier::Box<NDIM>::grow(box, 1)); b; b++)
            {
                const hier::Index<NDIM>& i = b();
                // Keep the entries next to the box along exactly one direction.
                int outside = 0;
                bool keep = true;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (i(d) < box.lower(d) || i(d) > box.upper(d)) ++outside;
                    if (shift(d) == 0 && (i(d) < domain_box.lower(d) || i(d) > domain_box.upper(d))) keep = false;
                }
                // The coarse cell that contains the entry.
                for (int d = 0; d < NDIM && keep; ++d)
                {
                    const int r = ratio_to_coarser(d);
                    const int c = i(d) >= 0 ? i(d) / r : -((-i(d) + r - 1) / r);
                    if (shift(d) == 0 && (c < coarse_domain.lower(d) + 2 || c > coarse_domain.upper(d) - 2))
                    {
                        keep = false;
                    }
                }
                if (outside != 1 || !keep) continue;
                bool covered = false;
                for (const auto& offset : offsets)
                {
                    hier::Index<NDIM> image = i;
                    for (int d = 0; d < NDIM; ++d) image(d) += offset(d);
                    for (int k = 0; k < level_boxes.getNumberOfBoxes() && !covered; ++k)
                    {
                        covered = index_box(level_boxes[k], comp).contains(image);
                    }
                }
                if (covered) continue;
                const int depth_count = u_cc ? u_cc->getDepth() : u_sc->getDepth();
                for (int depth = 0; depth < depth_count; ++depth)
                {
                    const double u = u_cc ?
                                         (*u_cc)(pdat::CellIndex<NDIM>(i), depth) :
                                         (*u_sc)(pdat::SideIndex<NDIM>(i, comp, pdat::SideIndex<NDIM>::Lower), depth);
                    const double w = u_cc ?
                                         (*w_cc)(pdat::CellIndex<NDIM>(i), depth) :
                                         (*w_sc)(pdat::SideIndex<NDIM>(i, comp, pdat::SideIndex<NDIM>::Lower), depth);
                    ++count;
                    if (std::memcmp(&u, &w, sizeof(double)) != 0)
                    {
                        std::ostringstream os;
                        os << std::setprecision(17) << "index (" << i(0);
                        for (int d = 1; d < NDIM; ++d) os << "," << i(d);
                        os << ") depth " << depth << " of patch " << p() << " has " << u << " and " << w;
                        TBOX_ERROR("compare_cf_ghost_values():\n"
                                   << "  " << name << ": coarse-fine ghost value of "
                                   << (u_cc ? std::string("cell data") : "side data of axis " + std::to_string(comp))
                                   << " differs at " << os.str() << "\n");
                    }
                }
            }
        }
    }
    return IBTK::IBTK_MPI::sumReduction(count);
}

#endif // #ifndef included_ibamr_tests_cf_ghost_values_h
