// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/MarkerPatchHierarchy.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/BoxArrayIndex.h>

#include <tbox/RestartManager.h>

#include <BoxArray.h>
#include <CartesianGridGeometry.h>

// At this point SAMRAI has already defined these macros, but older versions
// hdf5.h will unconditionally define them too. Get around it by undefining and
// redefining as needed.

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif

#include <hdf5.h>

#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <utility>

#include <ibtk/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace
{
static Timer* t_reinit;
static Timer* t_collect_all_markers;
static Timer* t_forward_euler_step;
static Timer* t_midpoint_step;
static Timer* t_trapezoidal_step;
static Timer* t_prune_and_redistribute;

std::vector<std::vector<hier::Box<NDIM>>>
compute_nonoverlapping_patch_boxes(const Pointer<BasePatchLevel<NDIM>>& c_level,
                                   const Pointer<BasePatchLevel<NDIM>>& f_level)
{
    const Pointer<PatchLevel<NDIM>> coarse_level = c_level;
    const Pointer<PatchLevel<NDIM>> fine_level = f_level;
    TBOX_ASSERT(coarse_level);
    TBOX_ASSERT(fine_level);
    TBOX_ASSERT(coarse_level->getLevelNumber() + 1 == fine_level->getLevelNumber());

    const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();

    // Get all (including those not on this processor) fine-level boxes:
    BoxList<NDIM> finer_box_list;
    long combined_size = 0;
    for (int i = 0; i < fine_level->getNumberOfPatches(); ++i)
    {
        Box<NDIM> patch_box = fine_level->getBoxForPatch(i);
        patch_box.coarsen(ratio);
        combined_size += patch_box.size();
        finer_box_list.addItem(patch_box);
    }
    finer_box_list.simplifyBoxes();

    // Remove said boxes from each coarse-level patch, including those that are not on this process. Removing a box that
    // does not intersect the patch changes nothing, so the loop removes from each patch only the boxes that intersect
    // it, in the order of the list:
    const detail::BoxArrayIndex finer_boxes(BoxArray<NDIM>(finer_box_list),
                                            std::vector<IntVector<NDIM>>(1, IntVector<NDIM>(0)));
    std::vector<int> finer_box_indices;
    std::vector<std::vector<Box<NDIM>>> result;
    long coarse_size = 0;
    for (int i = 0; i < coarse_level->getNumberOfPatches(); ++i)
    {
        BoxList<NDIM> coarse_box_list;
        coarse_box_list.addItem(coarse_level->getBoxForPatch(i));
        coarse_size += coarse_box_list.getFirstItem().size();
        finer_boxes.findOverlapIndices(finer_box_indices, coarse_level->getBoxForPatch(i));
        for (const int j : finer_box_indices)
        {
            coarse_box_list.removeIntersections(finer_boxes.getBoxes()[j]);
        }

        result.emplace_back();
        typename tbox::List<Box<NDIM>>::Iterator it(coarse_box_list);
        while (it)
        {
            result.back().push_back(*it);
            combined_size += (*it).size();
            it++;
        }
    }

    TBOX_ASSERT(coarse_size == combined_size);

    return result;
}

// Gather the markers of all processes.
std::tuple<std::vector<double>, std::vector<double>, std::vector<int>>
collect_markers(const std::vector<double>& local_positions,
                const std::vector<double>& local_velocities,
                const std::vector<int>& local_indices)
{
    const auto n_procs = IBTK_MPI::getNodes();
    const int num_indices = local_indices.size();
    std::vector<int> counts(n_procs);
    int ierr = MPI_Allgather(&num_indices, 1, MPI_INT, counts.data(), 1, MPI_INT, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    // We want the first entry to be zero and the last to be the sum
    std::vector<int> offsets(n_procs + 1);
    std::partial_sum(counts.begin(), counts.end(), offsets.begin() + 1);
    std::vector<int> vector_counts;
    std::vector<int> vector_offsets;
    for (int r = 0; r < n_procs; ++r)
    {
        vector_counts.push_back(counts[r] * NDIM);
        vector_offsets.push_back(offsets[r] * NDIM);
    }
    const auto num_total_indices = offsets.back();

    std::vector<int> new_indices(num_total_indices);
    std::vector<double> new_positions(num_total_indices * NDIM);
    std::vector<double> new_velocities(num_total_indices * NDIM);
    ierr = MPI_Allgatherv(local_indices.data(),
                          num_indices,
                          MPI_INT,
                          new_indices.data(),
                          counts.data(),
                          offsets.data(),
                          MPI_INT,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    ierr = MPI_Allgatherv(local_positions.data(),
                          num_indices * NDIM,
                          MPI_DOUBLE,
                          new_positions.data(),
                          vector_counts.data(),
                          vector_offsets.data(),
                          MPI_DOUBLE,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    ierr = MPI_Allgatherv(local_velocities.data(),
                          num_indices * NDIM,
                          MPI_DOUBLE,
                          new_velocities.data(),
                          vector_counts.data(),
                          vector_offsets.data(),
                          MPI_DOUBLE,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);

    return std::make_tuple(std::move(new_positions), std::move(new_velocities), std::move(new_indices));
}

void
do_interpolation(const int data_idx,
                 const std::vector<double>& positions,
                 const Pointer<Patch<NDIM>> patch,
                 const std::string& kernel,
                 std::vector<double>& velocities)
{
    Pointer<PatchData<NDIM>> data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double>> cc_data = data;
    Pointer<SideData<NDIM, double>> sc_data = data;
    const bool is_cc_data = cc_data;
    const bool is_sc_data = sc_data;
    // Only interpolate things within 1 cell of the patch box - we aren't
    // guaranteed to have more ghost data than that
    Box<NDIM> interp_box = data->getBox();
    interp_box.grow(1);
#ifndef NDEBUG
    std::fill(velocities.begin(), velocities.end(), std::numeric_limits<double>::signaling_NaN());
#endif

    if (is_cc_data)
    {
        LEInteractor::interpolate(velocities, NDIM, positions, NDIM, cc_data, patch, interp_box, kernel);
    }
    else if (is_sc_data)
    {
        LEInteractor::interpolate(velocities, NDIM, positions, NDIM, sc_data, patch, interp_box, kernel);
    }
    else
    {
        TBOX_ERROR("not implemented");
    }
#ifndef NDEBUG
    for (const double& v : velocities)
    {
        if (std::isnan(v))
        {
            TBOX_ERROR(
                "One of the marker points inside a patch did not have its velocity set by interpolation. "
                "This most likely means that the marker point is outside of its assigned patch, which "
                "should not happen at this point.");
        }
    }
#endif
}

// Tags of the messages between neighboring processes, on the communicator of the MarkerPatchHierarchy.
constexpr int PRUNE_COUNT_TAG = 1;
constexpr int PRUNE_INT_TAG = 2;
constexpr int PRUNE_DOUBLE_TAG = 3;
constexpr int REQUEST_COUNT_TAG = 4;
constexpr int REQUEST_INT_TAG = 5;
constexpr int REQUEST_DOUBLE_TAG = 6;
constexpr int REPLY_TAG = 7;

inline MPI_Datatype
mpi_type(const int*)
{
    return MPI_INT;
}

inline MPI_Datatype
mpi_type(const double*)
{
    return MPI_DOUBLE;
}

// Send send_counts[i] to ranks[i] and return the counts that the ranks send to this process.
std::vector<int>
exchange_counts(const MPI_Comm comm, const std::vector<int>& ranks, const int tag, const std::vector<int>& send_counts)
{
    const int n = static_cast<int>(ranks.size());
    std::vector<int> recv_counts(n, 0);
    std::vector<MPI_Request> requests(2 * n);
    for (int i = 0; i < n; ++i)
    {
        int ierr = MPI_Irecv(&recv_counts[i], 1, MPI_INT, ranks[i], tag, comm, &requests[i]);
        TBOX_ASSERT(ierr == 0);
        ierr = MPI_Isend(&send_counts[i], 1, MPI_INT, ranks[i], tag, comm, &requests[n + i]);
        TBOX_ASSERT(ierr == 0);
    }
    const int ierr = MPI_Waitall(2 * n, requests.data(), MPI_STATUSES_IGNORE);
    TBOX_ASSERT(ierr == 0);
    return recv_counts;
}

// Send the nonempty arrays send[i] to ranks[i], and receive recv_sizes[i] values from each rank with a nonzero size.
template <typename T>
std::vector<std::vector<T>>
exchange_arrays(const MPI_Comm comm,
                const std::vector<int>& ranks,
                const int tag,
                const std::vector<std::vector<T>>& send,
                const std::vector<int>& recv_sizes)
{
    const int n = static_cast<int>(ranks.size());
    std::vector<std::vector<T>> recv(n);
    std::vector<MPI_Request> requests;
    requests.reserve(2 * n);
    for (int i = 0; i < n; ++i)
    {
        if (recv_sizes[i] > 0)
        {
            recv[i].resize(recv_sizes[i]);
            requests.emplace_back();
            const int ierr = MPI_Irecv(
                recv[i].data(), recv_sizes[i], mpi_type(recv[i].data()), ranks[i], tag, comm, &requests.back());
            TBOX_ASSERT(ierr == 0);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        if (!send[i].empty())
        {
            requests.emplace_back();
            const int ierr = MPI_Isend(send[i].data(),
                                       static_cast<int>(send[i].size()),
                                       mpi_type(send[i].data()),
                                       ranks[i],
                                       tag,
                                       comm,
                                       &requests.back());
            TBOX_ASSERT(ierr == 0);
        }
    }
    const int ierr = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    TBOX_ASSERT(ierr == 0);
    return recv;
}

// Return true if the position X, which lies in the bounding box of the physical domain, is in a notch of the domain:
// the level 0 cell that contains X, in the index frame of the bounding box, is in no box of the physical domain.
bool
in_domain_notch(const double* const X, const Pointer<CartesianGridGeometry<NDIM>>& grid_geom)
{
    const BoxArray<NDIM>& domain = grid_geom->getPhysicalDomain();
    if (domain.getNumberOfBoxes() == 1)
    {
        return false;
    }
    const Box<NDIM> bounding_box = domain.getBoundingBox();
    const hier::Index<NDIM> index = IndexUtilities::getCellIndex(X,
                                                                 grid_geom->getXLower(),
                                                                 grid_geom->getXUpper(),
                                                                 grid_geom->getDx(),
                                                                 bounding_box.lower(),
                                                                 bounding_box.upper());
    for (int b = 0; b < domain.getNumberOfBoxes(); ++b)
    {
        if (domain[b].contains(index))
        {
            return false;
        }
    }
    return true;
}

// Return true if the position is outside the bounding box of the physical domain or in a notch of the domain.
bool
is_outside_domain(const IBTK::Point& position, const Pointer<CartesianGridGeometry<NDIM>>& grid_geom)
{
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (position[d] < x_lower[d] || x_upper[d] <= position[d])
        {
            return true;
        }
    }
    return in_domain_notch(position.data(), grid_geom);
}

// If the position X, which lies in the bounding box of the physical domain, is in a notch, move it into the box of
// the domain that is nearest in the maximum norm of the distance in cells (the lowest numbered box among equals), to
// the nearest point of that box that is at least sqrt(eps) times the box length from each of its faces, as
// IndexUtilities::clampToDomain(). Measure in cells, not lengths, so that the chosen box is at most one cell from the
// marker's patch when dx differs by direction.
void
move_out_of_domain_notch(double* const X, const Pointer<CartesianGridGeometry<NDIM>>& grid_geom)
{
    if (!in_domain_notch(X, grid_geom))
    {
        return;
    }
    const BoxArray<NDIM>& domain = grid_geom->getPhysicalDomain();
    const Box<NDIM> bounding_box = domain.getBoundingBox();
    const double* const x_lower = grid_geom->getXLower();
    const double* const dx = grid_geom->getDx();
    auto face = [&](const int b, const unsigned int d, const bool upper)
    {
        const int i = upper ? domain[b].upper(d) + 1 : domain[b].lower(d);
        return x_lower[d] + (i - bounding_box.lower(d)) * dx[d];
    };
    int nearest_box = -1;
    double nearest_distance = std::numeric_limits<double>::max();
    for (int b = 0; b < domain.getNumberOfBoxes(); ++b)
    {
        double distance = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            distance = std::max({ distance, (face(b, d, false) - X[d]) / dx[d], (X[d] - face(b, d, true)) / dx[d] });
        }
        if (distance < nearest_distance)
        {
            nearest_box = b;
            nearest_distance = distance;
        }
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const double lower = face(nearest_box, d, false);
        const double upper = face(nearest_box, d, true);
        const double offset = std::sqrt(std::numeric_limits<double>::epsilon()) * (upper - lower);
        const double x_min = std::max(lower + offset, std::nextafter(lower, upper));
        const double x_max = std::min(upper - offset, std::nextafter(upper, lower));
        X[d] = std::max(x_min, std::min(X[d], x_max));
    }
#ifndef NDEBUG
    TBOX_ASSERT(domain[nearest_box].contains(IndexUtilities::getCellIndex(X,
                                                                          grid_geom->getXLower(),
                                                                          grid_geom->getXUpper(),
                                                                          grid_geom->getDx(),
                                                                          bounding_box.lower(),
                                                                          bounding_box.upper())));
#endif
}
} // namespace

MarkerPatch::MarkerPatch(const Box<NDIM>& patch_box,
                         const std::vector<Box<NDIM>>& nonoverlapping_patch_boxes,
                         const Pointer<CartesianGridGeometry<NDIM>>& grid_geom,
                         const IntVector<NDIM>& ratio)
    : d_patch_box(patch_box), d_nonoverlapping_patch_boxes(nonoverlapping_patch_boxes)
{
    for (const Box<NDIM>& box : nonoverlapping_patch_boxes)
    {
        TBOX_ASSERT(patch_box.contains(box));
    }
    d_domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain().getBoundingBox(), ratio);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_x_lo[d] = grid_geom->getXLower()[d];
        d_x_up[d] = grid_geom->getXUpper()[d];
        d_dx[d] = grid_geom->getDx()[d] / ratio(d);
    }
}

void
MarkerPatch::insert(const int& index, const IBTK::Point& position, const IBTK::Vector& velocity)
{
    TBOX_ASSERT(index >= 0);
    const auto p = std::lower_bound(d_indices.begin(), d_indices.end(), index) - d_indices.begin();
    d_indices.insert(d_indices.begin() + p, index);
    d_positions.insert(d_positions.begin() + NDIM * p, position.data(), position.data() + position.size());
    d_velocities.insert(d_velocities.begin() + NDIM * p, velocity.data(), velocity.data() + velocity.size());
}

bool
MarkerPatch::contains(const IBTK::Point& position) const
{
    const auto index = IndexUtilities::getCellIndex(
        position.data(), d_x_lo.data(), d_x_up.data(), d_dx.data(), d_domain_box.lower(), d_domain_box.upper());
    if (!d_patch_box.contains(index)) return false;

    for (const auto& box : d_nonoverlapping_patch_boxes)
        if (box.contains(index)) return true;
    return false;
}

std::tuple<std::vector<int>, std::vector<IBTK::Point>, std::vector<IBTK::Vector>>
MarkerPatch::prune()
{
    std::vector<int> indices;
    std::vector<IBTK::Point> positions;
    std::vector<IBTK::Vector> velocities;
    // Prune from the back to the front to avoid moving stored markers
    for (int index = size() - 1; index >= 0; index--)
    {
        IBTK::Point point;
        IBTK::Point velocity;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            point[d] = d_positions[index * NDIM + d];
            velocity[d] = d_velocities[index * NDIM + d];
        }
        if (!contains(point))
        {
            // emplace so that things remain sorted
            indices.insert(indices.begin(), d_indices[index]);
            d_indices.erase(d_indices.begin() + index);
            auto p0 = d_positions.begin() + NDIM * index;
            auto p1 = d_positions.begin() + NDIM * (index + 1);
            positions.emplace(positions.begin(), &*p0);
            d_positions.erase(p0, p1);
            auto v0 = d_velocities.begin() + NDIM * index;
            auto v1 = d_velocities.begin() + NDIM * (index + 1);
            velocities.emplace(velocities.begin(), &*v0);
            d_velocities.erase(v0, v1);
        }
    }

    return std::make_tuple(std::move(indices), std::move(positions), std::move(velocities));
}

std::tuple<int, IBTK::Point, IBTK::Vector>
MarkerPatch::operator[](const unsigned int local_index) const
{
#ifndef NDEBUG
    TBOX_ASSERT(local_index < size());
#endif
    IBTK::Point position(&d_positions[local_index * NDIM]);
    IBTK::Point velocity(&d_velocities[local_index * NDIM]);
    return std::make_tuple(d_indices[local_index], position, velocity);
}

std::size_t
MarkerPatch::size() const
{
#ifndef NDEBUG
    TBOX_ASSERT(d_positions.size() == d_indices.size() * NDIM);
    TBOX_ASSERT(d_velocities.size() == d_indices.size() * NDIM);
#endif
    return d_indices.size();
}

MarkerPatchHierarchy::MarkerPatchHierarchy(const std::string& name,
                                           Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                           const std::vector<IBTK::Point>& positions,
                                           const std::vector<IBTK::Point>& velocities,
                                           const bool register_for_restart)
    : d_object_name(name),
      d_register_for_restart(register_for_restart),
      d_num_markers(positions.size()),
      d_hierarchy(patch_hierarchy),
      d_markers_outside_domain(
          Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()),
          std::vector<Box<NDIM>>{ Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()) },
          d_hierarchy->getGridGeometry(),
          IntVector<NDIM>(1))
{
    const int ierr = MPI_Comm_dup(IBTK_MPI::getCommunicator(), &d_comm);
    TBOX_ASSERT(ierr == 0);

    auto set_timer = [&](const char* name) { return TimerManager::getManager()->getTimer(name); };
    t_reinit = set_timer("IBTK::MarkerPatchHierarchy::reinit()");
    t_collect_all_markers = set_timer("IBTK::MarkerPatchHierarchy::collectAllMarkers()");
    t_forward_euler_step = set_timer("IBTK::MarkerPatchHierarchy::forwardEulerStep()");
    t_midpoint_step = set_timer("IBTK::MarkerPatchHierarchy::midpointStep()");
    t_trapezoidal_step = set_timer("IBTK::MarkerPatchHierarchy::trapezoidalStep()");
    t_prune_and_redistribute = set_timer("IBTK::MarkerPatchHierarchy::pruneAndRedistribute()");

    // Markers are special in that they are not always set up - in particular,
    // someone might want to add them later in a simulation. One possibility
    // is that we restart at step N and then we add markers at step N + M: we
    // have a restart file but we don't have a marker database in it.
    auto* restart_manager = RestartManager::getManager();
    auto restart_db = restart_manager->getRootDatabase();
    if (register_for_restart)
    {
        restart_manager->registerRestartItem(d_object_name, this);
    }
    if (restart_manager->isFromRestart() && restart_db->keyExists(d_object_name))
    {
        getFromDatabase(restart_db->getDatabase(d_object_name));
    }
    else
    {
        reinit(positions, velocities);
    }
}

MarkerPatchHierarchy::~MarkerPatchHierarchy()
{
    if (d_register_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    // An object that outlives MPI_Finalize() must not call MPI.
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized && d_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&d_comm);
    }
}

void
MarkerPatchHierarchy::buildExchangeData()
{
    const auto rank = IBTK_MPI::getRank();
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_x_lo[d] = grid_geom->getXLower()[d];
        d_x_up[d] = grid_geom->getXUpper()[d];
    }

    // The ownership of cells, from the parts of the patch boxes that no finer level covers:
    d_ownership.clear();
    d_ownership.resize(finest_ln + 1);
    for (int ln = finest_ln; ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM>> current_level = d_hierarchy->getPatchLevel(ln);
        LevelOwnership& ownership = d_ownership[ln];
        const IntVector<NDIM>& ratio = current_level->getRatio();
        if (ln == finest_ln)
        {
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                ownership.patch_boxes.emplace_back(1, current_level->getBoxForPatch(i));
            }
        }
        else
        {
            ownership.patch_boxes =
                compute_nonoverlapping_patch_boxes(current_level, d_hierarchy->getPatchLevel(ln + 1));
        }
        ownership.domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain().getBoundingBox(), ratio);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ownership.dx[d] = grid_geom->getDx()[d] / ratio(d);
        }

        int num_boxes = 0;
        for (const auto& boxes : ownership.patch_boxes)
        {
            num_boxes += static_cast<int>(boxes.size());
        }
        BoxArray<NDIM> all_boxes(num_boxes);
        int box_number = 0;
        int local_patch_number = 0;
        ownership.local_patch.assign(current_level->getNumberOfPatches(), -1);
        for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
        {
            if (rank == current_level->getMappingForPatch(i))
            {
                ownership.local_patch[i] = local_patch_number++;
            }
            for (const Box<NDIM>& box : ownership.patch_boxes[i])
            {
                all_boxes[box_number++] = box;
                ownership.box_rank.push_back(current_level->getMappingForPatch(i));
                ownership.box_patch.push_back(i);
            }
        }
        ownership.boxes = detail::BoxArrayIndex(all_boxes, std::vector<IntVector<NDIM>>(1, IntVector<NDIM>(0)));
#ifndef NDEBUG
        // Each cell has at most one owner on a level:
        std::vector<int> overlaps;
        for (int k = 0; k < num_boxes; ++k)
        {
            ownership.boxes.findOverlapIndices(overlaps, all_boxes[k]);
            TBOX_ASSERT(overlaps.size() == 1);
        }
#endif
    }

    // The neighbors: the ranks of the patches, on any level and including periodic images, within two cells (of the
    // coarser of the two levels) of a patch of this process.
    std::vector<detail::BoxArrayIndex> level_indices;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level_indices.push_back(detail::periodic_level_box_index(*level));
    }
    d_neighbor_ranks.clear();
    std::vector<int> overlaps;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (int i = 0; i < level->getNumberOfPatches(); ++i)
        {
            if (rank != level->getMappingForPatch(i))
            {
                continue;
            }
            for (int other_ln = 0; other_ln <= finest_ln; ++other_ln)
            {
                Pointer<PatchLevel<NDIM>> other_level = d_hierarchy->getPatchLevel(other_ln);
                Box<NDIM> box = level->getBoxForPatch(i);
                IntVector<NDIM> growth(2);
                if (other_ln > ln)
                {
                    IntVector<NDIM> ratio = other_level->getRatio();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        ratio(d) /= level->getRatio()(d);
                    }
                    box.refine(ratio);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        growth(d) = 2 * ratio(d);
                    }
                }
                else if (other_ln < ln)
                {
                    IntVector<NDIM> ratio = level->getRatio();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        ratio(d) /= other_level->getRatio()(d);
                    }
                    box.coarsen(ratio);
                }
                box.grow(growth);
                level_indices[other_ln].findOverlapIndices(overlaps, box);
                for (const int j : overlaps)
                {
                    const int other_rank = other_level->getMappingForPatch(j);
                    if (other_rank != rank)
                    {
                        d_neighbor_ranks.push_back(other_rank);
                    }
                }
            }
        }
    }
    std::sort(d_neighbor_ranks.begin(), d_neighbor_ranks.end());
    d_neighbor_ranks.erase(std::unique(d_neighbor_ranks.begin(), d_neighbor_ranks.end()), d_neighbor_ranks.end());

#ifndef NDEBUG
    // The neighbor relation must be symmetric, or the processes would not post matching receives:
    const int num_ranks = IBTK_MPI::getNodes();
    std::vector<int> is_neighbor(num_ranks, 0), neighbors_of(num_ranks, 0);
    for (const int r : d_neighbor_ranks)
    {
        is_neighbor[r] = 1;
    }
    const int ierr = MPI_Alltoall(is_neighbor.data(), 1, MPI_INT, neighbors_of.data(), 1, MPI_INT, d_comm);
    TBOX_ASSERT(ierr == 0);
    for (int r = 0; r < num_ranks; ++r)
    {
        if (is_neighbor[r] != neighbors_of[r])
        {
            TBOX_ERROR(d_object_name << "::reinit():\n  the neighbor relation of ranks " << rank << " and " << r
                                     << " is not symmetric\n");
        }
    }
#endif
}

MarkerPatchHierarchy::Owner
MarkerPatchHierarchy::findOwner(const double* const X, std::vector<int>& found) const
{
    for (int ln = static_cast<int>(d_ownership.size()) - 1; ln >= 0; --ln)
    {
        const LevelOwnership& ownership = d_ownership[ln];
        const hier::Index<NDIM> index = IndexUtilities::getCellIndex(X,
                                                                     d_x_lo.data(),
                                                                     d_x_up.data(),
                                                                     ownership.dx.data(),
                                                                     ownership.domain_box.lower(),
                                                                     ownership.domain_box.upper());
        ownership.boxes.findOverlapIndices(found, Box<NDIM>(index, index));
        if (!found.empty())
        {
            Owner owner;
            owner.level = ln;
            owner.rank = ownership.box_rank[found[0]];
            owner.patch = ownership.box_patch[found[0]];
            return owner;
        }
    }
    return Owner();
}

void
MarkerPatchHierarchy::checkHierarchyIsCurrent(const char* const function) const
{
    if (!d_layout_key.matches(*d_hierarchy, 0, d_hierarchy->getFinestLevelNumber()))
    {
        TBOX_ERROR(d_object_name << "::" << function
                                 << "():\n  the layout of the patch hierarchy changed since the last call to reinit()\n"
                                    "  call reinit() after every regrid\n");
    }
}

void
MarkerPatchHierarchy::reinit(const std::vector<IBTK::Point>& positions, const std::vector<IBTK::Point>& velocities)
{
    TBOX_ASSERT(positions.size() == velocities.size());
    IBTK_TIMER_START(t_reinit);
    d_num_markers = positions.size();
    const auto rank = IBTK_MPI::getRank();
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(positions.size());

    // The ownership data depend only on the layout of the levels:
    if (d_layout_key.changed(*d_hierarchy, 0, finest_ln))
    {
        buildExchangeData();
    }

    d_marker_patches.clear();
    d_marker_patches.resize(finest_ln + 1);
    for (int ln = finest_ln; ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM>> current_level = d_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = current_level->getRatio();
        for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
        {
            if (rank == current_level->getMappingForPatch(i))
            {
                d_marker_patches[ln].emplace_back(
                    current_level->getBoxForPatch(i), d_ownership[ln].patch_boxes[i], grid_geom, ratio);
            }
        }
    }

    d_markers_outside_domain.d_indices.resize(0);
    d_markers_outside_domain.d_positions.resize(0);
    d_markers_outside_domain.d_velocities.resize(0);

    // Every process finds the owner of every marker. The owner stores the marker, and process 0 stores the markers that
    // are outside of the domain.
    std::vector<int> found;
    for (unsigned int k = 0; k < positions.size(); ++k)
    {
        const Owner owner = findOwner(positions[k].data(), found);
        if (owner.level >= 0)
        {
            if (owner.rank == rank)
            {
                marker_emplaced[k] = true;
                d_marker_patches[owner.level][d_ownership[owner.level].local_patch[owner.patch]].insert(
                    k, positions[k], velocities[k]);
                ++num_emplaced_markers;
            }
        }
        else if (rank == 0 && is_outside_domain(positions[k], grid_geom))
        {
            marker_emplaced[k] = true;
            IBTK::Vector v;
            v.fill(0);
            d_markers_outside_domain.insert(k, positions[k], v);
            ++num_emplaced_markers;
        }
    }

    num_emplaced_markers = IBTK_MPI::sumReduction(num_emplaced_markers);
    TBOX_ASSERT(num_emplaced_markers == positions.size());
    // Do one more expensive check in debug mode:
#ifndef NDEBUG
    std::vector<unsigned int> marker_check(positions.size());
    for (unsigned int k = 0; k < marker_check.size(); ++k)
        if (marker_emplaced[k]) marker_check[k] += 1;
    const int ierr = MPI_Allreduce(
        MPI_IN_PLACE, marker_check.data(), marker_check.size(), MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    for (unsigned int k = 0; k < getNumberOfMarkers(); ++k)
    {
        if (marker_check[k] != 1)
        {
            TBOX_ERROR(d_object_name
                       << ": Marker point " << k << " is presently owned by " << marker_check[k]
                       << " patches. The most likely cause of this error is that the CFL number is greater than 1.");
        }
    }
#endif
    IBTK_TIMER_STOP(t_reinit);
}

void
MarkerPatchHierarchy::putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    TBOX_ASSERT(d_num_markers <= std::numeric_limits<int>::max());
    db->putInteger("num_markers", int(d_num_markers));

    auto put_marker_patch = [&](const MarkerPatch& marker_patch, const std::string& prefix)
    {
        db->putInteger(prefix + "_num_markers", static_cast<int>(marker_patch.d_indices.size()));
        // Yet another SAMRAI bug: we are not allowed to store zero-length
        // arrays in the database so we have to special-case that here ourselves
        if (marker_patch.d_indices.size() > 0)
        {
            db->putIntegerArray(
                prefix + "_indices", marker_patch.d_indices.data(), static_cast<int>(marker_patch.d_indices.size()));
            db->putDoubleArray(prefix + "_positions",
                               marker_patch.d_positions.data(),
                               static_cast<int>(marker_patch.d_positions.size()));
            db->putDoubleArray(prefix + "_velocities",
                               marker_patch.d_velocities.data(),
                               static_cast<int>(marker_patch.d_velocities.size()));
        }
    };

    int marker_patch_num = 0;
    for (const auto& level_marker_patches : d_marker_patches)
    {
        for (const auto& marker_patch : level_marker_patches)
        {
            const std::string key_prefix = "patch_" + std::to_string(marker_patch_num);
            put_marker_patch(marker_patch, key_prefix);
            ++marker_patch_num;
        }
    }

    put_marker_patch(d_markers_outside_domain, "outside_domain_");
}

void
MarkerPatchHierarchy::writeH5Part(const std::string& filename,
                                  const int /*time_step*/,
                                  const double simulation_time,
                                  const bool write_velocities) const
{
    const auto pair = collectAllMarkers();
    const auto& positions = pair.first;
    const auto& velocities = pair.second;

    const std::array<std::string, NDIM> position_datasets{ { "x",
                                                             "y",
#if NDIM == 3
                                                             "z"
#endif
    } };

    // technically momentum, but marker points don't have mass so this is the
    // nearest equivalent for us
    const std::array<std::string, NDIM> velocity_datasets{ { "px",
                                                             "py",
#if NDIM == 3
                                                             "pz"
#endif
    } };

    // Make these files also compatible with SAMRAI by encoding their types
    // in the manner perscribed by HDFDatabase.C
    auto set_samrai_attribute = [](const hid_t dataset_id, const int type_key)
    {
        const hid_t attribute_id = H5Screate(H5S_SCALAR);
        TBOX_ASSERT(attribute_id >= 0);
        const auto samrai_attribute_type = H5T_STD_I8BE;
        const hid_t attribute =
            H5Acreate(dataset_id, "Type", samrai_attribute_type, attribute_id, H5P_DEFAULT, H5P_DEFAULT);
        TBOX_ASSERT(attribute >= 0);
        herr_t status = H5Awrite(attribute, H5T_NATIVE_INT, &type_key);
        TBOX_ASSERT(status == 0);
        status = H5Aclose(attribute);
        TBOX_ASSERT(status == 0);
        status = H5Sclose(attribute_id);
        TBOX_ASSERT(status == 0);
    };
    // These values are defined in HDFDatabase.C
    constexpr int samrai_int_array = 7;
    constexpr int samrai_double_array = 5;

    const hsize_t dims[1]{ getNumberOfMarkers() };
    const hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);

    if (IBTK_MPI::getRank() == 0)
    {
        hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
        {
            TBOX_ERROR("An error occurred when calling H5Fcreate() with filename " << filename << std::endl);
        }
        hid_t group_id = H5Gcreate2(file_id, "/Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        TBOX_ASSERT(group_id >= 0);
        // Save the current time in a way VisIt can understand:
        {
            const hid_t attribute_id = H5Screate(H5S_SCALAR);
            TBOX_ASSERT(attribute_id >= 0);
            const hid_t attribute =
                H5Acreate(group_id, "time", H5T_NATIVE_DOUBLE, attribute_id, H5P_DEFAULT, H5P_DEFAULT);
            TBOX_ASSERT(attribute >= 0);
            herr_t status = H5Awrite(attribute, H5T_NATIVE_DOUBLE, &simulation_time);
            TBOX_ASSERT(status == 0);
            status = H5Aclose(attribute);
            TBOX_ASSERT(status == 0);
            status = H5Sclose(attribute_id);
            TBOX_ASSERT(status == 0);
        }

        const hid_t id_dataset_id =
            H5Dcreate2(group_id, "id", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        TBOX_ASSERT(id_dataset_id >= 0);
        set_samrai_attribute(id_dataset_id, samrai_int_array);
        std::vector<long> ids(getNumberOfMarkers());
        std::iota(ids.begin(), ids.end(), 0l);
        herr_t status = H5Dwrite(id_dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
        TBOX_ASSERT(status == 0);
        status = H5Dclose(id_dataset_id);
        TBOX_ASSERT(status == 0);

        // H5Part assumes all properties are in 1D arrays so we need to unpack first
        std::vector<double> unpacked_positions(getNumberOfMarkers());
        std::vector<double> unpacked_velocities(write_velocities ? getNumberOfMarkers() : 0u);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            for (unsigned int k = 0; k < getNumberOfMarkers(); ++k)
            {
                unpacked_positions[k] = positions[k][d];
                if (write_velocities)
                {
                    unpacked_velocities[k] = velocities[k][d];
                }
            }

            const hid_t position_dataset_id = H5Dcreate2(group_id,
                                                         position_datasets[d].c_str(),
                                                         H5T_NATIVE_DOUBLE,
                                                         dataspace_id,
                                                         H5P_DEFAULT,
                                                         H5P_DEFAULT,
                                                         H5P_DEFAULT);
            set_samrai_attribute(position_dataset_id, samrai_double_array);
            TBOX_ASSERT(position_dataset_id >= 0);
            status = H5Dwrite(
                position_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unpacked_positions.data());
            TBOX_ASSERT(status == 0);
            status = H5Dclose(position_dataset_id);
            TBOX_ASSERT(status == 0);
            if (write_velocities)
            {
                const hid_t velocity_dataset_id = H5Dcreate2(group_id,
                                                             velocity_datasets[d].c_str(),
                                                             H5T_NATIVE_DOUBLE,
                                                             dataspace_id,
                                                             H5P_DEFAULT,
                                                             H5P_DEFAULT,
                                                             H5P_DEFAULT);
                status = H5Dwrite(
                    velocity_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unpacked_velocities.data());
                set_samrai_attribute(velocity_dataset_id, samrai_double_array);
                TBOX_ASSERT(status == 0);
                status = H5Dclose(velocity_dataset_id);
                TBOX_ASSERT(status == 0);
            }
        }

        status = H5Gclose(group_id);
        TBOX_ASSERT(status == 0);
        status = H5Fclose(file_id);
        TBOX_ASSERT(status == 0);
    }
}

void
MarkerPatchHierarchy::getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    // the database does not store information present in the patch hierarchy,
    // so reconstruct the marker patches first:
    reinit({}, {});
    d_num_markers = static_cast<std::size_t>(db->getInteger("num_markers"));

    int num_loaded_markers = 0;
    auto get_marker_patch = [&](MarkerPatch& marker_patch, const std::string& prefix)
    {
        const auto num_markers = static_cast<std::size_t>(db->getInteger(prefix + "_num_markers"));
        // No arrays are saved if num_markers == 0
        if (num_markers > 0)
        {
            marker_patch.d_indices.resize(num_markers);
            marker_patch.d_positions.resize(num_markers * NDIM);
            marker_patch.d_velocities.resize(num_markers * NDIM);
            db->getIntegerArray(
                prefix + "_indices", marker_patch.d_indices.data(), static_cast<int>(marker_patch.d_indices.size()));
            db->getDoubleArray(prefix + "_positions",
                               marker_patch.d_positions.data(),
                               static_cast<int>(marker_patch.d_positions.size()));
            db->getDoubleArray(prefix + "_velocities",
                               marker_patch.d_velocities.data(),
                               static_cast<int>(marker_patch.d_velocities.size()));
        }
        num_loaded_markers += static_cast<int>(num_markers);
    };

    int marker_patch_num = 0;
    for (auto& level_marker_patches : d_marker_patches)
    {
        for (auto& marker_patch : level_marker_patches)
        {
            const std::string key_prefix = "patch_" + std::to_string(marker_patch_num);
            get_marker_patch(marker_patch, key_prefix);
            ++marker_patch_num;
        }
    }

    get_marker_patch(d_markers_outside_domain, "outside_domain_");
    num_loaded_markers = IBTK_MPI::sumReduction(num_loaded_markers);
    TBOX_ASSERT(num_loaded_markers == static_cast<int>(d_num_markers));
}

const MarkerPatch&
MarkerPatchHierarchy::getMarkerPatch(const int ln, const int local_patch_num) const
{
    TBOX_ASSERT(ln < static_cast<int>(d_marker_patches.size()));
    TBOX_ASSERT(local_patch_num < static_cast<int>(d_marker_patches[ln].size()));
    return d_marker_patches[ln][local_patch_num];
}

std::size_t
MarkerPatchHierarchy::getNumberOfMarkers() const
{
    return d_num_markers;
}

std::pair<std::vector<IBTK::Point>, std::vector<IBTK::Vector>>
MarkerPatchHierarchy::collectAllMarkers() const
{
    IBTK_TIMER_START(t_collect_all_markers);
    std::vector<double> local_positions;
    std::vector<double> local_velocities;
    std::vector<int> local_indices;

    auto extract_markers = [&](const MarkerPatch& marker_patch)
    {
        for (unsigned int k = 0; k < marker_patch.size(); ++k)
        {
            const auto marker_point = marker_patch[k];
            local_indices.push_back(std::get<0>(marker_point));
            local_positions.insert(
                local_positions.end(), std::get<1>(marker_point).data(), std::get<1>(marker_point).data() + NDIM);
            local_velocities.insert(
                local_velocities.end(), std::get<2>(marker_point).data(), std::get<2>(marker_point).data() + NDIM);
        }
    };

    for (const auto& level_marker_patches : d_marker_patches)
    {
        for (const auto& marker_patch : level_marker_patches)
        {
            extract_markers(marker_patch);
        }
    }
    extract_markers(d_markers_outside_domain);

    auto result = collect_markers(local_positions, local_velocities, local_indices);
    const auto& global_positions = std::get<0>(result);
    const auto& global_velocities = std::get<1>(result);
    const auto& global_indices = std::get<2>(result);
    TBOX_ASSERT(global_positions.size() == getNumberOfMarkers() * NDIM);
    TBOX_ASSERT(global_velocities.size() == getNumberOfMarkers() * NDIM);
    TBOX_ASSERT(global_indices.size() == getNumberOfMarkers());

    std::vector<IBTK::Point> positions(getNumberOfMarkers());
    std::vector<IBTK::Vector> velocities(getNumberOfMarkers());
    for (unsigned int k = 0; k < getNumberOfMarkers(); ++k)
    {
        positions[global_indices[k]] = IBTK::Point(&global_positions[k * NDIM]);
        velocities[global_indices[k]] = IBTK::Vector(&global_velocities[k * NDIM]);
    }

    IBTK_TIMER_STOP(t_collect_all_markers);
    return std::make_pair(std::move(positions), std::move(velocities));
}

void
MarkerPatchHierarchy::setVelocities(const int u_idx, const std::string& kernel)
{
    checkHierarchyIsCurrent("setVelocities");
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM>> current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                Pointer<Patch<NDIM>> patch = current_level->getPatch(p);
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                do_interpolation(u_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);
                ++local_patch_num;
            }
        }
    }
}

std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>>
MarkerPatchHierarchy::getLocalMarkerPatches()
{
    const auto rank = IBTK_MPI::getRank();
    std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>> local_marker_patches;
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM>> current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                local_marker_patches.emplace_back(&d_marker_patches[ln][local_patch_num], current_level->getPatch(p));
                ++local_patch_num;
            }
        }
    }
    return local_marker_patches;
}

void
MarkerPatchHierarchy::interpolateInCurrentPatches(
    const int u_idx,
    const std::string& kernel,
    const std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>>& local_marker_patches,
    const std::vector<std::vector<double>>& positions,
    std::vector<std::vector<double>>& velocities)
{
    // Interpolate each marker with the patch whose box contains its position. Another process can own that patch.
    const auto rank = IBTK_MPI::getRank();
    const Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    const IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift();
    const int num_neighbors = static_cast<int>(d_neighbor_ranks.size());

    // The index in local_marker_patches of the marker patch of the first local patch of each level:
    std::vector<unsigned int> first_slot(d_marker_patches.size(), 0);
    for (int ln = static_cast<int>(d_marker_patches.size()) - 2; ln >= 0; --ln)
    {
        first_slot[ln] = first_slot[ln + 1] + static_cast<unsigned int>(d_marker_patches[ln + 1].size());
    }

    // The positions to interpolate with each local patch, and where to put the velocities:
    struct Batch
    {
        std::vector<double> positions;
        std::vector<double*> velocities;
    };
    std::vector<Batch> batches(local_marker_patches.size());

    struct DepartingMarker
    {
        unsigned int patch = 0;
        unsigned int marker = 0;
    };
    std::vector<DepartingMarker> departing_markers;
    for (unsigned int i = 0; i < local_marker_patches.size(); ++i)
    {
        MarkerPatch& marker_patch = *local_marker_patches[i].first;
        velocities[i].resize(positions[i].size());
        for (unsigned int k = 0; k < marker_patch.size(); ++k)
        {
            const IBTK::Point position(&positions[i][NDIM * k]);
            if (marker_patch.contains(position))
            {
                batches[i].positions.insert(batches[i].positions.end(), position.data(), position.data() + NDIM);
                batches[i].velocities.push_back(&velocities[i][NDIM * k]);
            }
            else
            {
                departing_markers.push_back({ i, k });
            }
        }
    }

    // Send each departing marker to the process of the patch whose box contains it. Beyond a
    // physical boundary or in a notch of the domain, the owner is the patch whose box contains the nearest point of the
    // domain. Interpolation uses the ghost values of the owner.
    std::vector<int> found;
    std::vector<std::vector<int>> request_ints(num_neighbors);
    std::vector<std::vector<double>> request_positions(num_neighbors);
    std::vector<std::vector<DepartingMarker>> request_markers(num_neighbors);
    for (const DepartingMarker& departing_marker : departing_markers)
    {
        const unsigned int i = departing_marker.patch;
        double* const velocity = &velocities[i][NDIM * departing_marker.marker];
        std::fill_n(velocity, NDIM, 0.0);
        std::array<double, NDIM> position, clamped_position;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            position[d] = positions[i][NDIM * departing_marker.marker + d];
            if (periodic_shift[d])
            {
                position[d] = wrap_into_periodic_interval(position[d], d_x_lo[d], d_x_up[d]);
            }
            clamped_position[d] =
                periodic_shift[d] ? position[d] : IndexUtilities::clampToDomain(position[d], d_x_lo[d], d_x_up[d]);
        }
        move_out_of_domain_notch(clamped_position.data(), grid_geom);
        const Owner owner = findOwner(clamped_position.data(), found);
        if (owner.level < 0)
        {
            continue;
        }
        if (owner.rank == rank)
        {
            const unsigned int slot = first_slot[owner.level] + d_ownership[owner.level].local_patch[owner.patch];
#ifndef NDEBUG
            TBOX_ASSERT(local_marker_patches[slot].first ==
                        &d_marker_patches[owner.level][d_ownership[owner.level].local_patch[owner.patch]]);
#endif
            batches[slot].positions.insert(batches[slot].positions.end(), position.begin(), position.end());
            batches[slot].velocities.push_back(velocity);
        }
        else
        {
            const auto neighbor = std::lower_bound(d_neighbor_ranks.begin(), d_neighbor_ranks.end(), owner.rank);
            const bool is_neighbor = neighbor != d_neighbor_ranks.end() && *neighbor == owner.rank;
#ifndef NDEBUG
            TBOX_ASSERT(is_neighbor);
#endif
            // In an optimized build, a marker whose owner is not a neighbor has zero velocity:
            if (!is_neighbor)
            {
                continue;
            }
            const auto n = neighbor - d_neighbor_ranks.begin();
            request_ints[n].push_back(owner.level);
            request_ints[n].push_back(owner.patch);
            request_positions[n].insert(request_positions[n].end(), position.begin(), position.end());
            request_markers[n].push_back(departing_marker);
        }
    }

    // Send the requests to the owners, and interpolate the requests that this process owns:
    std::vector<int> request_counts(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        request_counts[n] = static_cast<int>(request_markers[n].size());
    }
    const std::vector<int> owner_counts = exchange_counts(d_comm, d_neighbor_ranks, REQUEST_COUNT_TAG, request_counts);
    std::vector<int> int_sizes(num_neighbors), double_sizes(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        int_sizes[n] = 2 * owner_counts[n];
        double_sizes[n] = NDIM * owner_counts[n];
    }
    const std::vector<std::vector<int>> owner_ints =
        exchange_arrays(d_comm, d_neighbor_ranks, REQUEST_INT_TAG, request_ints, int_sizes);
    const std::vector<std::vector<double>> owner_positions =
        exchange_arrays(d_comm, d_neighbor_ranks, REQUEST_DOUBLE_TAG, request_positions, double_sizes);
    std::vector<std::vector<double>> replies(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        replies[n].resize(double_sizes[n]);
        for (int j = 0; j < owner_counts[n]; ++j)
        {
            const int ln = owner_ints[n][2 * j];
            const int local_patch = d_ownership[ln].local_patch[owner_ints[n][2 * j + 1]];
            if (local_patch < 0)
            {
                TBOX_ERROR(d_object_name << "::interpolateInCurrentPatches():\n  received a request for a patch that "
                                            "this process does not own\n");
            }
            Batch& batch = batches[first_slot[ln] + local_patch];
            batch.positions.insert(
                batch.positions.end(), &owner_positions[n][NDIM * j], &owner_positions[n][NDIM * j] + NDIM);
            batch.velocities.push_back(&replies[n][NDIM * j]);
        }
    }
    for (unsigned int slot = 0; slot < batches.size(); ++slot)
    {
        const Batch& batch = batches[slot];
        if (batch.velocities.empty())
        {
            continue;
        }
        std::vector<double> batch_velocities(batch.positions.size());
        do_interpolation(u_idx, batch.positions, local_marker_patches[slot].second, kernel, batch_velocities);
        for (unsigned int m = 0; m < batch.velocities.size(); ++m)
        {
            std::copy_n(&batch_velocities[NDIM * m], NDIM, batch.velocities[m]);
        }
    }

    // Return the velocities to the processes that requested them:
    std::vector<int> reply_sizes(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        reply_sizes[n] = NDIM * request_counts[n];
    }
    const std::vector<std::vector<double>> request_velocities =
        exchange_arrays(d_comm, d_neighbor_ranks, REPLY_TAG, replies, reply_sizes);
    for (int n = 0; n < num_neighbors; ++n)
    {
        for (int j = 0; j < request_counts[n]; ++j)
        {
            const DepartingMarker& departing_marker = request_markers[n][j];
            std::copy_n(&request_velocities[n][NDIM * j],
                        NDIM,
                        &velocities[departing_marker.patch][NDIM * departing_marker.marker]);
        }
    }
}

void
MarkerPatchHierarchy::forwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    checkHierarchyIsCurrent("forwardEulerStep");
    IBTK_TIMER_START(t_forward_euler_step);
    const std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>> local_marker_patches = getLocalMarkerPatches();

    // 1. Do a forward Euler step:
    std::vector<std::vector<double>> positions;
    for (const auto& local_marker_patch : local_marker_patches)
    {
        MarkerPatch& marker_patch = *local_marker_patch.first;
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            marker_patch.d_positions[i] += dt * marker_patch.d_velocities[i];
        }
        positions.push_back(marker_patch.d_positions);
    }

    // 2. Update the velocities:
    std::vector<std::vector<double>> velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, positions, velocities);
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        local_marker_patches[n].first->d_velocities = velocities[n];
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_forward_euler_step);
}

void
MarkerPatchHierarchy::backwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    checkHierarchyIsCurrent("backwardEulerStep");
    const std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>> local_marker_patches = getLocalMarkerPatches();

    // 1. Do a forward Euler step:
    std::vector<std::vector<double>> new_positions;
    for (const auto& local_marker_patch : local_marker_patches)
    {
        const MarkerPatch& marker_patch = *local_marker_patch.first;
        new_positions.push_back(marker_patch.d_positions);
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            new_positions.back()[i] += dt * marker_patch.d_velocities[i];
        }
    }

    // 2. Compute new velocities:
    std::vector<std::vector<double>> new_velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, new_positions, new_velocities);

    // 3. Update the positions:
    std::vector<std::vector<double>> positions;
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        MarkerPatch& marker_patch = *local_marker_patches[n].first;
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            marker_patch.d_positions[i] += dt * new_velocities[n][i];
        }
        positions.push_back(marker_patch.d_positions);
    }

    // 4. Update the velocities:
    std::vector<std::vector<double>> velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, positions, velocities);
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        local_marker_patches[n].first->d_velocities = velocities[n];
    }

    pruneAndRedistribute();
}

void
MarkerPatchHierarchy::midpointStep(const double dt,
                                   const int u_half_idx,
                                   const int u_new_idx,
                                   const std::string& kernel)
{
    checkHierarchyIsCurrent("midpointStep");
    IBTK_TIMER_START(t_midpoint_step);
    const std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>> local_marker_patches = getLocalMarkerPatches();

    // 1. Do a half of a forward Euler step:
    std::vector<std::vector<double>> half_positions;
    for (const auto& local_marker_patch : local_marker_patches)
    {
        const MarkerPatch& marker_patch = *local_marker_patch.first;
        half_positions.push_back(marker_patch.d_positions);
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            half_positions.back()[i] += 0.5 * dt * marker_patch.d_velocities[i];
        }
    }

    // 2. Interpolate midpoint velocity:
    std::vector<std::vector<double>> half_velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_half_idx, kernel, local_marker_patches, half_positions, half_velocities);

    // 3. Do a midpoint step:
    std::vector<std::vector<double>> positions;
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        MarkerPatch& marker_patch = *local_marker_patches[n].first;
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            marker_patch.d_positions[i] += dt * half_velocities[n][i];
        }
        positions.push_back(marker_patch.d_positions);
    }

    // 4. Interpolate the velocity at the new time:
    std::vector<std::vector<double>> velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, positions, velocities);
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        local_marker_patches[n].first->d_velocities = velocities[n];
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_midpoint_step);
}

void
MarkerPatchHierarchy::trapezoidalStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    checkHierarchyIsCurrent("trapezoidalStep");
    IBTK_TIMER_START(t_trapezoidal_step);
    const std::vector<std::pair<MarkerPatch*, Pointer<Patch<NDIM>>>> local_marker_patches = getLocalMarkerPatches();

    // 1. Do a forward Euler step:
    std::vector<std::vector<double>> new_positions;
    for (const auto& local_marker_patch : local_marker_patches)
    {
        const MarkerPatch& marker_patch = *local_marker_patch.first;
        new_positions.push_back(marker_patch.d_positions);
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            new_positions.back()[i] += dt * marker_patch.d_velocities[i];
        }
    }

    // 2. Interpolate the velocity at the new time:
    std::vector<std::vector<double>> new_velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, new_positions, new_velocities);

    // 3. Do a trapezoidal step:
    std::vector<std::vector<double>> positions;
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        MarkerPatch& marker_patch = *local_marker_patches[n].first;
        for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
        {
            // This is left unfactored so it matches the implementations
            // elsewhere
            marker_patch.d_positions[i] += 0.5 * dt * marker_patch.d_velocities[i] + 0.5 * dt * new_velocities[n][i];
        }
        positions.push_back(marker_patch.d_positions);
    }

    // 4. Update the velocities:
    std::vector<std::vector<double>> velocities(local_marker_patches.size());
    interpolateInCurrentPatches(u_new_idx, kernel, local_marker_patches, positions, velocities);
    for (unsigned int n = 0; n < local_marker_patches.size(); ++n)
    {
        local_marker_patches[n].first->d_velocities = velocities[n];
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_trapezoidal_step);
}

void
MarkerPatchHierarchy::pruneAndRedistribute()
{
    IBTK_TIMER_START(t_prune_and_redistribute);
    const auto rank = IBTK_MPI::getRank();
    const int num_neighbors = static_cast<int>(d_neighbor_ranks.size());

    // 1. Collect all markers which have left their respective patches:
    std::vector<double> moved_positions;
    std::vector<double> moved_velocities;
    std::vector<int> moved_indices;
    for (auto& level_marker_patches : d_marker_patches)
    {
        for (auto& marker_patch : level_marker_patches)
        {
            const auto moved_markers = marker_patch.prune();
            const auto n_moved_markers = std::get<0>(moved_markers).size();
            moved_indices.reserve(moved_indices.size() + n_moved_markers);
            moved_positions.reserve(moved_indices.size() + n_moved_markers);
            moved_velocities.reserve(moved_indices.size() + n_moved_markers);
            for (unsigned int k = 0; k < std::get<0>(moved_markers).size(); ++k)
            {
                moved_indices.push_back(std::get<0>(moved_markers)[k]);
                moved_positions.insert(moved_positions.end(),
                                       std::get<1>(moved_markers)[k].data(),
                                       std::get<1>(moved_markers)[k].data() + NDIM);
                moved_velocities.insert(moved_velocities.end(),
                                        std::get<2>(moved_markers)[k].data(),
                                        std::get<2>(moved_markers)[k].data() + NDIM);
            }
        }
    }

    // 2. Apply periodicity constraints.
    const Pointer<CartesianGridGeometry<NDIM>> grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    const IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift();
    for (unsigned int k = 0; k < moved_indices.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift[d])
            {
                double& X = moved_positions[k * NDIM + d];
                X = wrap_into_periodic_interval(X, domain_x_lower[d], domain_x_upper[d]);
                TBOX_ASSERT(X >= domain_x_lower[d] && X < domain_x_upper[d]);
            }
        }
    }

    // 3. Insert each moved marker into the patch that owns it, or send it to the process of the owner. An unowned
    // marker
    //    outside the domain stays on this process with zero velocity.
    unsigned int num_emplaced_markers = 0;
    std::vector<int> found;
    std::vector<std::vector<int>> send_ints(num_neighbors);
    std::vector<std::vector<double>> send_doubles(num_neighbors);
    for (unsigned int k = 0; k < moved_indices.size(); ++k)
    {
        const IBTK::Point position(&moved_positions[k * NDIM]);
        const IBTK::Vector velocity(&moved_velocities[k * NDIM]);
        const Owner owner = findOwner(position.data(), found);
        if (owner.level < 0)
        {
            if (is_outside_domain(position, grid_geom))
            {
                IBTK::Vector zero_velocity;
                zero_velocity.fill(0.0);
                d_markers_outside_domain.insert(moved_indices[k], position, zero_velocity);
                ++num_emplaced_markers;
            }
        }
        else if (owner.rank == rank)
        {
            d_marker_patches[owner.level][d_ownership[owner.level].local_patch[owner.patch]].insert(
                moved_indices[k], position, velocity);
            ++num_emplaced_markers;
        }
        else
        {
            const auto neighbor = std::lower_bound(d_neighbor_ranks.begin(), d_neighbor_ranks.end(), owner.rank);
            const bool is_neighbor = neighbor != d_neighbor_ranks.end() && *neighbor == owner.rank;
#ifndef NDEBUG
            TBOX_ASSERT(is_neighbor);
#endif
            // In an optimized build, a marker whose owner is not a neighbor is lost:
            if (!is_neighbor)
            {
                continue;
            }
            const auto n = neighbor - d_neighbor_ranks.begin();
            send_ints[n].push_back(moved_indices[k]);
            send_ints[n].push_back(owner.level);
            send_ints[n].push_back(owner.patch);
            send_doubles[n].insert(send_doubles[n].end(), position.data(), position.data() + NDIM);
            send_doubles[n].insert(send_doubles[n].end(), velocity.data(), velocity.data() + NDIM);
        }
    }

    // 4. Exchange markers with the neighbors and insert the markers received.
    std::vector<int> send_counts(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        send_counts[n] = static_cast<int>(send_ints[n].size() / 3);
    }
    const std::vector<int> recv_counts = exchange_counts(d_comm, d_neighbor_ranks, PRUNE_COUNT_TAG, send_counts);
    std::vector<int> int_sizes(num_neighbors), double_sizes(num_neighbors);
    for (int n = 0; n < num_neighbors; ++n)
    {
        int_sizes[n] = 3 * recv_counts[n];
        double_sizes[n] = 2 * NDIM * recv_counts[n];
    }
    const std::vector<std::vector<int>> recv_ints =
        exchange_arrays(d_comm, d_neighbor_ranks, PRUNE_INT_TAG, send_ints, int_sizes);
    const std::vector<std::vector<double>> recv_doubles =
        exchange_arrays(d_comm, d_neighbor_ranks, PRUNE_DOUBLE_TAG, send_doubles, double_sizes);
    for (int n = 0; n < num_neighbors; ++n)
    {
        for (int j = 0; j < recv_counts[n]; ++j)
        {
            const int ln = recv_ints[n][3 * j + 1];
            const int local_patch = d_ownership[ln].local_patch[recv_ints[n][3 * j + 2]];
            if (local_patch < 0)
            {
                TBOX_ERROR(d_object_name << "::pruneAndRedistribute():\n  received a marker for a patch that this "
                                            "process does not own\n");
            }
            const IBTK::Point position(&recv_doubles[n][2 * NDIM * j]);
            const IBTK::Vector velocity(&recv_doubles[n][2 * NDIM * j + NDIM]);
            MarkerPatch& marker_patch = d_marker_patches[ln][local_patch];
#ifndef NDEBUG
            TBOX_ASSERT(marker_patch.contains(position));
#endif
            marker_patch.insert(recv_ints[n][3 * j], position, velocity);
            ++num_emplaced_markers;
        }
    }

    // 5. Check that every point is accounted for.
#ifndef NDEBUG
    {
        unsigned int num_moved_markers = static_cast<unsigned int>(moved_indices.size());
        int ierr =
            MPI_Allreduce(MPI_IN_PLACE, &num_emplaced_markers, 1, MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        ierr = MPI_Allreduce(MPI_IN_PLACE, &num_moved_markers, 1, MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        TBOX_ASSERT(num_emplaced_markers == num_moved_markers);
    }

    std::vector<unsigned int> marker_check(getNumberOfMarkers());
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        for (MarkerPatch& marker_patch : d_marker_patches[ln])
        {
            for (unsigned int k = 0; k < marker_patch.size(); ++k)
            {
                const auto index = std::get<0>(marker_patch[k]);
                TBOX_ASSERT(index < long(getNumberOfMarkers()));
                marker_check[index] += 1;
            }
        }
    }
    for (unsigned int k = 0; k < d_markers_outside_domain.size(); ++k)
    {
        const auto index = std::get<0>(d_markers_outside_domain[k]);
        TBOX_ASSERT(index < long(getNumberOfMarkers()));
        marker_check[index] += 1;
    }
    {
        const int ierr = MPI_Allreduce(
            MPI_IN_PLACE, marker_check.data(), marker_check.size(), MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
    }
    for (unsigned int i = 0; i < getNumberOfMarkers(); ++i)
    {
        if (marker_check[i] != 1)
        {
            TBOX_ERROR(d_object_name
                       << ": Marker point " << i << " is presently owned by " << marker_check[i]
                       << " patches. The most likely cause of this error is that the CFL number is greater than 1.");
        }
    }
#else
    NULL_USE(num_emplaced_markers);
#endif
    IBTK_TIMER_STOP(t_prune_and_redistribute);
}
} // namespace IBTK
