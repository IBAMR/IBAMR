// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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

#include <hdf5.h>
#include <tbox/RestartManager.h>

#include <CartesianGridGeometry.h>
#include <mpi.h>

#include <limits>
#include <numeric>

#include <ibtk/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace
{
static Timer *t_reinit;
static Timer *t_collect_all_markers;
static Timer *t_forward_euler_step;
static Timer *t_midpoint_step;
static Timer *t_trapezoidal_step;
static Timer *t_prune_and_redistribute;

std::vector<std::vector<hier::Box<NDIM> > >
compute_nonoverlapping_patch_boxes(const Pointer<BasePatchLevel<NDIM> >& c_level,
                                   const Pointer<BasePatchLevel<NDIM> >& f_level)
{
    const Pointer<PatchLevel<NDIM> > coarse_level = c_level;
    const Pointer<PatchLevel<NDIM> > fine_level = f_level;
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

    // Remove said boxes from each coarse-level patch:
    const auto rank = IBTK_MPI::getRank();
    std::vector<std::vector<Box<NDIM> > > result;
    long coarse_size = 0;
    for (int i = 0; i < coarse_level->getNumberOfPatches(); ++i)
    {
        BoxList<NDIM> coarse_box_list;
        coarse_box_list.addItem(coarse_level->getBoxForPatch(i));
        coarse_size += coarse_box_list.getFirstItem().size();
        coarse_box_list.removeIntersections(finer_box_list);

        const bool patch_is_local = rank == coarse_level->getMappingForPatch(i);
        if (patch_is_local) result.emplace_back();
        typename tbox::List<Box<NDIM> >::Iterator it(coarse_box_list);
        while (it)
        {
            if (patch_is_local) result.back().push_back(*it);
            combined_size += (*it).size();
            it++;
        }
    }

    TBOX_ASSERT(coarse_size == combined_size);

    return result;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<int> >
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
                 const Pointer<Patch<NDIM> > patch,
                 const std::string& kernel,
                 std::vector<double>& velocities)
{
    Pointer<PatchData<NDIM> > data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > cc_data = data;
    Pointer<SideData<NDIM, double> > sc_data = data;
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
} // namespace

MarkerPatch::MarkerPatch(const Box<NDIM>& patch_box,
                         const std::vector<Box<NDIM> >& nonoverlapping_patch_boxes,
                         const Pointer<CartesianGridGeometry<NDIM> >& grid_geom,
                         const IntVector<NDIM>& ratio)
    : d_patch_box(patch_box), d_nonoverlapping_patch_boxes(nonoverlapping_patch_boxes)
{
    for (const Box<NDIM>& box : nonoverlapping_patch_boxes)
    {
        TBOX_ASSERT(patch_box.contains(box));
    }
    d_domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
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

std::tuple<std::vector<int>, EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> >
MarkerPatch::prune()
{
    std::vector<int> indices;
    EigenAlignedVector<IBTK::Point> positions;
    EigenAlignedVector<IBTK::Vector> velocities;
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
                                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                           const EigenAlignedVector<IBTK::Point>& positions,
                                           const EigenAlignedVector<IBTK::Point>& velocities,
                                           const bool register_for_restart)
    : d_object_name(name),
      d_register_for_restart(register_for_restart),
      d_num_markers(positions.size()),
      d_hierarchy(patch_hierarchy),
      d_markers_outside_domain(
          Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()),
          std::vector<Box<NDIM> >{ Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()) },
          d_hierarchy->getGridGeometry(),
          IntVector<NDIM>(1))
{
    auto set_timer = [&](const char *name)
    { return TimerManager::getManager()->getTimer(name); };
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
}

void
MarkerPatchHierarchy::reinit(const EigenAlignedVector<IBTK::Point>& positions,
                             const EigenAlignedVector<IBTK::Point>& velocities)
{
    TBOX_ASSERT(positions.size() == velocities.size());
    IBTK_TIMER_START(t_reinit);
    d_num_markers = positions.size();
    const auto rank = IBTK_MPI::getRank();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(positions.size());

    auto insert_markers = [&](MarkerPatch& marker_patch)
    {
        for (unsigned int k = 0; k < positions.size(); ++k)
        {
            if (!marker_emplaced[k] && marker_patch.contains(positions[k]))
            {
                marker_emplaced[k] = true;
                marker_patch.insert(k, positions[k], velocities[k]);
                ++num_emplaced_markers;
            }
        }
    };

    d_marker_patches.clear();
    d_marker_patches.resize(d_hierarchy->getFinestLevelNumber() + 1);

    // Assign particles to levels in a top-down way:
    //
    // 1. Compute a set of boxes for the patch which do not intersect any
    //    finer patch level.
    // 2. Emplace particles in the vector corresponding to the local index of a
    //    patch in the present level.
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > finer_level =
            ln == d_hierarchy->getFinestLevelNumber() ? nullptr : d_hierarchy->getPatchLevel(ln + 1);
        const IntVector<NDIM>& ratio = current_level->getRatio();

        // If there is no finer level then each Patch has exactly one
        // nonoverlapping box
        if (!finer_level)
        {
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    const Box<NDIM> box = current_level->getPatch(i)->getBox();
                    d_marker_patches[ln].emplace_back(box, std::vector<Box<NDIM> >{ box }, grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                }
            }
        }
        // otherwise we need to subtract off the boxes on the finer level first.
        else
        {
            const std::vector<std::vector<hier::Box<NDIM> > > nonoverlapping_patch_boxes =
                compute_nonoverlapping_patch_boxes(current_level, finer_level);
            unsigned int local_num = 0;
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    d_marker_patches[ln].emplace_back(
                        current_level->getPatch(i)->getBox(), nonoverlapping_patch_boxes[local_num], grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                    ++local_num;
                }
            }
        }
    }

    // Handle any markers which may be outside of the domain:
    if (rank == 0)
    {
        d_markers_outside_domain.d_indices.resize(0);
        d_markers_outside_domain.d_positions.resize(0);
        d_markers_outside_domain.d_velocities.resize(0);

        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        const double* const domain_x_lower = grid_geom->getXLower();
        const double* const domain_x_upper = grid_geom->getXUpper();
        for (unsigned int k = 0; k < positions.size(); ++k)
        {
            const auto& position = positions[k];
            bool point_outside_domain = false;
            for (unsigned int d = 0; d < NDIM; ++d)
                point_outside_domain =
                    point_outside_domain || ((position[d] < domain_x_lower[d]) || (domain_x_upper[d] <= position[d]));

            if (point_outside_domain)
            {
                marker_emplaced[k] = true;
                IBTK::Vector v;
                v.fill(0);
                d_markers_outside_domain.insert(k, position, v);
                ++num_emplaced_markers;
            }
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
            H5Dcreate2(group_id, "id", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        TBOX_ASSERT(id_dataset_id >= 0);
        set_samrai_attribute(id_dataset_id, samrai_int_array);
        std::vector<int> ids(getNumberOfMarkers());
        herr_t status = H5Dwrite(id_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
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

std::pair<EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> >
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
    if (IBTK_MPI::getRank() == 0)
    {
        extract_markers(d_markers_outside_domain);
    }

    auto result = collect_markers(local_positions, local_velocities, local_indices);
    const auto& global_positions = std::get<0>(result);
    const auto& global_velocities = std::get<1>(result);
    const auto& global_indices = std::get<2>(result);
    TBOX_ASSERT(global_positions.size() == getNumberOfMarkers() * NDIM);
    TBOX_ASSERT(global_velocities.size() == getNumberOfMarkers() * NDIM);
    TBOX_ASSERT(global_indices.size() == getNumberOfMarkers());

    EigenAlignedVector<IBTK::Point> positions(getNumberOfMarkers());
    EigenAlignedVector<IBTK::Vector> velocities(getNumberOfMarkers());
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
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                do_interpolation(u_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);
                ++local_patch_num;
            }
        }
    }
}

void
MarkerPatchHierarchy::forwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    IBTK_TIMER_START(t_forward_euler_step);
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);

                // 1. Do a forward Euler step:
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    marker_patch.d_positions[i] += dt * marker_patch.d_velocities[i];
                }

                // 2. Update the velocities:
                do_interpolation(u_new_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);

                ++local_patch_num;
            }
        }
    }

    pruneAndRedistribute();
}

void
MarkerPatchHierarchy::backwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);

                // 1. Do a forward Euler step:
                std::vector<double> new_positions(marker_patch.d_positions);
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    new_positions[i] = marker_patch.d_positions[i] + dt * marker_patch.d_velocities[i];
                }

                // 2. Compute new velocities:
                std::vector<double> new_velocities(marker_patch.d_velocities);
                do_interpolation(u_new_idx, new_positions, patch, kernel, new_velocities);

                // 3. Update the positions:
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    marker_patch.d_positions[i] += dt * new_velocities[i];
                }

                // 4. Update the velocities:
                do_interpolation(u_new_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);

                ++local_patch_num;
            }
        }
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_forward_euler_step);
}

void
MarkerPatchHierarchy::midpointStep(const double dt,
                                   const int u_half_idx,
                                   const int u_new_idx,
                                   const std::string& kernel)
{
    IBTK_TIMER_START(t_midpoint_step);
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);

                // 1. Do a half of a forward Euler step:
                std::vector<double> half_positions(marker_patch.d_positions);
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    half_positions[i] = marker_patch.d_positions[i] + 0.5 * dt * marker_patch.d_velocities[i];
                }

                // 2. Interpolate midpoint velocity:
                std::vector<double> half_velocities(marker_patch.d_velocities);
                do_interpolation(u_half_idx, half_positions, patch, kernel, half_velocities);

                // 3. Do a midpoint step:
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    marker_patch.d_positions[i] += dt * half_velocities[i];
                }

                // 4. Interpolate the velocity at the new time:
                do_interpolation(u_new_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);

                ++local_patch_num;
            }
        }
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_midpoint_step);
}

void
MarkerPatchHierarchy::trapezoidalStep(const double dt, const int u_new_idx, const std::string& kernel)
{
    IBTK_TIMER_START(t_trapezoidal_step);
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);

                // 1. Do a forward Euler step:
                std::vector<double> new_positions(marker_patch.d_positions);
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    new_positions[i] = marker_patch.d_positions[i] + dt * marker_patch.d_velocities[i];
                }

                // 2. Interpolate the velocity at the new time:
                std::vector<double> new_velocities(marker_patch.d_velocities);
                do_interpolation(u_new_idx, new_positions, patch, kernel, new_velocities);

                // 3. Do a trapezoidal step:
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    // This is left unfactored so it matches the implementations
                    // elsewhere
                    marker_patch.d_positions[i] +=
                        0.5 * dt * marker_patch.d_velocities[i] + 0.5 * dt * new_velocities[i];
                }

                // 4. Update the velocities:
                do_interpolation(u_new_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);

                ++local_patch_num;
            }
        }
    }

    pruneAndRedistribute();
    IBTK_TIMER_STOP(t_trapezoidal_step);
}

void
MarkerPatchHierarchy::pruneAndRedistribute()
{
    IBTK_TIMER_START(t_prune_and_redistribute);
    const auto rank = IBTK_MPI::getRank();

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

    // 2. Communicate data so all processors have all errant marker points:
    std::vector<double> new_positions;
    std::vector<double> new_velocities;
    std::vector<int> new_indices;
    std::tie(new_positions, new_velocities, new_indices) =
        collect_markers(moved_positions, moved_velocities, moved_indices);

    // 3. Apply periodicity constraints.
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    const IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift();
    for (unsigned int k = 0; k < new_indices.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift[d])
            {
                double domain_length = domain_x_upper[d] - domain_x_lower[d];
                double& X = new_positions[k * NDIM + d];
                if (X < domain_x_lower[d]) X += domain_length;
                if (X >= domain_x_upper[d]) X -= domain_length;
                TBOX_ASSERT(X >= domain_x_lower[d] && X < domain_x_upper[d]);
            }
        }
    }

    // 4. Emplace the marker points.
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(new_indices.size());

    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        for (MarkerPatch& marker_patch : d_marker_patches[ln])
        {
            for (unsigned int k = 0; k < new_indices.size(); ++k)
            {
                IBTK::Point position(&new_positions[k * NDIM]);
                IBTK::Vector velocity(&new_velocities[k * NDIM]);
                if (!marker_emplaced[k] && marker_patch.contains(position))
                {
                    marker_emplaced[k] = true;
                    marker_patch.insert(new_indices[k], position, velocity);
                    ++num_emplaced_markers;
                }
            }
        }
    }

    // 5. Account for points outside the computational domain which could not
    //    re-enter the domain via a periodic boundary.
    if (rank == 0)
    {
        for (unsigned int k = 0; k < new_indices.size(); ++k)
        {
            if (!marker_emplaced[k])
            {
                IBTK::Point position(&new_positions[k * NDIM]);
                // Set external velocities to zero
                IBTK::Vector velocity;
                velocity.fill(0.0);
                bool point_outside_domain = false;
                for (unsigned int d = 0; d < NDIM; ++d)
                    point_outside_domain = point_outside_domain ||
                                           ((position[d] < domain_x_lower[d]) || (domain_x_upper[d] <= position[d]));

                if (point_outside_domain)
                {
                    marker_emplaced[k] = true;
                    d_markers_outside_domain.insert(new_indices[k], position, velocity);
                    ++num_emplaced_markers;
                }
            }
        }
    }

    // 6. Check that we accounted for all points.
#ifndef NDEBUG
    {
        const int ierr =
            MPI_Allreduce(MPI_IN_PLACE, &num_emplaced_markers, 1, MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);
        TBOX_ASSERT(num_emplaced_markers == new_indices.size());
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
    if (rank == 0)
    {
        for (unsigned int k = 0; k < d_markers_outside_domain.size(); ++k)
        {
            const auto index = std::get<0>(d_markers_outside_domain[k]);
            TBOX_ASSERT(index < long(getNumberOfMarkers()));
            marker_check[index] += 1;
        }
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
#endif
    IBTK_TIMER_STOP(t_prune_and_redistribute);
}
} // namespace IBTK
