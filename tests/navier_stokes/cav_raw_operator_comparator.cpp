// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "cav_raw_operator_comparator.h"

#include <ibamr/app_namespaces.h>

namespace IBAMR
{
namespace TestSupport
{
namespace
{
using EntryMap = std::map<std::pair<PetscInt, PetscInt>, double>;

std::ofstream
open_output(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("unable to open output file " + filename);
    out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    return out;
}

std::ifstream
open_input(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("unable to open input file " + filename);
    return in;
}

const CAVRawDofRecord*
find_dof_record(const CAVRawOperatorBundle& bundle, const PetscInt dof)
{
    const auto it = std::find_if(
        bundle.dofs.begin(), bundle.dofs.end(), [&](const CAVRawDofRecord& record) { return record.dof == dof; });
    return it == bundle.dofs.end() ? nullptr : &*it;
}

bool
valid_dof_records(const CAVRawOperatorBundle& bundle)
{
    if (bundle.dimension != NDIM || bundle.nrows != bundle.ncols || bundle.nrows < 0 ||
        bundle.dofs.size() != static_cast<std::size_t>(bundle.nrows))
    {
        return false;
    }
    std::vector<bool> found(static_cast<std::size_t>(bundle.nrows), false);
    for (const CAVRawDofRecord& record : bundle.dofs)
    {
        if (record.dof < 0 || record.dof >= bundle.nrows || found[static_cast<std::size_t>(record.dof)]) return false;
        found[static_cast<std::size_t>(record.dof)] = true;
    }
    return std::all_of(found.begin(), found.end(), [](const bool value) { return value; });
}

EntryMap
aggregate_entries(const std::vector<CAVRawMatrixEntry>& entries)
{
    EntryMap result;
    for (const CAVRawMatrixEntry& entry : entries)
    {
        result[{ entry.row, entry.column }] += entry.value;
    }
    for (auto it = result.begin(); it != result.end();)
    {
        if (it->second == 0.0)
            it = result.erase(it);
        else
            ++it;
    }
    return result;
}

std::pair<double, double>
comparison_errors(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size())
        return { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
    double max_abs_error = 0.0;
    double lhs_norm = 0.0;
    double rhs_norm = 0.0;
    for (std::size_t k = 0; k < lhs.size(); ++k)
    {
        max_abs_error = std::max(max_abs_error, std::abs(lhs[k] - rhs[k]));
        lhs_norm = std::max(lhs_norm, std::abs(lhs[k]));
        rhs_norm = std::max(rhs_norm, std::abs(rhs[k]));
    }
    return { max_abs_error, max_abs_error / std::max({ 1.0, lhs_norm, rhs_norm }) };
}

std::vector<PetscInt>
build_candidate_to_reference_map(const CAVRawOperatorBundle& candidate,
                                 const CAVRawOperatorBundle& reference,
                                 const double coordinate_tolerance)
{
    if (!valid_dof_records(candidate) || !valid_dof_records(reference) || candidate.nrows != reference.nrows) return {};

    std::vector<PetscInt> candidate_to_reference(static_cast<std::size_t>(candidate.nrows), -1);
    std::set<PetscInt> used_reference_dofs;
    for (const CAVRawDofRecord& candidate_record : candidate.dofs)
    {
        std::vector<PetscInt> matches;
        for (const CAVRawDofRecord& reference_record : reference.dofs)
        {
            if (candidate_record.kind != reference_record.kind || candidate_record.axis != reference_record.axis ||
                used_reference_dofs.count(reference_record.dof) != 0)
            {
                continue;
            }
            bool coordinates_match = true;
            for (int d = 0; d < NDIM; ++d)
            {
                coordinates_match = coordinates_match && std::abs(candidate_record.position[d] -
                                                                  reference_record.position[d]) <= coordinate_tolerance;
            }
            if (coordinates_match) matches.push_back(reference_record.dof);
        }
        if (matches.size() != 1) return {};
        candidate_to_reference[static_cast<std::size_t>(candidate_record.dof)] = matches.front();
        used_reference_dofs.insert(matches.front());
    }
    if (used_reference_dofs.size() != reference.dofs.size()) return {};
    return candidate_to_reference;
}

std::vector<double>
map_vector(const CAVRawOperatorBundle& candidate,
           const std::vector<double>& values,
           const std::vector<PetscInt>& candidate_to_reference,
           const bool equation_values,
           const double pressure_equation_scale)
{
    if (values.size() != static_cast<std::size_t>(candidate.nrows)) return {};
    std::vector<double> mapped(values.size(), 0.0);
    for (PetscInt candidate_dof = 0; candidate_dof < candidate.nrows; ++candidate_dof)
    {
        const CAVRawDofRecord* const record = find_dof_record(candidate, candidate_dof);
        if (!record) return {};
        const double scale = equation_values && record->kind == CAVRawDofKind::PRESSURE ? pressure_equation_scale : 1.0;
        mapped[static_cast<std::size_t>(candidate_to_reference[static_cast<std::size_t>(candidate_dof)])] =
            scale * values[static_cast<std::size_t>(candidate_dof)];
    }
    return mapped;
}

void
subtract_pressure_mean(std::vector<double>& values, const CAVRawOperatorBundle& reference)
{
    double sum = 0.0;
    std::size_t count = 0;
    for (const CAVRawDofRecord& record : reference.dofs)
    {
        if (record.kind == CAVRawDofKind::PRESSURE)
        {
            sum += values[static_cast<std::size_t>(record.dof)];
            ++count;
        }
    }
    if (count == 0) return;
    const double mean = sum / static_cast<double>(count);
    for (const CAVRawDofRecord& record : reference.dofs)
    {
        if (record.kind == CAVRawDofKind::PRESSURE) values[static_cast<std::size_t>(record.dof)] -= mean;
    }
}

CAVRawOperatorBundle
make_control_reference()
{
    CAVRawOperatorBundle bundle;
    bundle.dimension = NDIM;
    bundle.nrows = 4;
    bundle.ncols = 4;
    bundle.dofs.resize(4);

    bundle.dofs[0].dof = 0;
    bundle.dofs[0].kind = CAVRawDofKind::VELOCITY;
    bundle.dofs[0].axis = 0;
    bundle.dofs[0].position[0] = 0.0;
    bundle.dofs[0].position[1] = 0.5;

    bundle.dofs[1].dof = 1;
    bundle.dofs[1].kind = CAVRawDofKind::VELOCITY;
    bundle.dofs[1].axis = 1;
    bundle.dofs[1].position[0] = 0.5;
    bundle.dofs[1].position[1] = 0.0;

    bundle.dofs[2].dof = 2;
    bundle.dofs[2].kind = CAVRawDofKind::PRESSURE;
    bundle.dofs[2].axis = -1;
    bundle.dofs[2].position[0] = 0.5;
    bundle.dofs[2].position[1] = 0.5;

    bundle.dofs[3].dof = 3;
    bundle.dofs[3].kind = CAVRawDofKind::PRESSURE;
    bundle.dofs[3].axis = -1;
    bundle.dofs[3].position[0] = 1.5;
    bundle.dofs[3].position[1] = 0.5;

    bundle.matrix_entries = { { 0, 0, 2.0 }, { 0, 2, 1.0 },  { 1, 1, 3.0 },  { 1, 3, -1.0 },
                              { 2, 0, 4.0 }, { 2, 1, -4.0 }, { 3, 0, -2.0 }, { 3, 1, 2.0 } };
    bundle.equation_values = { 1.0, -2.0, 3.0, -4.0 };
    bundle.state_values = { 0.25, -0.5, 2.0, 5.0 };
    bundle.ordered_dofs = { 2, 3, 0, 1 };
    return bundle;
}

CAVRawOperatorBundle
permute_bundle(const CAVRawOperatorBundle& input, const std::vector<PetscInt>& old_to_new)
{
    CAVRawOperatorBundle output = input;
    for (CAVRawDofRecord& record : output.dofs) record.dof = old_to_new[static_cast<std::size_t>(record.dof)];
    std::sort(
        output.dofs.begin(), output.dofs.end(), [](const auto& lhs, const auto& rhs) { return lhs.dof < rhs.dof; });
    for (CAVRawMatrixEntry& entry : output.matrix_entries)
    {
        entry.row = old_to_new[static_cast<std::size_t>(entry.row)];
        entry.column = old_to_new[static_cast<std::size_t>(entry.column)];
    }
    output.equation_values.assign(input.equation_values.size(), 0.0);
    output.state_values.assign(input.state_values.size(), 0.0);
    for (std::size_t old_dof = 0; old_dof < old_to_new.size(); ++old_dof)
    {
        const std::size_t new_dof = static_cast<std::size_t>(old_to_new[old_dof]);
        output.equation_values[new_dof] = input.equation_values[old_dof];
        output.state_values[new_dof] = input.state_values[old_dof];
    }
    for (PetscInt& dof : output.ordered_dofs) dof = old_to_new[static_cast<std::size_t>(dof)];
    return output;
}

void
scale_pressure_equations(CAVRawOperatorBundle& bundle, const double scale)
{
    for (CAVRawMatrixEntry& entry : bundle.matrix_entries)
    {
        const CAVRawDofRecord* const row_record = find_dof_record(bundle, entry.row);
        if (row_record && row_record->kind == CAVRawDofKind::PRESSURE) entry.value *= scale;
    }
    for (const CAVRawDofRecord& record : bundle.dofs)
    {
        if (record.kind == CAVRawDofKind::PRESSURE)
            bundle.equation_values[static_cast<std::size_t>(record.dof)] *= scale;
    }
}

void
shift_state(CAVRawOperatorBundle& bundle, const CAVRawDofKind kind, const double shift)
{
    for (const CAVRawDofRecord& record : bundle.dofs)
    {
        if (record.kind == kind) bundle.state_values[static_cast<std::size_t>(record.dof)] += shift;
    }
}
} // namespace

CAVRawOperatorBundle
captureCAVRawOperatorBundle(const StaggeredStokesPETScLevelSolver::LiveOperatorStateView& view,
                            Pointer<PatchLevel<NDIM>> level,
                            const int u_dof_index_idx,
                            const int p_dof_index_idx,
                            Vec equation_values,
                            Vec state_values)
{
    // Export directly from the borrowed live matrix and the level's global velocity-pressure numbering. The
    // production matrix is never duplicated or translated into subdomain-local algebraic indices here.
    if (IBTK_MPI::getNodes() != 1) throw std::runtime_error("raw operator export currently requires one MPI rank");
    if (!view.initialized || !view.operator_mat || !view.locally_owned_velocity_dofs ||
        !view.locally_owned_pressure_dofs || !level)
    {
        throw std::runtime_error("raw operator export requires an initialized live operator view");
    }

    CAVRawOperatorBundle bundle;
    int ierr = MatGetSize(view.operator_mat, &bundle.nrows, &bundle.ncols);
    IBTK_CHKERRQ(ierr);
    PetscInt first_row = 0, one_past_row = 0;
    ierr = MatGetOwnershipRange(view.operator_mat, &first_row, &one_past_row);
    IBTK_CHKERRQ(ierr);
    if (first_row != 0 || one_past_row != bundle.nrows)
        throw std::runtime_error("serial raw operator export does not own every matrix row");

    for (PetscInt row = 0; row < bundle.nrows; ++row)
    {
        PetscInt ncols = 0;
        const PetscInt* columns = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(view.operator_mat, row, &ncols, &columns, &values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            bundle.matrix_entries.push_back({ row, columns[k], static_cast<double>(PetscRealPart(values[k])) });
        }
        ierr = MatRestoreRow(view.operator_mat, row, &ncols, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }

    std::map<PetscInt, CAVRawDofRecord> records;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();
        const hier::Index<NDIM>& patch_lower = patch_box.lower();
        const double* const x_lower = patch_geometry->getXLower();
        const double* const dx = patch_geometry->getDx();

        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const SideIndex<NDIM> side_index(b(), axis, SideIndex<NDIM>::Lower);
                const PetscInt dof = (*u_dof_data)(side_index);
                if (dof < 0 || records.count(dof) != 0) continue;
                CAVRawDofRecord record;
                record.dof = dof;
                record.kind = CAVRawDofKind::VELOCITY;
                record.axis = axis;
                for (int d = 0; d < NDIM; ++d)
                {
                    record.index[d] = side_index(d);
                    record.position[d] =
                        x_lower[d] +
                        (static_cast<double>(side_index(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5)) * dx[d];
                }
                records.emplace(dof, record);
            }
        }

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& cell_index = b();
            const PetscInt dof = (*p_dof_data)(cell_index);
            if (dof < 0 || records.count(dof) != 0) continue;
            CAVRawDofRecord record;
            record.dof = dof;
            record.kind = CAVRawDofKind::PRESSURE;
            record.axis = -1;
            for (int d = 0; d < NDIM; ++d)
            {
                record.index[d] = cell_index(d);
                record.position[d] = x_lower[d] + (static_cast<double>(cell_index(d) - patch_lower(d)) + 0.5) * dx[d];
            }
            records.emplace(dof, record);
        }
    }
    for (const auto& record : records) bundle.dofs.push_back(record.second);
    if (!valid_dof_records(bundle)) throw std::runtime_error("raw operator export found an invalid global DOF map");

    const auto copy_vec = [&](Vec vec, std::vector<double>& destination)
    {
        destination.assign(static_cast<std::size_t>(bundle.nrows), 0.0);
        if (!vec) return;
        PetscInt size = 0;
        int vec_ierr = VecGetSize(vec, &size);
        IBTK_CHKERRQ(vec_ierr);
        if (size != bundle.nrows) throw std::runtime_error("raw operator vector size does not match matrix size");
        const PetscScalar* values = nullptr;
        vec_ierr = VecGetArrayRead(vec, &values);
        IBTK_CHKERRQ(vec_ierr);
        for (PetscInt k = 0; k < size; ++k) destination[static_cast<std::size_t>(k)] = PetscRealPart(values[k]);
        vec_ierr = VecRestoreArrayRead(vec, &values);
        IBTK_CHKERRQ(vec_ierr);
    };
    copy_vec(equation_values, bundle.equation_values);
    copy_vec(state_values, bundle.state_values);
    bundle.ordered_dofs.reserve(bundle.dofs.size());
    for (const CAVRawDofRecord& record : bundle.dofs) bundle.ordered_dofs.push_back(record.dof);
    return bundle;
}

void
writeCAVRawOperatorBundle(const CAVRawOperatorBundle& bundle, const std::string& prefix)
{
    {
        std::ofstream out = open_output(prefix + ".meta");
        out << "cav-raw-operator-v1 " << bundle.dimension << " " << bundle.nrows << " " << bundle.ncols << "\n";
    }
    {
        std::ofstream out = open_output(prefix + ".dofs");
        out << bundle.dofs.size() << "\n";
        for (const CAVRawDofRecord& record : bundle.dofs)
        {
            out << record.dof << " " << (record.kind == CAVRawDofKind::VELOCITY ? "V" : "P") << " " << record.axis;
            for (const int index : record.index) out << " " << index;
            for (const double position : record.position) out << " " << position;
            out << "\n";
        }
    }
    {
        std::ofstream out = open_output(prefix + ".mtx");
        out << "%%MatrixMarket matrix coordinate real general\n";
        out << bundle.nrows << " " << bundle.ncols << " " << bundle.matrix_entries.size() << "\n";
        for (const CAVRawMatrixEntry& entry : bundle.matrix_entries)
        {
            out << entry.row + 1 << " " << entry.column + 1 << " " << entry.value << "\n";
        }
    }
    const auto write_vector = [&](const std::string& suffix, const std::vector<double>& values)
    {
        std::ofstream out = open_output(prefix + suffix);
        out << values.size() << "\n";
        for (std::size_t k = 0; k < values.size(); ++k) out << k << " " << values[k] << "\n";
    };
    write_vector(".equation", bundle.equation_values);
    write_vector(".state", bundle.state_values);
    {
        std::ofstream out = open_output(prefix + ".order");
        out << bundle.ordered_dofs.size() << "\n";
        for (const PetscInt dof : bundle.ordered_dofs) out << dof << "\n";
    }
}

CAVRawOperatorBundle
readCAVRawOperatorBundle(const std::string& prefix)
{
    CAVRawOperatorBundle bundle;
    {
        std::ifstream in = open_input(prefix + ".meta");
        std::string schema;
        in >> schema >> bundle.dimension >> bundle.nrows >> bundle.ncols;
        if (!in || schema != "cav-raw-operator-v1" || bundle.dimension != NDIM)
            throw std::runtime_error("invalid raw operator metadata");
    }
    {
        std::ifstream in = open_input(prefix + ".dofs");
        std::size_t count = 0;
        in >> count;
        bundle.dofs.resize(count);
        for (CAVRawDofRecord& record : bundle.dofs)
        {
            std::string kind;
            in >> record.dof >> kind >> record.axis;
            if (kind == "P")
                record.kind = CAVRawDofKind::PRESSURE;
            else if (kind == "V")
                record.kind = CAVRawDofKind::VELOCITY;
            else
                throw std::runtime_error("invalid raw operator DOF kind");
            for (int& index : record.index) in >> index;
            for (double& position : record.position) in >> position;
        }
        if (!in) throw std::runtime_error("invalid raw operator DOF map");
    }
    {
        std::ifstream in = open_input(prefix + ".mtx");
        std::string header;
        std::getline(in, header);
        if (header != "%%MatrixMarket matrix coordinate real general")
            throw std::runtime_error("invalid raw operator MatrixMarket header");
        PetscInt nrows = 0, ncols = 0;
        std::size_t count = 0;
        in >> nrows >> ncols >> count;
        if (nrows != bundle.nrows || ncols != bundle.ncols) throw std::runtime_error("raw operator size mismatch");
        bundle.matrix_entries.resize(count);
        for (CAVRawMatrixEntry& entry : bundle.matrix_entries)
        {
            in >> entry.row >> entry.column >> entry.value;
            --entry.row;
            --entry.column;
        }
        if (!in) throw std::runtime_error("invalid raw operator MatrixMarket entries");
    }
    const auto read_vector = [&](const std::string& suffix, std::vector<double>& values)
    {
        std::ifstream in = open_input(prefix + suffix);
        std::size_t count = 0;
        in >> count;
        values.assign(count, 0.0);
        std::vector<bool> found(count, false);
        for (std::size_t k = 0; k < count; ++k)
        {
            std::size_t index = 0;
            double value = 0.0;
            in >> index >> value;
            if (index >= count || found[index]) throw std::runtime_error("invalid raw operator vector index");
            values[index] = value;
            found[index] = true;
        }
        if (!in || !std::all_of(found.begin(), found.end(), [](const bool value) { return value; }))
            throw std::runtime_error("invalid raw operator vector");
    };
    read_vector(".equation", bundle.equation_values);
    read_vector(".state", bundle.state_values);
    {
        std::ifstream in = open_input(prefix + ".order");
        std::size_t count = 0;
        in >> count;
        bundle.ordered_dofs.resize(count);
        for (PetscInt& dof : bundle.ordered_dofs) in >> dof;
        if (!in) throw std::runtime_error("invalid raw operator order data");
    }
    return bundle;
}

CAVRawComparison
compareCAVRawOperatorBundles(const CAVRawOperatorBundle& candidate,
                             const CAVRawOperatorBundle& reference,
                             const CAVRawMappingSpec& mapping)
{
    CAVRawComparison result;
    if (candidate.equation_values.size() != static_cast<std::size_t>(candidate.nrows) ||
        candidate.state_values.size() != static_cast<std::size_t>(candidate.nrows) ||
        reference.equation_values.size() != static_cast<std::size_t>(reference.nrows) ||
        reference.state_values.size() != static_cast<std::size_t>(reference.nrows))
    {
        result.first_mismatch = "equation or state vector has the wrong size";
        return result;
    }
    const std::vector<PetscInt> candidate_to_reference =
        build_candidate_to_reference_map(candidate, reference, mapping.coordinate_tolerance);
    result.dof_bijection = candidate_to_reference.size() == static_cast<std::size_t>(candidate.nrows);
    if (!result.dof_bijection)
    {
        result.first_mismatch = "DOF mapping is not a bijection";
        return result;
    }

    EntryMap mapped_candidate;
    for (const CAVRawMatrixEntry& entry : candidate.matrix_entries)
    {
        if (entry.row < 0 || entry.row >= candidate.nrows || entry.column < 0 || entry.column >= candidate.ncols)
        {
            result.first_mismatch = "candidate matrix entry is out of range";
            return result;
        }
        const CAVRawDofRecord* const row_record = find_dof_record(candidate, entry.row);
        const double scale =
            row_record && row_record->kind == CAVRawDofKind::PRESSURE ? mapping.candidate_pressure_equation_scale : 1.0;
        const PetscInt mapped_row = candidate_to_reference[static_cast<std::size_t>(entry.row)];
        const PetscInt mapped_column = candidate_to_reference[static_cast<std::size_t>(entry.column)];
        mapped_candidate[{ mapped_row, mapped_column }] += scale * entry.value;
    }
    for (auto it = mapped_candidate.begin(); it != mapped_candidate.end();)
    {
        if (it->second == 0.0)
            it = mapped_candidate.erase(it);
        else
            ++it;
    }
    const EntryMap reference_entries = aggregate_entries(reference.matrix_entries);
    result.matrix_structure = mapped_candidate.size() == reference_entries.size();
    if (result.matrix_structure)
    {
        auto candidate_it = mapped_candidate.begin();
        auto reference_it = reference_entries.begin();
        for (; candidate_it != mapped_candidate.end(); ++candidate_it, ++reference_it)
        {
            if (candidate_it->first != reference_it->first)
            {
                result.matrix_structure = false;
                break;
            }
        }
    }

    std::vector<double> candidate_matrix_values;
    std::vector<double> reference_matrix_values;
    std::set<std::pair<PetscInt, PetscInt>> matrix_coordinates;
    for (const auto& entry : mapped_candidate) matrix_coordinates.insert(entry.first);
    for (const auto& entry : reference_entries) matrix_coordinates.insert(entry.first);
    for (const auto& coordinate : matrix_coordinates)
    {
        const auto candidate_it = mapped_candidate.find(coordinate);
        const auto reference_it = reference_entries.find(coordinate);
        candidate_matrix_values.push_back(candidate_it == mapped_candidate.end() ? 0.0 : candidate_it->second);
        reference_matrix_values.push_back(reference_it == reference_entries.end() ? 0.0 : reference_it->second);
    }
    std::tie(result.matrix_max_abs_error, result.matrix_scaled_error) =
        comparison_errors(candidate_matrix_values, reference_matrix_values);
    result.matrix_values = result.matrix_structure && result.matrix_scaled_error <= mapping.comparison_tolerance;

    const std::vector<double> mapped_equation = map_vector(
        candidate, candidate.equation_values, candidate_to_reference, true, mapping.candidate_pressure_equation_scale);
    std::tie(result.equation_max_abs_error, result.equation_scaled_error) =
        comparison_errors(mapped_equation, reference.equation_values);
    result.equation_values = result.equation_scaled_error <= mapping.comparison_tolerance;

    std::vector<double> mapped_state =
        map_vector(candidate, candidate.state_values, candidate_to_reference, false, 1.0);
    std::vector<double> reference_state = reference.state_values;
    if (mapping.align_pressure_gauge && mapped_state.size() == reference_state.size())
    {
        subtract_pressure_mean(mapped_state, reference);
        subtract_pressure_mean(reference_state, reference);
    }
    std::tie(result.state_max_abs_error, result.state_scaled_error) = comparison_errors(mapped_state, reference_state);
    result.state_values = result.state_scaled_error <= mapping.comparison_tolerance;

    std::vector<PetscInt> mapped_order;
    mapped_order.reserve(candidate.ordered_dofs.size());
    for (const PetscInt dof : candidate.ordered_dofs)
    {
        if (dof < 0 || dof >= candidate.nrows)
        {
            result.first_mismatch = "ordered DOF is out of range";
            return result;
        }
        mapped_order.push_back(candidate_to_reference[static_cast<std::size_t>(dof)]);
    }
    result.ordered_dofs = mapped_order == reference.ordered_dofs;

    result.matched = result.dof_bijection && result.matrix_structure && result.matrix_values &&
                     result.equation_values && result.state_values && result.ordered_dofs;
    if (!result.matched && result.first_mismatch.empty())
    {
        if (!result.matrix_structure)
            result.first_mismatch = "matrix structure differs after mapping";
        else if (!result.matrix_values)
            result.first_mismatch = "matrix values differ after mapping";
        else if (!result.equation_values)
            result.first_mismatch = "equation values differ after row-sign mapping";
        else if (!result.state_values)
            result.first_mismatch = "state values differ after pressure-gauge alignment";
        else if (!result.ordered_dofs)
            result.first_mismatch = "ordered artifact differs after mapping";
    }
    return result;
}

PetscInt
findCAVRawDof(const CAVRawOperatorBundle& bundle,
              const CAVRawDofKind kind,
              const int axis,
              const std::array<int, NDIM>& index)
{
    const auto it = std::find_if(bundle.dofs.begin(),
                                 bundle.dofs.end(),
                                 [&](const CAVRawDofRecord& record)
                                 { return record.kind == kind && record.axis == axis && record.index == index; });
    return it == bundle.dofs.end() ? -1 : it->dof;
}

double
getCAVRawMatrixValue(const CAVRawOperatorBundle& bundle, const PetscInt row, const PetscInt column)
{
    double value = 0.0;
    for (const CAVRawMatrixEntry& entry : bundle.matrix_entries)
    {
        if (entry.row == row && entry.column == column) value += entry.value;
    }
    return value;
}

std::size_t
countCAVRawMatrixRowNonzeros(const CAVRawOperatorBundle& bundle, const PetscInt row)
{
    const EntryMap entries = aggregate_entries(bundle.matrix_entries);
    return static_cast<std::size_t>(
        std::count_if(entries.begin(), entries.end(), [&](const auto& entry) { return entry.first.first == row; }));
}

CAVRawComparatorControlResults
runCAVRawComparatorControls()
{
    CAVRawComparatorControlResults results;
    const CAVRawOperatorBundle reference = make_control_reference();
    CAVRawMappingSpec identity_mapping;
    results.exact_control = compareCAVRawOperatorBundles(reference, reference, identity_mapping).matched;

    const std::vector<PetscInt> permutation = { 2, 0, 3, 1 };
    CAVRawOperatorBundle declared_candidate = permute_bundle(reference, permutation);
    scale_pressure_equations(declared_candidate, -1.0);
    shift_state(declared_candidate, CAVRawDofKind::PRESSURE, 7.0);
    CAVRawMappingSpec declared_mapping;
    declared_mapping.candidate_pressure_equation_scale = -1.0;
    declared_mapping.align_pressure_gauge = true;
    results.declared_permutation_sign_gauge =
        compareCAVRawOperatorBundles(declared_candidate, reference, declared_mapping).matched;

    CAVRawOperatorBundle undeclared_candidate = declared_candidate;
    undeclared_candidate.dofs = reference.dofs;
    results.undeclared_permutation_detected =
        !compareCAVRawOperatorBundles(undeclared_candidate, reference, declared_mapping).matched;

    CAVRawOperatorBundle velocity_sign_candidate = declared_candidate;
    velocity_sign_candidate.equation_values[static_cast<std::size_t>(permutation[0])] *= -1.0;
    results.velocity_sign_detected =
        !compareCAVRawOperatorBundles(velocity_sign_candidate, reference, declared_mapping).matched;

    CAVRawMappingSpec wrong_pressure_sign_mapping = declared_mapping;
    wrong_pressure_sign_mapping.candidate_pressure_equation_scale = 1.0;
    results.pressure_sign_detected =
        !compareCAVRawOperatorBundles(declared_candidate, reference, wrong_pressure_sign_mapping).matched;

    CAVRawOperatorBundle legal_gauge_candidate = permute_bundle(reference, permutation);
    shift_state(legal_gauge_candidate, CAVRawDofKind::PRESSURE, -4.0);
    CAVRawMappingSpec gauge_mapping;
    gauge_mapping.align_pressure_gauge = true;
    results.legal_pressure_gauge =
        compareCAVRawOperatorBundles(legal_gauge_candidate, reference, gauge_mapping).matched;

    CAVRawOperatorBundle illegal_velocity_candidate = legal_gauge_candidate;
    shift_state(illegal_velocity_candidate, CAVRawDofKind::VELOCITY, 1.0);
    results.illegal_velocity_shift_detected =
        !compareCAVRawOperatorBundles(illegal_velocity_candidate, reference, gauge_mapping).matched;

    CAVRawOperatorBundle omitted_dof_candidate = declared_candidate;
    omitted_dof_candidate.dofs.pop_back();
    results.omitted_dof_detected =
        !compareCAVRawOperatorBundles(omitted_dof_candidate, reference, declared_mapping).matched;

    CAVRawOperatorBundle omitted_row_candidate = declared_candidate;
    const PetscInt omitted_row = permutation[0];
    omitted_row_candidate.matrix_entries.erase(std::remove_if(omitted_row_candidate.matrix_entries.begin(),
                                                              omitted_row_candidate.matrix_entries.end(),
                                                              [&](const CAVRawMatrixEntry& entry)
                                                              { return entry.row == omitted_row; }),
                                               omitted_row_candidate.matrix_entries.end());
    results.omitted_row_detected =
        !compareCAVRawOperatorBundles(omitted_row_candidate, reference, declared_mapping).matched;

    CAVRawOperatorBundle ordered_omission_candidate = declared_candidate;
    ordered_omission_candidate.ordered_dofs.pop_back();
    results.ordered_omission_detected =
        !compareCAVRawOperatorBundles(ordered_omission_candidate, reference, declared_mapping).matched;

    CAVRawOperatorBundle reordered_candidate = declared_candidate;
    std::swap(reordered_candidate.ordered_dofs[0], reordered_candidate.ordered_dofs[1]);
    results.ordered_reordering_detected =
        !compareCAVRawOperatorBundles(reordered_candidate, reference, declared_mapping).matched;
    return results;
}

} // namespace TestSupport
} // namespace IBAMR
