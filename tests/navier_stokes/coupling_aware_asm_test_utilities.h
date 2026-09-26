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

#ifndef included_IBAMR_tests_coupling_aware_asm_test_utilities
#define included_IBAMR_tests_coupling_aware_asm_test_utilities

#include <ibamr/config.h>

#include <ibtk/IBTK_CHKERRQ.h>

#include <petscis.h>

#include <array>
#include <map>
#include <set>
#include <vector>

using CACell = std::array<int, NDIM>;
using CAFields = std::map<CACell, std::array<int, NDIM + 1>>;

/*! \brief Wrap a logical cell into the periodic test domain. */
inline CACell
ca_wrap(CACell cell, const int n)
{
    for (int& index : cell)
    {
        index = (index + n) % n;
    }
    return cell;
}

/*! \brief Return the pressure and incident velocities of one logical cell. */
inline std::set<int>
ca_cell_stencil(const CAFields& fields, const CACell& cell, const int n)
{
    const std::array<int, NDIM + 1>& local = fields.at(ca_wrap(cell, n));
    std::set<int> result(local.begin(), local.end());
    for (int axis = 0; axis < NDIM; ++axis)
    {
        CACell upper = cell;
        ++upper[axis];
        result.insert(fields.at(ca_wrap(upper, n))[axis]);
    }
    return result;
}

/*! \brief Return the two complete cells incident to a logical face. */
inline std::set<int>
ca_face_stencil(const CAFields& fields, CACell cell, const int axis, const int n)
{
    std::set<int> result = ca_cell_stencil(fields, cell, n);
    --cell[axis];
    const std::set<int> lower = ca_cell_stencil(fields, cell, n);
    result.insert(lower.begin(), lower.end());
    return result;
}

/*! \brief Read the existing PETSc subdomain query into ordered sets. */
inline std::vector<std::set<int>>
ca_read_sets(const std::vector<IS>& indices)
{
    std::vector<std::set<int>> result;
    for (IS is : indices)
    {
        PetscInt count = 0;
        const PetscInt* dofs = nullptr;
        int ierr = ISGetLocalSize(is, &count);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(is, &dofs);
        IBTK_CHKERRQ(ierr);
        result.emplace_back(dofs, dofs + count);
        ierr = ISRestoreIndices(is, &dofs);
        IBTK_CHKERRQ(ierr);
    }
    return result;
}

#endif
