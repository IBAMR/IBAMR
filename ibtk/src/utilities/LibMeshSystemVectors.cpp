// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <ibtk/LibMeshSystemVectors.h>

#include <tbox/Utilities.h>

#include <libmesh/equation_systems.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/system.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// CLASS DEFINITION /////////////////////////////

LibMeshSystemVectors::LibMeshSystemVectors(const std::vector<libMesh::EquationSystems*>& equation_systems,
                                           std::string system_name)
    : LibMeshSystemVectors(equation_systems, std::vector<bool>(equation_systems.size(), true), system_name)
{
}

LibMeshSystemVectors::LibMeshSystemVectors(const std::vector<libMesh::EquationSystems*>& equation_systems,
                                           const std::vector<bool>& part_mask,
                                           std::string system_name)
    : d_part_mask(part_mask), d_system_name(std::move(system_name))
{
    TBOX_ASSERT(part_mask.size() == equation_systems.size());
    for (unsigned int part = 0; part < equation_systems.size(); ++part)
    {
        if (d_part_mask[part])
            d_systems.push_back(&equation_systems[part]->get_system(d_system_name));
        else
            d_systems.push_back(nullptr);
    }
}

libMesh::PetscVector<double>&
LibMeshSystemVectors::get(const std::string& vec_name, const unsigned int part)
{
    TBOX_ASSERT(part < d_systems.size());
    TBOX_ASSERT(d_part_mask[part]);
    if (!vec_stored_in_map(vec_name))
    {
        if (vec_name == "solution") return dynamic_cast<libMesh::PetscVector<double>&>(*d_systems[part]->solution);
        if (vec_name == "current" || vec_name == "current_local_solution")
            return dynamic_cast<libMesh::PetscVector<double>&>(*d_systems[part]->current_local_solution);
        else
            TBOX_ASSERT(false);
    }
    maybeAdd(vec_name);
    return dynamic_cast<libMesh::PetscVector<double>&>(d_systems[part]->get_vector(vec_name));
}

std::vector<libMesh::PetscVector<double>*>
LibMeshSystemVectors::get(const std::string& vec_name)
{
    std::vector<libMesh::PetscVector<double>*> result;
    for (unsigned int part = 0; part < d_systems.size(); ++part)
    {
        if (d_part_mask[part])
            result.push_back(&get(vec_name, part));
        else
            result.push_back(nullptr);
    }

    return result;
}

void
LibMeshSystemVectors::reinit()
{
}

void
LibMeshSystemVectors::copy(const std::string& source, const std::vector<std::string>& dests)
{
    // The source vectors must already exist
    TBOX_ASSERT(!vec_stored_in_map(source) || d_systems[0]->request_vector(source));
    for (const std::string& dest : dests)
    {
        maybeAdd(dest);
        for (unsigned int part = 0; part < d_systems.size(); ++part)
        {
            if (d_part_mask[part])
            {
                auto& left = get(dest, part);
                auto& right = get(source, part);
                left = right;
            }
        }
    }
}

void
LibMeshSystemVectors::zero(const std::string& vec_name)
{
    // The source vectors must already exist
    bool vec_exists = false;
    for (unsigned int part = 0; part < d_systems.size(); ++part)
        if (d_part_mask[part] && d_systems[part]->request_vector(vec_name))
        {
            vec_exists = true;
            break;
        }

    TBOX_ASSERT(vec_exists);
    std::vector<libMesh::PetscVector<double>*> vecs = get(vec_name);
    for (unsigned int part = 0; part < d_systems.size(); ++part)
        if (d_part_mask[part]) vecs[part]->zero();
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
LibMeshSystemVectors::maybeAdd(const std::string& vec_name)
{
    bool vec_exists = false;
    for (unsigned int part = 0; part < d_systems.size(); ++part)
        if (d_part_mask[part] && d_systems[part]->request_vector(vec_name))
        {
            vec_exists = true;
            break;
        }
    if (vec_stored_in_map(vec_name) && !vec_exists)
        for (unsigned int part = 0; part < d_systems.size(); ++part)
            if (d_part_mask[part])
                d_systems[part]->add_vector(vec_name, /*projections*/ true, /*type*/ libMesh::GHOSTED);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
