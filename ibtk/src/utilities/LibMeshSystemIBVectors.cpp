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

#include <ibtk/LibMeshSystemIBVectors.h>

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

LibMeshSystemIBVectors::LibMeshSystemIBVectors(const std::vector<IBTK::FEDataManager*>& fe_data_managers,
                                               std::string system_name)
    : LibMeshSystemVectors(std::vector<bool>(fe_data_managers.size(), true), system_name),
      d_fe_data_managers(fe_data_managers)
{
    for (IBTK::FEDataManager* fe_data_manager : d_fe_data_managers)
    {
        libMesh::EquationSystems* es = fe_data_manager->getEquationSystems();
        d_systems.push_back(&es->get_system(d_system_name));
    }
}

LibMeshSystemIBVectors::LibMeshSystemIBVectors(const std::vector<IBTK::FEDataManager*>& fe_data_managers,
                                               const std::vector<bool>& part_mask,
                                               std::string system_name)
    : LibMeshSystemVectors(part_mask, system_name), d_fe_data_managers(fe_data_managers)
{
    TBOX_ASSERT(part_mask.size() == fe_data_managers.size());
    for (unsigned int part = 0; part < d_fe_data_managers.size(); ++part)
    {
        libMesh::EquationSystems* es = d_fe_data_managers[part]->getEquationSystems();
        if (d_part_mask[part])
            d_systems.push_back(&es->get_system(d_system_name));
        else
            d_systems.push_back(nullptr);
    }
}

libMesh::PetscVector<double>&
LibMeshSystemIBVectors::getIBGhosted(const std::string& vec_name, const unsigned int part)
{
    TBOX_ASSERT(part < d_systems.size());
    TBOX_ASSERT(d_part_mask[part]);
    return *maybeAddIBGhosted(vec_name)[part];
}

std::vector<libMesh::PetscVector<double>*>
LibMeshSystemIBVectors::getIBGhosted(const std::string& vec_name)
{
    std::vector<std::unique_ptr<libMesh::PetscVector<double> > >& stored_vectors = maybeAddIBGhosted(vec_name);
    std::vector<libMesh::PetscVector<double>*> result;
    for (std::unique_ptr<libMesh::PetscVector<double> >& ptr : stored_vectors) result.push_back(ptr.get());
    return result;
}

void
LibMeshSystemIBVectors::reinit()
{
    d_ib_ghosted_vectors.clear();
}

/////////////////////////////// PROTECTED ////////////////////////////////////

std::vector<std::unique_ptr<libMesh::PetscVector<double> > >&
LibMeshSystemIBVectors::maybeAddIBGhosted(const std::string& vec_name)
{
    TBOX_ASSERT(d_fe_data_managers.size() == d_systems.size());
    std::vector<std::unique_ptr<libMesh::PetscVector<double> > >& stored_vectors = d_ib_ghosted_vectors[vec_name];
    if (stored_vectors.empty())
        for (unsigned int part = 0; part < d_systems.size(); ++part)
            if (d_part_mask[part])
                stored_vectors.emplace_back(d_fe_data_managers[part]->buildIBGhostedVector(d_system_name));
    return stored_vectors;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
