// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/LMesh.h"

#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace IBTK
{
class LNode;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMesh::LMesh(std::string object_name, std::vector<LNode*> local_nodes, std::vector<LNode*> ghost_nodes)
    : d_object_name(std::move(object_name)),
      d_local_nodes(std::move(local_nodes)),
      d_ghost_nodes(std::move(ghost_nodes))
{
    // intentionally blank
    return;
} // LMesh

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
