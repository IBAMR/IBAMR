// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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


#include <ibamr/SurfaceTensionForceFunction.h>

#include <BasePatchLevel.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace IBAMR
{

SurfaceTensionForceFunction::SurfaceTensionForceFunction(
    const std::string& object_name,
    Pointer<Database> input_db,
    const AdvDiffHierarchyIntegrator* adv_diff_solver,
    const Pointer<Variable<NDIM>> level_set_var)
    : CartGridFunction(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(level_set_var)
{
    // Intentionally blank.
}

} // namespace IBAMR
