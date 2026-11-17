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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/ViscosityModelLib.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <cmath>
#include <string>
#include <utility>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

CarreauYasudaViscosity::CarreauYasudaViscosity(std::string object_name, Pointer<Database> input_db)
    : d_object_name(std::move(object_name))
{
    const auto get_parameter = [&](const std::string& key, const bool positive)
    {
        const double value = input_db->getDouble(key);
        if (!(positive ? value > 0.0 : value >= 0.0))
        {
            TBOX_ERROR(d_object_name << "::CarreauYasudaViscosity():\n"
                                     << "  invalid " << key << " = " << value << "\n"
                                     << "  " << key << " must be " << (positive ? "positive" : "nonnegative") << "\n");
        }
        return value;
    };
    d_zero_shear_viscosity = get_parameter("zero_shear_viscosity", true);
    d_infinite_shear_viscosity = get_parameter("infinite_shear_viscosity", false);
    d_relaxation_time = get_parameter("relaxation_time", false);
    d_power_law_index = get_parameter("power_law_index", true);
    d_yasuda_exponent = get_parameter("yasuda_exponent", true);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_EE_var = new CellVariable<NDIM, double>(d_object_name + "::EE", NDIM * (NDIM + 1) / 2);
    d_EE_idx = var_db->registerVariableAndContext(
        d_EE_var, var_db->getContext(d_object_name + "::SCRATCH"), IntVector<NDIM>(0));
    return;
} // CarreauYasudaViscosity

int
CarreauYasudaViscosity::getVelocityGhostWidth() const
{
    return 1;
} // getVelocityGhostWidth

void
CarreauYasudaViscosity::setDataOnPatchHierarchy(const int mu_idx,
                                                const int U_idx,
                                                HierarchyMathOps& hier_math_ops,
                                                const double time)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = hier_math_ops.getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_EE_idx, time);
    }

    Pointer<Variable<NDIM>> U_var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(U_idx, U_var);
    const Pointer<SideVariable<NDIM, double>> U_sc_var = U_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(U_sc_var);
#endif
    hier_math_ops.strain_rate(d_EE_idx, d_EE_var, U_idx, U_sc_var, /*src_ghost_fill*/ nullptr, time);

    const double exponent = (d_power_law_index - 1.0) / d_yasuda_exponent;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
            Pointer<CellData<NDIM, double>> EE_data = patch->getPatchData(d_EE_idx);
            for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                // E : E counts each off-diagonal component twice.
                double EE_norm_sq = 0.0;
                for (int d = 0; d < NDIM * (NDIM + 1) / 2; ++d)
                {
                    EE_norm_sq += (d < NDIM ? 1.0 : 2.0) * (*EE_data)(i, d) * (*EE_data)(i, d);
                }
                const double shear_rate = std::sqrt(2.0 * EE_norm_sq);
                (*mu_data)(i) =
                    d_infinite_shear_viscosity +
                    (d_zero_shear_viscosity - d_infinite_shear_viscosity) *
                        std::pow(1.0 + std::pow(d_relaxation_time * shear_rate, d_yasuda_exponent), exponent);
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_EE_idx);
    }
    return;
} // setDataOnPatchHierarchy

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
