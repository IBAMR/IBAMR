// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/CFGiesekusStrategy.h"

#include "ibtk/ibtk_utilities.h"

#include "CellData.h"
#include "CellIterator.h"
#include "Patch.h"
#include "tbox/Database.h"

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

// Namespace
namespace IBAMR
{
CFGiesekusStrategy::CFGiesekusStrategy(std::string object_name, Pointer<Database> input_db)
    : CFStrategy(std::move(object_name))
{
    d_relaxation_time = input_db->getDouble("relaxation_time");
    d_viscosity = input_db->getDouble("viscosity");
    d_alpha = input_db->getDouble("alpha");
    return;
} // Constructor

void
CFGiesekusStrategy::computeStress(int sig_idx,
                                  Pointer<CellVariable<NDIM, double> > /*sig_var*/,
                                  Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  double /*data_time*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > sig_data = patch->getPatchData(sig_idx);

            for (CellIterator<NDIM> ci(sig_data->getGhostBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
#if (NDIM == 2)
                (*sig_data)(idx, 0) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 0) - 1.0);
                (*sig_data)(idx, 1) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 1) - 1.0);
                (*sig_data)(idx, 2) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 2));
#endif
#if (NDIM == 3)
                (*sig_data)(idx, 0) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 0) - 1.0);
                (*sig_data)(idx, 1) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 1) - 1.0);
                (*sig_data)(idx, 2) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 2) - 1.0);
                (*sig_data)(idx, 3) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 3));
                (*sig_data)(idx, 4) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 4));
                (*sig_data)(idx, 5) = d_viscosity / d_relaxation_time * ((*sig_data)(idx, 5));
#endif
            }
        }
    }
}

void
CFGiesekusStrategy::computeRelaxation(const int R_idx,
                                      Pointer<CellVariable<NDIM, double> > /*R_var*/,
                                      int C_idx,
                                      Pointer<CellVariable<NDIM, double> > /*C_var*/,
                                      TensorEvolutionType evolve_type,
                                      Pointer<PatchHierarchy<NDIM> > hierarchy,
                                      double /*data_time*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const double l_inv = 1.0 / d_relaxation_time;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<CellData<NDIM, double> > R_data = patch->getPatchData(R_idx);
            Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(C_idx);

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                MatrixNd mat = convert_to_conformation_tensor(*C_data, idx, evolve_type);
                MatrixNd eye = MatrixNd::Identity();
                mat = -l_inv * ((mat - eye) + d_alpha * (mat - eye) * (mat - eye));
#if (NDIM == 2)
                (*R_data)(idx, 0) = mat(0, 0);
                (*R_data)(idx, 1) = mat(1, 1);
                (*R_data)(idx, 2) = mat(0, 1);
#endif
#if (NDIM == 3)
                (*R_data)(idx, 0) = mat(0, 0);
                (*R_data)(idx, 1) = mat(1, 1);
                (*R_data)(idx, 2) = mat(2, 2);
                (*R_data)(idx, 3) = mat(1, 2);
                (*R_data)(idx, 4) = mat(0, 2);
                (*R_data)(idx, 5) = mat(0, 1);
#endif
            }
        }
    }
} // setDataOnPatch

} // namespace IBAMR
