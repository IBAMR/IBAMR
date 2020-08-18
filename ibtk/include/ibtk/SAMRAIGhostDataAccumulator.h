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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SAMRAIGhostDataAccumulator
#define included_IBTK_SAMRAIGhostDataAccumulator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <petscvec.h>

#include <BasePatchHierarchy.h>
#include <RefineSchedule.h>
#include <Variable.h>

#include <vector>

namespace IBTK
{
/*!
 * \brief Class that can accumulate data summed into ghost regions on a patch
 * hierarchy into their correct locations.
 */
class SAMRAIGhostDataAccumulator
{
public:
    /*!
     * Constructor.
     */
    SAMRAIGhostDataAccumulator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > patch_hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                               const SAMRAI::hier::IntVector<NDIM> gcw,
                               const int coarsest_ln,
                               const int finest_ln);

    /*!
     * Accumulate data by summing values in ghost positions into the entry on
     * the owning processor associated with the same degree of freedom.
     */
    void accumulateGhostData(const int idx);

    /*!
     * Destructor.
     */
    ~SAMRAIGhostDataAccumulator();

protected:
    /*!
     * Pointer to the patch hierarchy under consideration.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > d_hierarchy;

    /*!
     * Pointer to the variable whose data layout we copied.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_var;

    /*!
     * Ghost cell width.
     */
    const SAMRAI::hier::IntVector<NDIM> d_gcw;

    /*!
     * Coarsest level of the patch hierarchy on which we work.
     */
    const int d_coarsest_ln = -1;

    /*!
     * Finest level of the patch hierarchy on which we work.
     */
    const int d_finest_ln = -1;

    /*!
     * Boolean indicating whether or not we have cell-centered data.
     */
    bool d_cc_data = true;

    /*!
     * Index into d_hierarchy that contains the dof numbering.
     */
    int d_global_dof_idx = IBTK::invalid_index;

    /*!
     * Index into d_hierarchy that contains the local (i.e., the indices
     * directly into the buffer owned by each Vec in d_vecs) dof numbering.
     */
    int d_local_dof_idx = IBTK::invalid_index;

    /*!
     * PETSc Vec objects storing the global ordering on each level.
     */
    std::vector<Vec> d_vecs;
};
} // namespace IBTK

#endif
