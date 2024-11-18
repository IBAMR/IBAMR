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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_CFGiesekusStrategy
#define included_CFGiesekusStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CFStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFGiesekusStrategy is a concrete CFStrategy that computes the relaxation function for the
 * Giesekus fluid model.
 *
 * The input database is searched for the following parameters:
 * <ul>
 *   <li><code>relaxation_time</code>: Relaxation time of fluid</li>
 *   <li><code>alpha</code>: Nonlinear parameter</li>
 *   <li><code>viscosity</code>: Polymeric contribution to the viscosity</li>
 * </ul>
 */
class CFGiesekusStrategy : public CFStrategy
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CFGiesekusStrategy(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    void computeRelaxation(int R_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > R_var,
                           int C_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > C_var,
                           TensorEvolutionType evolve_type,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                           double data_time) override;

    void computeStress(int sig_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > sig_var,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                       double data_time) override;

private:
    double d_alpha = std::numeric_limits<double>::quiet_NaN(),
           d_relaxation_time = std::numeric_limits<double>::quiet_NaN(),
           d_viscosity = std::numeric_limits<double>::quiet_NaN();
};

} // Namespace IBAMR
#endif
