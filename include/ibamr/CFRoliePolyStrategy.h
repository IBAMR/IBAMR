// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#ifndef included_CFRoliePolyStrategy
#define included_CFRoliePolyStrategy

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

#include <cmath>
#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFRoliePolyStrategy is a concrete CFStrategy that computes the relaxation function for the
 * Rolie-Poly fluid model.
 *
 * The input database is searched for the following parameters:
 * <ul>
 *   <li><code>lambda_d</code>: Reptation time</li>
 *   <li><code>lambda_R</code>: Rouse time</li>
 *   <li><code>beta</code>: CCR coefficient</li>
 *   <li><code>delta</code></li>
 *   <li><code>viscosity</code>: Polymeric contribution to the viscosity</li>
 * </ul>
 * We let the reptation time <code>lambda_d</code> be the relaxation time used to convert the conformation tensor to the
 * stress tensor.
 */
class CFRoliePolyStrategy : public CFStrategy
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CFRoliePolyStrategy(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

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
    double d_lambda_d = std::numeric_limits<double>::quiet_NaN(), d_lambda_R = std::numeric_limits<double>::quiet_NaN(),
           d_beta = std::numeric_limits<double>::quiet_NaN(), d_delta = std::numeric_limits<double>::quiet_NaN(),
           d_viscosity = std::numeric_limits<double>::quiet_NaN();
};

} // Namespace IBAMR
#endif
