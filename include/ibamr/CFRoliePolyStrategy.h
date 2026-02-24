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

#ifndef included_CFRoliePolyStrategy
#define included_CFRoliePolyStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/CFStrategy.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIDatabase.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchGeometry.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIUtilities.h>
#include <SAMRAIVariable.h>

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
    CFRoliePolyStrategy(const std::string& object_name, SAMRAIPointer<SAMRAIDatabase> input_db);

    void computeRelaxation(int R_idx,
                           SAMRAIPointer<SAMRAICellVariable<double>> R_var,
                           int C_idx,
                           SAMRAIPointer<SAMRAICellVariable<double>> C_var,
                           TensorEvolutionType evolve_type,
                           SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                           double data_time) override;

    void computeStress(int sig_idx,
                       SAMRAIPointer<SAMRAICellVariable<double>> sig_var,
                       SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                       double data_time) override;

private:
    double d_lambda_d = std::numeric_limits<double>::quiet_NaN(), d_lambda_R = std::numeric_limits<double>::quiet_NaN(),
           d_beta = std::numeric_limits<double>::quiet_NaN(), d_delta = std::numeric_limits<double>::quiet_NaN(),
           d_viscosity = std::numeric_limits<double>::quiet_NaN();
};

} // Namespace IBAMR
#endif
