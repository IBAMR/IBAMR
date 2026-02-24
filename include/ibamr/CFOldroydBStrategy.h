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

#ifndef included_CFOldroydBStrategy
#define included_CFOldroydBStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CFStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchGeometry.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"

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
 * \brief Class CFOldroydBStrategy is a concrete CFStrategy that computes the relaxation function for an
 * Oldroyd-B fluid model.
 *
 * The input database is searched for the following parameters:
 * <ul>
 *   <li><code>relaxation_time</code>: Relaxation time of fluid</li>
 *   <li><code>viscosity</code>: Polymeric contribution to the viscosity</li>
 * </ul>
 */
class CFOldroydBStrategy : public CFStrategy
{
public:
    /*!
     * \brief This constructor reads in the parameters for the model from the input database.
     */
    CFOldroydBStrategy(std::string object_name, SAMRAIPointer<SAMRAIDatabase> input_db);

    void computeRelaxation(int R_idx,
                           SAMRAIPointer<SAMRAICellVariable<double> > R_var,
                           int C_idx,
                           SAMRAIPointer<SAMRAICellVariable<double> > C_var,
                           TensorEvolutionType evolve_type,
                           SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                           double data_time) override;

    void computeStress(int sig_idx,
                       SAMRAIPointer<SAMRAICellVariable<double> > sig_var,
                       SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                       double data_time) override;

private:
    double d_relaxation_time = std::numeric_limits<double>::quiet_NaN();
    double d_viscosity = std::numeric_limits<double>::quiet_NaN();
};

} // Namespace IBAMR
#endif
