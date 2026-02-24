// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
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

#ifndef included_CFUpperConvectiveOperator
#define included_CFUpperConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/CFStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIFaceData.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableContext.h"
#include "SAMRAIVariableDatabase.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Class CFUpperConvectiveOperator is a concrete ConvectiveOperator that implements the upper convective
 * derivative. This uses second order finite differences to compute the velocity gradients. The transport component can
 * be chosen from any current convective operators. If the advected quantity allows for a square root or logarithmic
 * decomposition, this class can advect the symmetric square root or logarithm of the tensor. Note that this class
 * requires the registration of a source funtion before it can be applied.
 */

class CFUpperConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * Constructor that takes in a string for the convective operator. The convective operator type must be one of
     * "DEFAULT", "CENTERED", "CUI", "PPM", or "WAVE_PROP"
     */
    CFUpperConvectiveOperator(const std::string& object_name,
                              SAMRAIPointer<SAMRAICellVariable<double> > Q_var,
                              SAMRAIPointer<SAMRAIDatabase> input_db,
                              const std::string& convective_op_type,
                              ConvectiveDifferencingType difference_type,
                              const std::vector<SAMRAIRobinBcCoefStrategy*>& Q_bc_coefs,
                              const std::vector<SAMRAIRobinBcCoefStrategy*>& u_bc_coefs);

    /*!
     * Constructor that takes in a convective operator
     */
    CFUpperConvectiveOperator(const std::string& object_name,
                              SAMRAIPointer<SAMRAICellVariable<double> > Q_var,
                              SAMRAIPointer<SAMRAIDatabase> input_db,
                              SAMRAIPointer<ConvectiveOperator> convective_op,
                              ConvectiveDifferencingType difference_type,
                              const std::vector<SAMRAIRobinBcCoefStrategy*>& Q_bc_coefs,
                              const std::vector<SAMRAIRobinBcCoefStrategy*>& u_bc_coefs);

    CFUpperConvectiveOperator() = delete;

    CFUpperConvectiveOperator(const CFUpperConvectiveOperator& from) = delete;

    CFUpperConvectiveOperator& operator=(const CFUpperConvectiveOperator& that) = delete;

    // Destructor
    ~CFUpperConvectiveOperator();

    void applyConvectiveOperator(int Q_idx, int Y_idx) override;

    void initializeOperatorState(const SAMRAISAMRAIVectorReal<double>& in,
                                 const SAMRAISAMRAIVectorReal<double>& out) override;

    void deallocateOperatorState() override;

    /*!
     * \brief: Registers a source function with the convective operator. Note that this source function is given the
     * evolved version of the tensor, and therefore must first be converted from the square root or logarithm to the
     * tensor.
     */
    void registerCFStrategy(SAMRAIPointer<IBAMR::CFStrategy> cf_strategy);

private:
    void commonConstructor(SAMRAIPointer<SAMRAIDatabase> input_db);
    // Hierarchy configuration.
    SAMRAIPointer<SAMRAIPatchHierarchy> d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Scratch data.
    SAMRAIPointer<SAMRAICellVariable<double> > d_Q_var;
    SAMRAIPointer<SAMRAISideVariable<double> > d_u_adv_var;
    int d_u_scratch_idx = IBTK::invalid_index;

    // Source function data.
    SAMRAIPointer<IBAMR::CFStrategy> d_cf_strategy;
    int d_s_idx = IBTK::invalid_index;

    // Convective Operator
    SAMRAIPointer<IBAMR::ConvectiveOperator> d_convec_oper;
    int d_Q_convec_idx = IBTK::invalid_index;
    const std::vector<SAMRAIRobinBcCoefStrategy*> d_Q_bc_coefs;
    const std::vector<SAMRAIRobinBcCoefStrategy*> d_u_bc_coefs;
    TensorEvolutionType d_evolve_type = STANDARD;
    std::string d_interp_type = "LINEAR";
};
} // namespace IBAMR

#endif
