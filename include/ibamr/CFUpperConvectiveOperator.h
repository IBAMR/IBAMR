// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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
#include "ibamr/CFRelaxationOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

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
    // Constructor
    CFUpperConvectiveOperator(const std::string& object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              const std::string& difference_form,
                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);

    CFUpperConvectiveOperator() = delete;

    CFUpperConvectiveOperator(const CFUpperConvectiveOperator& from) = delete;

    CFUpperConvectiveOperator& operator=(const CFUpperConvectiveOperator& that) = delete;

    // Destructor
    ~CFUpperConvectiveOperator();

    void applyConvectiveOperator(int Q_idx, int Y_idx) override;

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    void deallocateOperatorState() override;

    /*!
     * \brief: Registers a source function with the convective operator. Note that this source function is given the
     * evolved version of the tensor, and therefore must first be converted from the square root or logarithm to the
     * tensor.
     */
    void registerSourceFunction(SAMRAI::tbox::Pointer<IBAMR::CFRelaxationOperator> source_fcn);

private:
    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_u_adv_var;
    int d_u_scratch_idx = IBTK::invalid_index;

    // Source function data.
    SAMRAI::tbox::Pointer<IBAMR::CFRelaxationOperator> d_s_fcn;
    int d_s_idx = IBTK::invalid_index;

    // Convective Operator
    std::string d_difference_form;
    SAMRAI::tbox::Pointer<IBAMR::ConvectiveOperator> d_convec_oper;
    int d_Q_convec_idx = IBTK::invalid_index;
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_Q_bc_coefs;
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;
    TensorEvolutionType d_evolve_type = STANDARD;
    std::string d_interp_type = "LINEAR";
};
} // namespace IBAMR

#endif
