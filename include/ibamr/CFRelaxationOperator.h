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

#ifndef included_CFRelaxationOperator
#define included_CFRelaxationOperator
/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"

#include "CellVariable.h"
#include "HierarchyDataOpsManager.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFRelaxationOperator is an abstract class that provides an interface for specifying a relaxation
 * function for the extra stress for Oldroyd-B type viscoelastic fluid models.
 *
 * The patch data index specifying either the conformation tensor or the square root or logarithm of the conformation
 * tensor is given to the class. The function setDataOnPatchHierarchy is called by the advection diffusion integrator
 * and expects the relaxation tensor.
 */
class CFRelaxationOperator : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor does nothing interesting.
     */
    CFRelaxationOperator(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CFRelaxationOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CFRelaxationOperator(const CFRelaxationOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     */
    CFRelaxationOperator& operator=(const CFRelaxationOperator& that) = delete;

    /*!
     * \brief Empty destructor.
     */
    virtual ~CFRelaxationOperator() = default;

    /*!
     * \brief Sets the index for the transported quantity.
     */
    void setPatchDataIndex(int);

    /*!
     * \brief Indicates whether the concrete CFRelaxationOperator object is
     * time-dependent. Returns true.
     */
    bool isTimeDependent() const override;

protected:
    /*!
     * \brief This function converts the data stored in the patch data index to the conformation tensor. This has a
     * default implementation that converts from the logarithm or square root to the full conformation tensor.
     */
    virtual IBTK::MatrixNd convertToConformation(const IBTK::MatrixNd& mat);

    int d_W_cc_idx = IBTK::invalid_index;

private:
    TensorEvolutionType d_evolve_type = STANDARD;
};

} // Namespace IBAMR
#endif
