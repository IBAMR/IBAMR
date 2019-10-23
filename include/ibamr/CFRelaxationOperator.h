// Filename: CFRelaxationOperator.h
// Created on 23 Oct 2019 by Aaron Barrett
//
// Copyright (c) 2002-2019, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_CFRelaxationOperator
#define included_CFRelaxationOperator
/////////////////////////////// INCLUDES /////////////////////////////////////
#include "IBAMR_config.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "HierarchyDataOpsManager.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <unsupported/Eigen/MatrixFunctions>

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
