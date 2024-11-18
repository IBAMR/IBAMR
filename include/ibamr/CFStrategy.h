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

#ifndef included_CFStrategy
#define included_CFStrategy
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

IBTK_DISABLE_EXTRA_WARNINGS
#include <unsupported/Eigen/MatrixFunctions>
IBTK_ENABLE_EXTRA_WARNINGS

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFStrategy is an abstract class that provides an interface for specifying the details of a complex fluid
 * model.
 *
 * Derived classes must implement the <code>computeRelaxation</code> and <code>computeStress</code> functions.
 *
 * The patch data index specifying either the conformation tensor or the square root or logarithm of the conformation
 * tensor is given to the class. The function setDataOnPatchHierarchy is called by the advection diffusion integrator
 * and expects the relaxation tensor.
 */
class CFStrategy : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief This constructor does nothing interesting.
     */
    CFStrategy(std::string object_name);

    /*!
     * \brief Compute the action of the relaxation operator, given the conformation tensor.
     *
     * The parameter <code>evolve_type</code> determines the form of the tensor provided in <code>C_idx</code>. A helper
     * function to convert from the data stored in <code>C_idx</code> to the conformation tensor is provided by
     * <code>convert_to_conformation_tensor</code>.
     *
     * The symmetric tensor in patch data index <code>C_idx</code> is stored using Voigt notation.
     */
    virtual void computeRelaxation(int R_idx,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > R_var,
                                   int C_idx,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > C_var,
                                   TensorEvolutionType evolve_type,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                   double data_time) = 0;

    /*!
     * \brief Convert the conformation tensor to the stress tensor.
     *
     * This conversion should be done in-place. The conformation tensor is provided in the patch index
     * <code>sig_idx</code> with ghost cells filled. Implementations of this function MUST provide correct ghost cell
     * information in <code>sig_idx</code>.
     */
    virtual void computeStress(int sig_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > sig_var,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               double data_time) = 0;

private:
    std::string d_object_name;
};

/*!
 * Convert the provided matrix in the form given by evolve_type to the conformation tensor.
 */
inline void
convert_to_conformation_tensor(IBTK::MatrixNd& mat, TensorEvolutionType evolve_type)
{
    switch (evolve_type)
    {
    case SQUARE_ROOT:
        mat *= mat;
        break;
    case LOGARITHM:
        mat = mat.exp();
        break;
    case STANDARD:
    default:
        // Do nothing
        break;
    }
}

/*!
 * Return the conformation tensor, given the evolved data and the evolved type.
 */
inline IBTK::MatrixNd
convert_to_conformation_tensor(const SAMRAI::pdat::CellData<NDIM, double>& C_data,
                               const SAMRAI::pdat::CellIndex<NDIM>& idx,
                               TensorEvolutionType evolve_type)
{
    IBTK::MatrixNd mat;
#if (NDIM == 2)
    mat(0, 0) = C_data(idx, 0);
    mat(1, 1) = C_data(idx, 1);
    mat(0, 1) = mat(1, 0) = C_data(idx, 2);
#endif
#if (NDIM == 3)
    mat(0, 0) = C_data(idx, 0);
    mat(1, 1) = C_data(idx, 1);
    mat(2, 2) = C_data(idx, 2);
    mat(1, 2) = mat(2, 1) = C_data(idx, 3);
    mat(0, 2) = mat(2, 0) = C_data(idx, 4);
    mat(0, 1) = mat(1, 0) = C_data(idx, 5);
#endif
    convert_to_conformation_tensor(mat, evolve_type);
    return mat;
}

} // Namespace IBAMR
#endif
