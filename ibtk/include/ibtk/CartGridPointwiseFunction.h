// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartGridPointwiseFunction
#define included_IBTK_CartGridPointwiseFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Function signatures for CartGridPointwiseFunction
 */
namespace PointwiseFunctions
{
using ScalarFcn = std::function<double(double, const VectorNd&, double)>;
using VectorFcn = std::function<VectorNd(const VectorNd&, const VectorNd&, double)>;
using OtherFcn = std::function<VectorXd(const VectorXd&, const VectorNd&, double)>;
using StaggeredFcn = std::function<double(double, const VectorNd&, double, int)>;
} // namespace PointwiseFunctions

/*!
 * CartGridPointwiseFunction provides a lightweight class that can evaluate a function pointwise on the Cartesian grid.
 * The constructor takes a function argument that will be applied to each index in the hierarchy.
 *
 * The valid function signature depends on the data centering and depth of the variable:
 *
 * For cell or node centered scalar quantities, the valid function signature must match PointwiseFunctions::ScalarFcn.
 * For cell or node centered vector quantities, the valid function signature must match PointwiseFunctions::VectorFcn.
 * For other cell or node centered quantities, the valid function signature must match PointwiseFunctions::OtherFcn.
 * For side, face, or edge centered quantities, the valid function signature must match
 * PointwiseFunctions::StaggeredFcn.
 *
 * The arguments given to the function include the current value stored on that index, the physical location of the
 * point, and the time. For side, face, or edge centered quantities, the axis is also provided.
 *
 */
template <typename F>
class CartGridPointwiseFunction : public CartGridFunction
{
public:
    /*!
     * \brief Constructor. The function f must match on the signatures provided by PointwiseFunctions.
     */
    CartGridPointwiseFunction(std::string object_name, F f);

    /*!
     * \brief Does this function depend on time?
     */
    bool isTimeDependent() const override
    {
        return true;
    }

    /*!
     * \brief Evaluate the function on the provided patch.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

private:
    F d_f;
};

} // namespace IBTK

#endif // #ifndef included_IBTK_CartGridPointwiseFunction
