// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartCellDoubleLinearCoarsen
#define included_IBTK_CartCellDoubleLinearCoarsen

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <CoarsenOperator.h>
#include <IntVector.h>

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
} // namespace SAMRAI

namespace IBTK
{
class CartCellDoubleLinearCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    CartCellDoubleLinearCoarsen() = default;
    ~CartCellDoubleLinearCoarsen() = default;

    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                             const std::string& op_name) const override;
    const std::string& getOperatorName() const override;
    int getOperatorPriority() const override;
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const override;
    void coarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                 const SAMRAI::hier::Patch<NDIM>& fine,
                 int dst_component,
                 int src_component,
                 const SAMRAI::hier::Box<NDIM>& coarse_box,
                 const SAMRAI::hier::IntVector<NDIM>& ratio) const override;

private:
    CartCellDoubleLinearCoarsen(const CartCellDoubleLinearCoarsen& from) = delete;
    CartCellDoubleLinearCoarsen& operator=(const CartCellDoubleLinearCoarsen& that) = delete;

    static const std::string s_op_name;
};
} // namespace IBTK

#endif
