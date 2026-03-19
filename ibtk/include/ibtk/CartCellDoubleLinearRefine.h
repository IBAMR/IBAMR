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

#ifndef included_IBTK_CartCellDoubleLinearRefine
#define included_IBTK_CartCellDoubleLinearRefine

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <IntVector.h>
#include <RefineOperator.h>

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
class CartCellDoubleLinearRefine : public SAMRAI::xfer::RefineOperator<NDIM>
{
public:
    CartCellDoubleLinearRefine() = default;
    ~CartCellDoubleLinearRefine() = default;

    bool findRefineOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                            const std::string& op_name) const override;
    const std::string& getOperatorName() const override;
    int getOperatorPriority() const override;
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const override;
    void refine(SAMRAI::hier::Patch<NDIM>& fine,
                const SAMRAI::hier::Patch<NDIM>& coarse,
                int dst_component,
                int src_component,
                const SAMRAI::hier::Box<NDIM>& fine_box,
                const SAMRAI::hier::IntVector<NDIM>& ratio) const override;

private:
    CartCellDoubleLinearRefine(const CartCellDoubleLinearRefine& from) = delete;
    CartCellDoubleLinearRefine& operator=(const CartCellDoubleLinearRefine& that) = delete;

    static const std::string s_op_name;
};
} // namespace IBTK

#endif
