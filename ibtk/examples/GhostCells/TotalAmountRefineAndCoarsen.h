// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_TotalAmountRefineAndCoarsen
#define included_TotalAmountRefineAndCoarsen

#include <ibtk/samrai_compatibility_names.h>

#include <IBTK_config.h>
#include <SAMRAIBox.h>
#include <SAMRAICoarsenOperator.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRefineOperator.h>
#include <SAMRAIVariable.h>

#include <string>

namespace IBTK
{
/*!
 * \brief Class TotalAmountRefine is a concrete SAMRAI::xfer::RefineOperator for refining a cell centered quantity that
 * represents the total amount of stuff in a cell. We assume a constant profile across the entire coarse cell.
 *
 * Refine operators must implement five methods to be valid extensions of RefineOperator. We must also register the new
 * refine operator with the GridGeometry before we can use it.
 */
class TotalAmountRefine : public SAMRAIRefineOperator
{
public:
    /*!
     * \brief Default constructor.
     */
    TotalAmountRefine() = default;

    /*!
     * \brief Default destructor.
     */
    ~TotalAmountRefine() = default;

    /*!
     * \brief Returns true if the refinement operation matches the given variable and operation name.
     */
    bool findRefineOperator(const SAMRAIPointer<SAMRAIVariable>& var, const std::string& op_name) const override;
    /*!
     * \brief Returns the operator name.
     */
    const std::string& getOperatorName() const override;

    /*!
     * \brief Returns the priority of the operator. SAMRAI guarantees that operators with lower priority will be
     * performed first.
     */
    int getOperatorPriority() const override;

    /*!
     * \brief Return the stencil width of the operator.
     */
    SAMRAIIntVector getStencilWidth() const override;

    /*!
     * \brief Refine the source component into the destination component in the box fine_box. Note the coarse patch is
     * gaurenteed to have sufficient data for the stencil width of the operator.
     */
    void refine(SAMRAIPatch& fine,
                const SAMRAIPatch& coarse,
                const int dst_component,
                const int src_component,
                const SAMRAIBox& fine_box,
                const SAMRAIIntVector& ratio) const override;

private:
    static std::string s_object_name;
};

/*!
 * \brief Class TotalAmountCoarsen is a concrete SAMRAI::xfer::CoarsenOperator for coarsening a cell centered quantity
 * that represents the total amount of stuff in a cell.
 *
 * Coarsen operators must implement five methods to be valid extensions of CoarsenOperator. We must also register the
 * new coarsen operator with the GridGeometry before we can use it.
 */
class TotalAmountCoarsen : public SAMRAICoarsenOperator
{
public:
    /*!
     * \brief Default constructor.
     */
    TotalAmountCoarsen() = default;

    /*!
     * \brief Default destructor.
     */
    ~TotalAmountCoarsen() = default;

    /*!
     * \brief Returns true if the coarsening operation matches the given variable and operation name.
     */
    bool findCoarsenOperator(const SAMRAIPointer<SAMRAIVariable>& var, const std::string& op_name) const override;

    /*!
     * \brief Returns the operator name.
     */
    const std::string& getOperatorName() const override;

    /*!
     * \brief Returns the priority of the operator. SAMRAI guarantees that operators with lower priority will be
     * performed first.
     */
    int getOperatorPriority() const override;

    /*!
     * \brief Return the stencil width of the operator.
     */
    SAMRAIIntVector getStencilWidth() const override;

    /*!
     * \brief Coarsen the source component on the fine patch to the destination component on the coarse patch. The fine
     * patch is guaranteed to contain sufficient ghost cell width for the stencil width of the operator.
     */
    void coarsen(SAMRAIPatch& coarse,
                 const SAMRAIPatch& fine,
                 const int dst_component,
                 const int src_component,
                 const SAMRAIBox& coarse_box,
                 const SAMRAIIntVector& ratio) const override;

private:
    static std::string s_object_name;
};

} // namespace IBTK
#endif
