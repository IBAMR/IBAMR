#ifndef included_IBMarkerRefineOperator
#define included_IBMarkerRefineOperator

// Filename: IBMarkerRefineOperator.h
// Last modified: <04.Oct.2007 23:26:42 griffith@box221.cims.nyu.edu>
// Created on 04 Oct 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <RefineOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMarkerRefineOperator is a concrete
 * SAMRAI::xfer::RefineOperator for refining IB marker data from coarser levels
 * to finer levels in the patch hierarchy.
 */
class IBMarkerRefineOperator
    : public SAMRAI::xfer::RefineOperator<NDIM>
{
public:
    /*!
     * \brief Defaultonstructor.
     */
    IBMarkerRefineOperator();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBMarkerRefineOperator();

    /*!
     * \name Implementation of SAMRAI::xfer::RefineOperator interface.
     */
    //\{

    /*!
     * Return true if the refining operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    virtual bool
    findRefineOperator(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
        const std::string& op_name) const;

    /*!
     * Return name string identifier of the refining operation.
     */
    virtual const std::string&
    getOperatorName() const;

    /*!
     * Return the priority of this operator relative to other refining
     * operators.  The SAMRAI transfer routines guarantee that refining using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    virtual int
    getOperatorPriority() const;

    /*!
     * Return the stencil width associated with the refining operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each refining operator.
     */
    virtual SAMRAI::hier::IntVector<NDIM>
    getStencilWidth() const;

    /*!
     * Refine the source component on the fine patch to the destination
     * component on the coarse patch. The refining operation is performed on the
     * intersection of the destination patch and the coarse box.  The fine patch
     * is guaranteed to contain sufficient data for the stencil width of the
     * refining operator.
     */
    virtual void
    refine(
        SAMRAI::hier::Patch<NDIM>& fine,
        const SAMRAI::hier::Patch<NDIM>& coarse,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::Box<NDIM>& fine_box,
        const SAMRAI::hier::IntVector<NDIM>& ratio) const;

    //\}

protected:

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMarkerRefineOperator(
        const IBMarkerRefineOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMarkerRefineOperator&
    operator=(
        const IBMarkerRefineOperator& that);

    /*!
     * The operator name.
     */
    static const std::string s_op_name;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBMarkerRefineOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMarkerRefineOperator
