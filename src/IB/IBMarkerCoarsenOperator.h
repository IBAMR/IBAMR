#ifndef included_IBMarkerCoarsenOperator
#define included_IBMarkerCoarsenOperator

// Filename: IBMarkerCoarsenOperator.h
// Last modified: <21.Sep.2007 21:43:41 griffith@box221.cims.nyu.edu>
// Created on 13 Sep 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CoarsenOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMarkerCoarsenOperator is a concrete
 * SAMRAI::xfer::CoarsenPatchStrategy for coarsening IB marker data from finer
 * levels to coarser levels in the patch hierarchy.
 */
class IBMarkerCoarsenOperator
    : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    IBMarkerCoarsenOperator();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBMarkerCoarsenOperator();

    /*!
     * \name Implementation of SAMRAI::xfer::CoarsenOperator interface.
     */
    //\{

    /*!
     * Return true if the coarsening operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    virtual bool
    findCoarsenOperator(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
        const std::string& op_name) const;

    /*!
     * Return name string identifier of the coarsening operation.
     */
    virtual const std::string&
    getOperatorName() const;

    /*!
     * Return the priority of this operator relative to other coarsening
     * operators.  The SAMRAI transfer routines guarantee that coarsening using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    virtual int
    getOperatorPriority() const;

    /*!
     * Return the stencil width associated with the coarsening operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each coarsening operator.
     */
    virtual SAMRAI::hier::IntVector<NDIM>
    getStencilWidth() const;

    /*!
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch. The coarsening operation is performed on
     * the intersection of the destination patch and the coarse box.  The fine
     * patch is guaranteed to contain sufficient data for the stencil width of
     * the coarsening operator.
     */
    virtual void
    coarsen(
        SAMRAI::hier::Patch<NDIM>& coarse,
        const SAMRAI::hier::Patch<NDIM>& fine,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::Box<NDIM>& coarse_box,
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
    IBMarkerCoarsenOperator(
        const IBMarkerCoarsenOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMarkerCoarsenOperator&
    operator=(
        const IBMarkerCoarsenOperator& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <stools/IBMarkerCoarsenOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMarkerCoarsenOperator
