#ifndef included_LNodeIndexVariable2
#define included_LNodeIndexVariable2

// Filename: LNodeIndexVariable2.h
// Last modified: <04.Jun.2007 12:37:29 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <Variable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeIndexVariable2 provides a SAMRAI::hier::Variable class
 * corresponding to patch data of type LNodeIndexData2.
 */
class LNodeIndexVariable2
    : public SAMRAI::hier::Variable<NDIM>
{
public:
    /*!
     * Create an LNodeIndexVariable2 object with the specified name.
     */
    LNodeIndexVariable2(
        const std::string& name);

    /*!
     * Virtual destructor for LNodeIndexVariable2 objects.
     */
    virtual
    ~LNodeIndexVariable2();

    /*!
     * Return false since the LNodeIndex data index space matches the
     * cell-centered index space for AMR patches.  Thus, LNodeIndex data does
     * not live on patch borders.
     */
    bool
    dataLivesOnPatchBorder() const;

    /*!
     * Return true so that the LNodeIndex data quantities will always be treated
     * as though fine values represent them on coarse-fine interfaces.  Note
     * that this is really artificial since the LNodeIndex data index space
     * matches the cell-centered index space for AMR patches.  Thus, LNodeIndex
     * data does not live on patch borders and so there is no ambiguity
     * reagrding coarse-fine interface values.
     */
    bool
    fineBoundaryRepresentsVariable() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexVariable2();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexVariable2(
        const LNodeIndexVariable2& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexVariable2&
    operator=(
        const LNodeIndexVariable2& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LNodeIndexVariable2.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexVariable2
