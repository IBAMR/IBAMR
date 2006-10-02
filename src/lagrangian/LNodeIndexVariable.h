//
// LNodeIndexVariable.h
//
// Created on 01 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <07.Mar.2005 17:56:31 boyce@trasnaform.cims.nyu.edu>
//

#ifndef included_LNodeIndexVariable
#define included_LNodeIndexVariable

// SAMRAI-tools INCLUDES
//
#include "LNodeIndexSet.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "Variable.h"

using namespace SAMRAI;
using namespace std;

// CLASS DEFINITION
//

/*!
 * @brief Class LNodeIndexVariable provides a hier::Variable<NDIM> class
 * corresponding to patch data of type LNodeIndexData.
 */
class LNodeIndexVariable
    : public hier::Variable<NDIM>
{
public:
    /*!
     * Create an LNodeIndexVariable object with the specified name.
     */
    LNodeIndexVariable(
        const string& name);
    
    /*!
     * Virtual destructor for LNodeIndexVariable objects.
     */
    virtual ~LNodeIndexVariable();

    /*!
     * Return false since the LNodeIndex data index space matches the
     * cell-centered index space for AMR patches.  Thus, LNodeIndex
     * data does not live on patch borders.
     */
    bool dataLivesOnPatchBorder() const;

    /*!
     * Return true so that the LNodeIndex data quantities will always
     * be treated as though fine values represent them on coarse-fine
     * interfaces.  Note that this is really artificial since the
     * LNodeIndex data index space matches the cell-centered index
     * space for AMR patches.  Thus, LNodeIndex data does not live on
     * patch borders and so there is no ambiguity reagrding
     * coarse-fine interface values.
     */
    bool fineBoundaryRepresentsVariable() const;
    
private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    LNodeIndexVariable();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    LNodeIndexVariable(
        const LNodeIndexVariable& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LNodeIndexVariable& operator=(
        const LNodeIndexVariable& that);
};

#endif //#ifndef included_LNodeIndexVariable

//////////////////////////////////////////////////////////////////////////////
