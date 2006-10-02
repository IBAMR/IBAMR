//
// LNodeIndexData.h
//
// Created on 01 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <08.Mar.2005 12:41:14 boyce@trasnaform.cims.nyu.edu>
//

#ifndef included_LNodeIndexData
#define included_LNodeIndexData

// STL INCLUDES
//
#include <vector>

// SAMRAI-tools INCLUDES
//
#include "LNodeIndexSet.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "Box.h"
#include "IndexData.h"
#include "IntVector.h"

using namespace SAMRAI;
using namespace std;

// FORWARD DECLARATIONS
//
class LDataManager;

// CLASS DEFINITION
//

/*!
 * @brief Class LNodeIndexData is a specialization of the templated
 * pdat::IndexData<NDIM> class.  It extends the pdat::IndexData<NDIM> class by providing the
 * IDs of the Lagrangian nodes in the patch interior and in the ghost
 * cell region of the patch.
 */
class LNodeIndexData
    : public pdat::IndexData<NDIM,LNodeIndexSet>
{
public:
    //@{ @name Friend declarations.
    friend class pdat::IndexIterator<NDIM,LNodeIndexSet>;
    friend class LDataManager;
    //@}
    
    /*!
     * The constructor for an pdat::IndexData<NDIM> object.  The box describes the
     * interior of the index space and the ghosts vector describes the
     * ghost nodes in each coordinate direction.
     */
    LNodeIndexData(
        const hier::Box<NDIM>& box,
        const hier::IntVector<NDIM>& ghosts);
    
    /*!
     * The virtual destructor for an LNodeIndexData object.
     */
    virtual ~LNodeIndexData();
    
    /*!
     * The hier::PatchData<NDIM> functionality of LNodeIndexData is provided by
     * the pdat::IndexData<NDIM,LNodeIndexSet> class.
     */
    using pdat::IndexData<NDIM,LNodeIndexSet>::copy;
    using pdat::IndexData<NDIM,LNodeIndexSet>::copy2;
    using pdat::IndexData<NDIM,LNodeIndexSet>::canEstimateStreamSizeFromBox;
    using pdat::IndexData<NDIM,LNodeIndexSet>::getDataStreamSize;
    using pdat::IndexData<NDIM,LNodeIndexSet>::packStream;
    using pdat::IndexData<NDIM,LNodeIndexSet>::unpackStream;
    using pdat::IndexData<NDIM,LNodeIndexSet>::appendItem;
    using pdat::IndexData<NDIM,LNodeIndexSet>::addItem;
    using pdat::IndexData<NDIM,LNodeIndexSet>::removeItem;
    using pdat::IndexData<NDIM,LNodeIndexSet>::removeInsideBox;
    using pdat::IndexData<NDIM,LNodeIndexSet>::removeOutsideBox;
    using pdat::IndexData<NDIM,LNodeIndexSet>::removeGhostItems;
    using pdat::IndexData<NDIM,LNodeIndexSet>::removeAllItems;
    using pdat::IndexData<NDIM,LNodeIndexSet>::isElement;
    using pdat::IndexData<NDIM,LNodeIndexSet>::getItem;
    using pdat::IndexData<NDIM,LNodeIndexSet>::getSpecializedFromDatabase;
    using pdat::IndexData<NDIM,LNodeIndexSet>::putSpecializedToDatabase;

    /*!
     * @return A constant refrence to the set of local PETSc indices
     * to Lagrangian nodes which lie in the patch interior.
     *
     * NOTE: This should be an ordered, contiguous set of indices.
     */
    const vector<int>& getInteriorLocalIndices() const;
    
    /*!
     * @return A constant refrence to the set of local PETSc indices
     * to Lagrangian nodes which lie in ghost cell region of this
     * hier::PatchData<NDIM> object.
     */
    const vector<int>& getGhostLocalIndices() const;
    
private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    LNodeIndexData();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    LNodeIndexData(
        const LNodeIndexData& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LNodeIndexData& operator=(
        const LNodeIndexData& that);

    /*!
     * Each LNodeIndexData object maintains collections of local
     * (PETSc) indices corresponding to the local and ghost nodes
     * which lie in the patch data ghost box.
     */
    vector<int> d_interior_local_indices;
    vector<int> d_ghost_local_indices;
};

// INLINED FUNCTION DEFINITIONS
//
#ifndef DEBUG_NO_INLINE
#include "LNodeIndexData.I"
#endif

#endif //#ifndef included_LNodeIndexData

//////////////////////////////////////////////////////////////////////////////
