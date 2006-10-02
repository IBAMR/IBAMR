//
// LNodeIndexData.I
//
// Created on 17 Apr 2004
//         by Boyce Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <11.Jul.2004 21:38:52 boyce@trasnaform.speakeasy.net>
//

//////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#define inline
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
const vector<int>& LNodeIndexData::getInteriorLocalIndices() const
{
    return(d_interior_local_indices);
}// getInteriorLocalIndices

inline
const vector<int>& LNodeIndexData::getGhostLocalIndices() const
{
    return(d_ghost_local_indices);
}// getGhostLocalIndices

//////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

//////////////////////////////////////////////////////////////////////////////
