// Filename: LNodeIndexData2.I
// Last modified: <04.Jun.2007 15:20:05 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline bool
LNodeIndexData2::canEstimateStreamSizeFromBox() const
{
    return false;
}// canEstimateStreamSizeFromBox

inline const std::vector<int>&
LNodeIndexData2::getInteriorLocalIndices() const
{
    return d_interior_local_indices;
}// getInteriorLocalIndices

inline const std::vector<int>&
LNodeIndexData2::getGhostLocalIndices() const
{
    return d_ghost_local_indices;
}// getGhostLocalIndices

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
