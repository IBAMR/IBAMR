//
// LNodeIndexSet.I
//
// Created on 29 Feb 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <29.Jun.2005 21:42:43 boyce@mstu1.cims.nyu.edu>
//

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#define inline
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNodeIndexSet::LNodeIndexSet()
    : d_set(),
      d_offset(0)
{
    // intentionally blank
    return;
}// LNodeIndexSet

inline
LNodeIndexSet::LNodeIndexSet(
    const LNodeIndexSet& from)
{
    d_set    = from.d_set   ;
    d_offset = from.d_offset;
    return;
}// LNodeIndexSet

inline
LNodeIndexSet::~LNodeIndexSet()
{
    d_set.clear();
    return;
}// ~LNodeIndexSet

inline
LNodeIndexSet& LNodeIndexSet::operator=(
    const LNodeIndexSet& that)
{
    if (this == &that) return(*this);  // check for self-assignment

    d_set    = that.d_set   ;
    d_offset = that.d_offset;
    return(*this);
}// operator=

inline
LNodeIndexSet::const_iterator LNodeIndexSet::begin() const
{
    return(d_set.begin());
}// begin

inline
LNodeIndexSet::iterator LNodeIndexSet::begin()
{
    return(d_set.begin());
}// begin

inline
LNodeIndexSet::const_iterator LNodeIndexSet::end() const
{
    return(d_set.end());
}// end

inline
LNodeIndexSet::iterator LNodeIndexSet::end()
{
    return(d_set.end());
}// end

inline
LNodeIndexSet::size_type LNodeIndexSet::size() const
{
    return(d_set.size());
}// size

inline
bool LNodeIndexSet::empty() const
{
    return(d_set.empty());
}// empty

inline
void LNodeIndexSet::push_back(
    const LNodeIndexSet::value_type& value)
{
    d_set.push_back(value);
    return;
}// push_back

inline
const hier::IntVector<NDIM>& LNodeIndexSet::getPeriodicOffset() const
{
    return(d_offset);
}// getPeriodicOffset

inline
void LNodeIndexSet::setPeriodicOffset(
    const hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
}// setPeriodicOffset

inline
void LNodeIndexSet::copySourceItem(
    const hier::Index<NDIM>& src_index,
    const hier::IntVector<NDIM>& src_offset,
    const LNodeIndexSet& src_item)
{
    d_set    = src_item.d_set;
    d_offset = src_offset;
    return;
}// copySourceItem

/////////////////////////////// PRIVATE //////////////////////////////////////

inline
void LNodeIndexSet::trimToFit()
{
    if (d_set.size() < d_set.capacity())
    {
        IndexSet(d_set).swap(d_set);
    }
    return;
}// trimToFit

//////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

//////////////////////////////////////////////////////////////////////////////
