#ifndef BZ_RANKS_H
#define BZ_RANKS_H

BZ_NAMESPACE(blitz)

/*
 * Rank numbers start with zero, which may be confusing to users coming
 * from Fortran.  To make code more readable, the following constants
 * may help.  Example: instead of
 *
 * int firstRankExtent = A.extent(0);
 *
 * One can write:
 *
 * int firstRankExtent = A.extent(firstRank);
 */

const int firstRank    = 0;
const int secondRank   = 1;
const int thirdRank    = 2;
const int fourthRank   = 3;
const int fifthRank    = 4;
const int sixthRank    = 5;
const int seventhRank  = 6;
const int eighthRank   = 7;
const int ninthRank    = 8;
const int tenthRank    = 9;
const int eleventhRank = 10;

const int firstDim    = 0;
const int secondDim   = 1;
const int thirdDim    = 2;
const int fourthDim   = 3;
const int fifthDim    = 4;
const int sixthDim    = 5;
const int seventhDim  = 6;
const int eighthDim   = 7;
const int ninthDim    = 8;
const int tenthDim    = 9;
const int eleventhDim = 10;

BZ_NAMESPACE_END

#endif
