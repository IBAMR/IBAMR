//
// LagDataIO.h
//
// Created on 09 Jun 2005
//         by Boyce Griffith (boyce@bigboy.verizon.net).
//
// Last modified: <13.Jun.2005 21:53:51 boyce@mstu1.cims.nyu.edu>
//

#ifndef included_LagDataIO
#define included_LagDataIO

// STL INCLUDES
//
#include <string>

// BLITZ++ INCLUDES
//
#include <blitz/array.h>

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "tbox/Utilities.h"

using namespace SAMRAI;
using namespace std;

// CLASS DEFINITION
//

/*!
 * @brief Description of class.
 */
#if 0
class LagDataIO
{
public:
    /*!
     * @brief Default constructor.
     */
    LagDataIO();
    
    /*!
     * @brief Copy constructor.
     *
     * @param from The value to copy to this object.
     */
    LagDataIO(
        const LagDataIO& from);
    
    /*!
     * @brief Destructor.
     */
    ~LagDataIO();
    
    /*!
     * @brief Assignment operator.
     *
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LagDataIO& operator=(
        const LagDataIO& that);
    
    static void readFibers(
        const string& file_name,
        int& ngroups,
        blitz::Array<int,1>& nfiber_per_group,
        blitz::Array<int,1>& npoint_per_fiber,
        blitz::Array<blitz::Array<int,2>,1>& point_idxs,
        blitz::Array<blitz::Array<double,3>,1>& point_locs,
        const bool compute_idxs=false);
    
    static void readMarkers(
        const string& file_name,
        int& nclouds,
        blitz::Array<int,1>& nmarks,
        blitz::Array<blitz::Array<int,1>,1>& marker_idxs,
        blitz::Array<blitz::Array<double,2>,1>& marker_locs,
        const bool compute_idxs=false);

protected:
    
private:
};
#endif

// INLINED FUNCTION DEFINITIONS
//
//#ifndef DEBUG_NO_INLINE
//#include "LagDataIO.I"
//#endif

#endif //#ifndef included_LagDataIO

//////////////////////////////////////////////////////////////////////////////
