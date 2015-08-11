// Filename: IBMethodPostProcessStrategy.h
// Created on 24 Sep 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBMethodPostProcessStrategy
#define included_IBMethodPostProcessStrategy

#include <vector>

#include "boost/shared_ptr.hpp"

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{

class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

namespace IBAMR
{
class IBMethod;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMethodPostProcessStrategy provides a generic interface for
 * specifying post-processing code for use in an IB computation.
 */
class IBMethodPostProcessStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    IBMethodPostProcessStrategy();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBMethodPostProcessStrategy();

    /*!
     * \brief Post-process data on the patch hierarchy.
     */
    virtual void postprocessData(int u_idx,
                                 int p_idx,
                                 int f_idx,
                                 const std::vector<boost::shared_ptr<IBTK::LData>>& F_data,
                                 const std::vector<boost::shared_ptr<IBTK::LData>>& X_data,
                                 const std::vector<boost::shared_ptr<IBTK::LData>>& U_data,
                                 const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                 int coarsest_level_number,
                                 int finest_level_number,
                                 double data_time,
                                 IBMethod* ib_method) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMethodPostProcessStrategy(const IBMethodPostProcessStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMethodPostProcessStrategy& operator=(const IBMethodPostProcessStrategy& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMethodPostProcessStrategy
