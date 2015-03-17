// Filename: AdvDiffPhysicalBoundaryUtilities.h
// Created on 24 Aug 2012 by Boyce Griffith
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

#ifndef included_AdvDiffPhysicalBoundaryUtilities
#define included_AdvDiffPhysicalBoundaryUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

namespace SAMRAI
{
namespace hier
{
class Patch;
} // namespace hier
namespace pdat
{
template <class TYPE>
class CellData;
template <class TYPE>
class FaceData;
} // namespace pdat
namespace solv
{
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffPhysicalBoundaryUtilities is a utility class that
 * provides functions useful for dealing with physical boundary conditions for
 * advection-diffusion solvers.
 */
class AdvDiffPhysicalBoundaryUtilities
{
public:
    /*!
     * \brief Set physical boundary conditions at physical boundaries.
     */
    static void setPhysicalBoundaryConditions(const boost::shared_ptr<SAMRAI::pdat::CellData<double> >& Q_data,
                                              const boost::shared_ptr<SAMRAI::pdat::FaceData<double> >& u_ADV_data,
                                              const boost::shared_ptr<SAMRAI::hier::Patch>& patch,
                                              const std::vector<boost::shared_ptr<SAMRAI::solv::RobinBcCoefStrategy>>& bc_coefs,
                                              const double fill_time,
                                              const bool inflow_boundaries_only,
                                              const bool homogeneous_bc);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPhysicalBoundaryUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPhysicalBoundaryUtilities(const AdvDiffPhysicalBoundaryUtilities& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPhysicalBoundaryUtilities& operator=(const AdvDiffPhysicalBoundaryUtilities& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffPhysicalBoundaryUtilities
