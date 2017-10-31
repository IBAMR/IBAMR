// Filename: RelaxationLSMethod.h
// Created on 10 Oct 2017 by Amneet Bhalla and Nishant Nangia
//
// Copyright (c) 2002-2014, Amneet Bhalla and Nishant Nangia
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
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HinitER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_RelaxationLSMethod
#define included_IBAMR_RelaxationLSMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "ibamr/LSInitStrategy.h"
#include "ibamr/ibamr_enums.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
} // namespace pdat
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class RelaxationLSMethod provides a relaxation algorithm implementation
 * of the level set method. Specifically, this class iterates (to steady-state) the PDE
 * \f$\frac{\partial Q}{\partial \tau}+sgn(Q^0)(|\nabla Q | - 1) = 0\f$, which produces a solution
 * to the Eikonal equation \f$ |\nabla Q | = 1 \f$. The solution of the Eikonal equation
 * produces the signed distance away from an interface.
 *
 *
 * Reference
 * Min, C., <A HREF="http://www.sciencedirect.com/science/article/pii/S0021999109007189">
 * On reinitializing level set functions</A>
 */
class RelaxationLSMethod : public IBAMR::LSInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    RelaxationLSMethod(const std::string& object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = NULL,
                       bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~RelaxationLSMethod();

    /*!
     * \brief Initialize level set data using the fast-sweeping method.
     */
    void initializeLSData(int D_idx,
                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hierarchy_math_ops,
                          double time,
                          bool initial_time);

    /////////////////////////////// PROTECTED ////////////////////////////////////

protected:

    /////////////////////////////// PRIVATE //////////////////////////////////////

private:
    /*!
     * \brief Do one relaxation step over the hierarchy.
     */
    void relax(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
               int dist_idx,
               int dist_init_idx,
               const int iter) const;

    /*!
     * \brief Do one relaxation step over a patch.
     */
    void relax(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_data,
               const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_init_idx,
               const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
               const int iter) const;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RelaxationLSMethod(const RelaxationLSMethod& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RelaxationLSMethod& operator=(const RelaxationLSMethod& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_RelaxationLSMethod
