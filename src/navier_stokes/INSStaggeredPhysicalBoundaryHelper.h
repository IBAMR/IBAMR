// Filename: INSStaggeredPhysicalBoundaryHelper.h
// Created on 22 Jul 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef included_INSStaggeredPhysicalBoundaryHelper
#define included_INSStaggeredPhysicalBoundaryHelper

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <Variable.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredPhysicalBoundaryHelper provides various helper
 * functions required to specify physical boundary conditions for a staggered
 * grid discretization of the incompressible Navier-Stokes equations.
 */
class INSStaggeredPhysicalBoundaryHelper
    : SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    INSStaggeredPhysicalBoundaryHelper();

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredPhysicalBoundaryHelper();

    /*!
     * \brief Set values located on the physical boundary to zero on the
     * specified range of levels in the patch hierarchy using the cached
     * boundary data.
     *
     * \note By default, boundary conditions are cached over the complete range
     * of levels of the patch hierarchy.
     */
    void
    zeroValuesAtDirichletBoundaries(
        const int patch_data_idx,
        const int coarsest_level_number=-1,
        const int finest_ln=-1) const;

    /*!
     * \brief Reset boundary values located on the physical boundary to zero on
     * the specified range of levels in the patch hierarchy using the cached
     * boundary data.
     *
     * \note By default, boundary conditions are cached over the complete range
     * of levels of the patch hierarchy.
     */
    void
    resetValuesAtDirichletBoundaries(
        const int patch_data_idx,
        const int coarsest_level_number=-1,
        const int finest_ln=-1) const;

    /*!
     * \brief Cache boundary coefficient data.
     */
    void
    cacheBcCoefData(
        const int u_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& u_var,
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
        const double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& gcw_to_fill,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >& hierarchy);

    /*!
     * \brief Clear cached boundary coefficient data.
     */
    void
    clearBcCoefData();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredPhysicalBoundaryHelper(
        const INSStaggeredPhysicalBoundaryHelper& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredPhysicalBoundaryHelper&
    operator=(
        const INSStaggeredPhysicalBoundaryHelper& that);

    /*!
     * Cached hierarchy-related information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    std::vector<std::map<int,std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,bool> > > > > d_dirichlet_bdry_locs;
    std::vector<std::map<int,std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > > > > d_dirichlet_bdry_vals;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredPhysicalBoundaryHelper.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredPhysicalBoundaryHelper
