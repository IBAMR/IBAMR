#ifndef included_INSStaggeredPhysicalBoundaryHelper
#define included_INSStaggeredPhysicalBoundaryHelper

// Filename: INSStaggeredPhysicalBoundaryHelper.h
// Last modified: <24.Nov.2008 16:15:05 griffith@box230.cims.nyu.edu>
// Created on 22 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

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
     * specified range of levels in the patch hierarchy.
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
     * \brief Cache boundary coefficient data.
     */
    void
    cacheBcCoefData(
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
    std::vector<std::map<int,std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,bool> > > > > d_dirichlet_bdry;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredPhysicalBoundaryHelper.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredPhysicalBoundaryHelper
