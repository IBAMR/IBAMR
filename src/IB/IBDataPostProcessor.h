#ifndef included_IBDataPostProcessor
#define included_IBDataPostProcessor

// Filename: IBDataPostProcessor.h
// Last modified: <24.Sep.2008 17:32:46 griffith@box230.cims.nyu.edu>
// Created on 24 Sep 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeLevelData.h>
#include <ibtk/LDataManager.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBDataPostProcessor provides a generic interface for specifying
 * post-processing code for use in an IB computation.
 */
class IBDataPostProcessor
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBDataPostProcessor();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBDataPostProcessor();

    /*!
     * \brief Setup the data needed to post-process the data on the specified
     * level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool initial_time,
        IBTK::LDataManager* const lag_manager);

    /*!
     * \brief Post-process data on the patch hierarchy.
     */
    virtual void
    postProcessData(
        const int u_idx,
        const int p_idx,
        const int f_idx,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > F_data,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > X_data,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > U_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level_number,
        const int finest_level_number,
        const double data_time,
        IBTK::LDataManager* const lag_manager) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBDataPostProcessor(
        const IBDataPostProcessor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBDataPostProcessor&
    operator=(
        const IBDataPostProcessor& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBDataPostProcessor.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBDataPostProcessor
