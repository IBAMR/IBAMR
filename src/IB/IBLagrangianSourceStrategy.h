#ifndef included_IBLagrangianSourceStrategy
#define included_IBLagrangianSourceStrategy

// Filename: IBLagrangianSourceStrategy.h
// Last modified: <05.Apr.2007 18:44:03 griffith@box221.cims.nyu.edu>
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LDataManager.h>
#include <ibamr/LNodeLevelData.h>

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
 * \brief Class IBLagrangianSourceStrategy provides a mechanism for
 * specifying fluid source/sinks at arbitrary point locations.
 */
class IBLagrangianSourceStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBLagrangianSourceStrategy();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBLagrangianSourceStrategy();

    /*!
     * \brief Set the current and new times for the present timestep.
     *
     * \note A default empty implementation is provided.
     */
    virtual void setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \brief Setup the data needed to compute source/sink data on the
     * specified level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool initial_time,
        LDataManager* const lag_manager);

    /*!
     * \brief Specify the number of distributed internal sources or
     * sinks.
     *
     * \note The return value must be the \em total number of internal
     * sources/sinks in the \em entire computational domain.  This
     * implies that the return value must be \em identical on each MPI
     * process.
     */
    virtual int getNumSources(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        LDataManager* const lag_manager) = 0;

    /*!
     * \brief Compute the source locations for each of the distributed
     * internal sources or sinks.
     *
     * \note Implementations of this method \em must compute the same
     * values for \a X_src on \em each MPI process.  That is to say,
     * \a X_src must provide the location of all of the distributed
     * sources/sinks.
     */
    virtual void getSourceLocations(
        std::vector<std::vector<double> >& X_src,
        vector<double>& r_src,
        SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        LDataManager* const lag_manager) = 0;

    /*!
     * \brief Set the normalized pressures at the sources.
     */
    virtual void setSourcePressures(
        const std::vector<double>& P_src,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        LDataManager* const lag_manager) = 0;

    /*!
     * \brief Compute the source strengths for each of the distributed
     * internal sources or sinks.
     *
     * \note Implementations of this method \em must compute the same
     * values for \a Q_src on \em each MPI process.  That is to say,
     * \a Q_src must provide the strengths of all of the distributed
     * sources/sinks.
     */
    virtual void computeSourceStrengths(
        std::vector<double>& Q_src,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        LDataManager* const lag_manager) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IBLagrangianSourceStrategy(
        const IBLagrangianSourceStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBLagrangianSourceStrategy& operator=(
        const IBLagrangianSourceStrategy& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBLagrangianSourceStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBLagrangianSourceStrategy
