#ifndef included_IBLagrangianSourceStrategy
#define included_IBLagrangianSourceStrategy

// Filename: IBLagrangianSourceStrategy.h
// Last modified: <23.Jan.2007 21:23:10 griffith@box221.cims.nyu.edu>
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
 * @brief Class IBLagrangianSourceStrategy provides a mechanism for
 * specifying fluid source/sinks at arbitrary point locations.
 */
class IBLagrangianSourceStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * @brief Default constructor.
     */
    IBLagrangianSourceStrategy();

    /*!
     * @brief Destructor.
     */
    virtual ~IBLagrangianSourceStrategy();

    virtual int getNumSources(
        const int level_number) const = 0;

    virtual void getSourceLocations(
        std::vector<std::vector<double> >& X_src,
        vector<double>& r_src,
        SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        const LDataManager* const lag_manager) = 0;

    virtual void computeSourceStrengths(
        std::vector<double>& Q_src,
        const std::vector<double>& P_src,
        const std::vector<std::vector<double> >& X_src,
        const vector<double>& r_src,
        const int level_number,
        const double data_time) = 0;

    /*!
     * @brief Setup the data needed to compute source/sink data on the
     * specified level of the patch hierarchy.
     *
     * NOTE: A default empty implementation is provided.
     */
    virtual void initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool initial_time,
        const LDataManager* const lag_manager);

private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    IBLagrangianSourceStrategy(
        const IBLagrangianSourceStrategy& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    IBLagrangianSourceStrategy& operator=(
        const IBLagrangianSourceStrategy& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBLagrangianSourceStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBLagrangianSourceStrategy
