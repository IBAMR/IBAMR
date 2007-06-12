#ifndef included_IBInstrumentPanel
#define included_IBInstrumentPanel

// Filename: IBInstrumentPanel.h
// Last modified: <12.Jun.2007 01:55:30 boyce@bigboy.nyconnect.com>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LDataManager.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <Index.h>
#include <PatchHierarchy.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInstrumentPanel provides support for flow meters and pressure
 * gauges.
 */
class IBInstrumentPanel
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBInstrumentPanel(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentPanel();

    /*!
     * \brief Compute the flow rates and pressures in the various distributed
     * internal meters.
     */
    void
    readMeters(
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_var,
        const int U_data_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
        const int P_data_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager);

    /*!
     * \brief Compute the positions of all of the flow meters and pressure
     * gauges.
     */
    void
    initializeMeterData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager);

    /*!
     * \brief Initialize hierarchy-dependent data.
     */
    void
    initializeHierarchyData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInstrumentPanel(
        const IBInstrumentPanel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentPanel&
    operator=(
        const IBInstrumentPanel& that);

    /*!
     * \brief Hierarchy independent data.
     */
    int d_interp_level, d_num_meters;
    std::vector<int> d_num_perimeter_nodes;

    /*!
     * \brief Hierarchy and configuration dependent data.
     */
    std::vector<std::vector<double> > d_X_perimeter, d_X_centroid;
    std::vector<std::vector<std::vector<double> > > d_X_web, d_dA_web;

    /*!
     * \brief Data structures employed to manage the mapping between web patch
     * data (i.e., centroids and area-weigthed normals) and cell indices.
     */
    struct IndexFortranOrder
        : public std::binary_function<SAMRAI::hier::Index<NDIM>,SAMRAI::hier::Index<NDIM>,bool>
    {
        inline bool
        operator()(
            const SAMRAI::hier::Index<NDIM>& lhs,
            const SAMRAI::hier::Index<NDIM>& rhs) const
            {

                return (lhs(0) < rhs(0)
#if (NDIM>1)
                        || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM>2)
                        || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                        );
            }
    };

    struct WebPatch
    {
        int meter_num;
        double* const X;
        double* const dA;
    };

    typedef std::map<SAMRAI::hier::Index<NDIM>,WebPatch,IndexFortranOrder> WebPatchMap;
    std::vector<WebPatchMap> d_web_patch_map, d_centroid_map;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBInstrumentPanel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentPanel
