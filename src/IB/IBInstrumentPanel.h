#ifndef included_IBInstrumentPanel
#define included_IBInstrumentPanel

// Filename: IBInstrumentPanel.h
// Last modified: <13.Jun.2007 15:08:35 griffith@box221.cims.nyu.edu>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LDataManager.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <Index.h>
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>

// BLITZ++ INLCUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInstrumentPanel provides support for flow meters and pressure
 * gauges.
 *
 * \note Use of class IBInstrumentPanel requires the Blitz++ array library.
 */
class IBInstrumentPanel
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name          String used for error reporting.
     * \param dump_directory_name  String indicating the directory where visualization data is to be written.
     */
    IBInstrumentPanel(
        const std::string& object_name,
        const std::string& dump_directory_name);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentPanel();

    /*!
     * \return A const reference to the vector of instrument names.
     */
    const std::vector<std::string>&
    getInstrumentNames() const;

    /*!
     * \return A const reference to the most recent time that the flow meter and
     * pressure gauge values were set.
     *
     * \note This value is not initialized until the first call is made to
     * readInstrumentData().
     */
    const double&
    getInstrumentDataReadTime() const;

    /*!
     * \return A const reference to the vector of flow meter values.
     *
     * \note This vector is not initialized until the first call is made to
     * readInstrumentData().
     */
    const std::vector<double>&
    getFlowValues() const;

    /*!
     * \return A const reference to the vector of pressure gauge values.
     *
     * \note This vector is not initialized until the first call is made to
     * readInstrumentData().
     */
    const std::vector<double>&
    getPressureValues() const;

    /*!
     * \brief Initialize hierarchy-independent data.
     *
     * The data initialized by this method is assumed \em not to change during
     * the course of a simulation.
     */
    void
    initializeHierarchyIndependentData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager);

    /*!
     * \brief Initialize hierarchy- and configuration-dependent data.
     */
    void
    initializeHierarchyDependentData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager);

    /*!
     * \brief Compute the flow rates and pressures in the various distributed
     * internal flow meters and pressure gauges.
     */
    void
    readInstrumentData(
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_var,
        const int U_data_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
        const int P_data_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager,
        const double data_time);

    /*!
     * \brief Write the plot data to disk.
     */
    void
    writePlotData(
        const int time_step_number,
        const double simulation_time);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBInstrumentPanel();

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

    /*
     * The object name is used for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * The directory where data is to be dumped.
     */
    std::string d_dump_directory_name;

    /*
     * Time step number (passed in by user).
     */
    int d_time_step_number;

    /*!
     * \brief Instrumentation data.
     */
    int d_num_meters;
    std::vector<int> d_num_perimeter_nodes;
    std::vector<blitz::TinyVector<double,NDIM> > d_X_centroid;
    std::vector<blitz::Array<blitz::TinyVector<double,NDIM>,1> > d_X_perimeter;
    std::vector<blitz::Array<blitz::TinyVector<double,NDIM>,2> > d_X_web, d_dA_web;

    double d_instrument_read_time;
    std::vector<std::string> d_instrument_names;
    std::vector<double> d_flow_values, d_pressure_values;

    /*!
     * \brief Data structures employed to manage mappings between cell indices
     * and web patch data (i.e., patch centroids and area-weigthed normals) and
     * meter centroid data.
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
            }// operator()
    };

    struct WebPatch
    {
        int meter_num;
        const blitz::TinyVector<double,NDIM>* X;
        const blitz::TinyVector<double,NDIM>* dA;
    };

    typedef std::multimap<SAMRAI::hier::Index<NDIM>,WebPatch,IndexFortranOrder> WebPatchMap;
    std::vector<WebPatchMap> d_web_patch_map;

    struct MeterCentroid
    {
        int meter_num;
        const blitz::TinyVector<double,NDIM>* X;
    };

    typedef std::multimap<SAMRAI::hier::Index<NDIM>,MeterCentroid,IndexFortranOrder> MeterCentroidMap;
    std::vector<MeterCentroidMap> d_meter_centroid_map;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBInstrumentPanel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentPanel
