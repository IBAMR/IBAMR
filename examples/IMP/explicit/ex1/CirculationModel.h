// Filename: CirculationModel.h
// Created on 04 May 2007 by Boyce Griffith

#ifndef included_CirculationModel
#define included_CirculationModel

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel
 */
class CirculationModel
    : public CartGridFunction
{
public:
    /*!
     * \brief Constructor
     */
    CirculationModel(
        const string& object_name,
        Pointer<Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~CirculationModel();

    /*!
     * \brief Advance time-dependent data.
     */
    void
    advanceTimeDependentData(
        const double data_time,
        const double dt,
        const Pointer<PatchHierarchy<NDIM> > hierarchy,
        const int U_idx,
        const int P_idx,
        const int wgt_cc_idx,
        const int wgt_sc_idx);

    /*!
     * \name Implementation of CartGridFunction interface.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool
    isTimeDependent() const;

    /*!
     * \brief Set data on the specified patch interior.
     */
    void
    setDataOnPatch(
        int data_idx,
        Pointer<hier::Variable<NDIM> > var,
        Pointer<Patch<NDIM> > patch,
        double data_time,
        bool initial_time=false,
        Pointer<PatchLevel<NDIM> > patch_level=Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CirculationModel();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel(
        const CirculationModel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel&
    operator=(
        const CirculationModel& that);

    /*
     * Mean pressure and resulting flow rate.
     */
    double d_P_mean, d_radius;
    double d_Q;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <CirculationModel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CirculationModel
