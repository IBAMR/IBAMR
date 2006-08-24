#ifndef included_QInit
#define included_QInit

// Filename: QInit.h
// Last modified: <24.Aug.2006 00:22:39 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <tbox/Database.h>
#include <GridGeometry.h>
#include <Patch.h>
#include <tbox/Pointer.h>
#include <Variable.h>
#include <tbox/Array.h>

// C++ STDLIB INCLUDES
#include <string>

// NAMESPACE
using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * @brief Class to initialize the value of the advected scalar Q.
 */
class QInit
    : public SetDataStrategy
{
public:
    /*!
     * @brief Constructor.
     */
    QInit(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);
    
    /*!
     * @brief Destructor.
     */
    ~QInit();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    bool isTimeDependent() const { return false; }
    
    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);
    
protected:
    
private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    QInit(
        const QInit& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    QInit& operator=(
        const QInit& that);

    /*!
     * Read input values, indicated above, from given database.
     */
    void getFromInput(
        tbox::Pointer<tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The grid geometry.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * The initialization type.
     */
    string d_init_type;

    /*
     * The center of the initial data.
     */
    tbox::Array<double> d_X;
    
    /*
     * Parameters for Gaussian initial conditions.
     */
    double d_gaussian_kappa;

    /*
     * Parameters for the Zalesak slotted cylinder.
     */
    double d_zalesak_r;
    double d_zalesak_slot_w;
    double d_zalesak_slot_l;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifndef DEBUG_NO_INLINE
//#include "QInit.I"
//#endif

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_QInit
