#ifndef included_USet
#define included_USet

// Filename: USet.h
// Last modified: <24.Aug.2006 00:22:46 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <GridGeometry.h>
#include <Patch.h>
#include <Variable.h>
#include <VariableContext.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <string>

// NAMESPACE
using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * @brief Class to initialize the value of the advection velocity u.
 */
class USet
    : public SetDataStrategy
{
public:
    /*!
     * @brief Default constructor.
     */
    USet(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);
    
    /*!
     * @brief Destructor.
     */
    ~USet();
    
    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    bool isTimeDependent() const { return true; }
    
    /*!
     * Set the data on the patch interior to some values.
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
    USet(
        const USet& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    USet& operator=(
        const USet& that);
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
     * The amplification and frequency of the sin wave used in setting
     * velocities.
     */
    tbox::Array<double> d_kappa, d_omega;

    /*
     * Parameters for uniform constant velocity.
     */
    tbox::Array<double> d_uniform_u;    
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifndef DEBUG_NO_INLINE
//#include "USet.I"
//#endif

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_USet
