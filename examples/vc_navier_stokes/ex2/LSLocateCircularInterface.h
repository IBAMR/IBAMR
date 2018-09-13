// Filename LSLocateCircularInterface.h
// Created on Nov 15, 2017 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateCircularInterface
#define included_LSLocateCircularInterface

///////////////////////////// INCLUDES ///////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>
#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*
 * Pre processing call back function to be hooked into IBAMR:LInitStrategy
 */

void callLSLocateCircularInterfaceCallbackFunction(int D_idx,
                                                   SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                   double time,
                                                   bool initial_time,
                                                   void* ctx);

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    IBTK::Vector X0;
    double R;
};

class LSLocateCircularInterface
{
    /*!
     * \brief class LSLocateCircularInterface is a utility class which is used to identify
     * the circular interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateCircularInterface(const std::string& object_name,
                              SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                              CircularInterface init_circle);

    /*!
     * Destructor for this class.
     */
    ~LSLocateCircularInterface();

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double time,
                              const bool initial_time);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    LSLocateCircularInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateCircularInterface& operator=(const LSLocateCircularInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateCircularInterface(const LSLocateCircularInterface& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to the advection-diffusion solver
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;

    /*!
     * Initial level set information.
     */
    CircularInterface d_init_circle;
};

#endif // #ifndef included_LSLocateCircularInterface