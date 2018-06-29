// Filename LSLocateGasInterface.h
// Created on Nov 15, 2017 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateGasInterface
#define included_LSLocateGasInterface

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

void callLSLocateGasInterfaceCallbackFunction(int D_idx,
                                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                              double time,
                                              bool initial_time,
                                              void* ctx);

class LSLocateGasInterface
{
    /*!
     * \brief class LSLocateGasInterface is a utility class which is used to identify
     * the circular interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateGasInterface(const std::string& object_name,
                         SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                         const double init_height);

    /*!
     * Destructor for this class.
     */
    ~LSLocateGasInterface();

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
    LSLocateGasInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateGasInterface& operator=(const LSLocateGasInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateGasInterface(const LSLocateGasInterface& from);

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
    double d_init_height;
};

#endif // #ifndef included_LSLocateGasInterface
