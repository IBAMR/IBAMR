// Filename LSLocateColumnInterface.h
// Created on Jul 5, 2018 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateColumnInterface
#define included_LSLocateColumnInterface

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

void callLSLocateColumnInterfaceCallbackFunction(int D_idx,
                                                 SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                 double time,
                                                 bool initial_time,
                                                 void* ctx);

// Struct to maintain the properties of the column interface, in terms the upper right corner
// Assumes the lower left corner of the column is the lower left corner of the domain, for simplicity
struct ColumnInterface
{
    IBTK::Vector X_UR;
};

class LSLocateColumnInterface
{
    /*!
     * \brief class LSLocateColumnInterface is a utility class which is used to identify
     * the circular interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateColumnInterface(const std::string& object_name,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                            ColumnInterface init_column);

    /*!
     * Destructor for this class.
     */
    ~LSLocateColumnInterface();

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
    LSLocateColumnInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateColumnInterface& operator=(const LSLocateColumnInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateColumnInterface(const LSLocateColumnInterface& from);

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
    ColumnInterface d_init_column;
};

#endif // #ifndef included_LSLocateColumnInterface
