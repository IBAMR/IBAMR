// Filename LSLocateLayerInterface.h
// Created on Jul 5, 2018 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateLayerInterface
#define included_LSLocateLayerInterface

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

void callLSLocateLayerInterfaceCallbackFunction(int D_idx,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                double time,
                                                bool initial_time,
                                                void* ctx);

// Struct to maintain the properties of the fluid layer interface, in terms of y- (z-) coordinate.
struct LayerInterface
{
    double height;
};

class LSLocateLayerInterface
{
    /*!
     * \brief class LSLocateLayerInterface is a utility class which is used to identify
     * the circular interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateLayerInterface(const std::string& object_name,
                           SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                           LayerInterface init_layer);

    /*!
     * Destructor for this class.
     */
    ~LSLocateLayerInterface();

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
    LSLocateLayerInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateLayerInterface& operator=(const LSLocateLayerInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateLayerInterface(const LSLocateLayerInterface& from);

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
    LayerInterface d_init_layer;
};

#endif // #ifndef included_LSLocateLayerInterface
