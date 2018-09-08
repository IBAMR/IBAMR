// Filename LSLocateFluidInterface.h
// Created on Jul 5, 2018 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateFluidInterface
#define included_LSLocateFluidInterface

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

void callLSLocateFluidInterfaceCallbackFunction(int D_idx,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                double time,
                                                bool initial_time,
                                                void* ctx);

// Struct to maintain the height of the thin liquid film
struct FilmInterface
{
    double height;
};

// Struct to maintain the bubble interface
struct CircularInterface
{
    IBTK::Vector X0;
    double R;
};

class LSLocateFluidInterface
{
    /*!
     * \brief class LSLocateFluidInterface is a utility class which is used to identify
     * the fluid
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateFluidInterface(const std::string& object_name,
                           SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                           FilmInterface init_film,
                           CircularInterface init_circle);

    /*!
     * Destructor for this class.
     */
    ~LSLocateFluidInterface();

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
    LSLocateFluidInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateFluidInterface& operator=(const LSLocateFluidInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateFluidInterface(const LSLocateFluidInterface& from);

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
    FilmInterface d_init_film;
    CircularInterface d_init_circle;
};

#endif // #ifndef included_LSLocateFluidInterface
