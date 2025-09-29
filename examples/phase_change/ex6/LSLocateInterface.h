// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateInterface
#define included_LSLocateInterface

///////////////////////////// INCLUDES ///////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*
 * Pre processing call back function to be hooked into IBAMR:LSInitStrategy
 */

void callLSLocateInterfaceCallbackFunction(int D_idx,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           double time,
                                           bool initial_time,
                                           void* ctx);

class LSLocateInterface
{
    /*!
     * \brief class LSLocateInterface is a utility class which is used to identify
     * the interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateInterface(const std::string& object_name,
                      SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                      const double initial_horizontal_interface_position,
                      std::vector<std::pair<double, IBTK::Vector> > bubbles_position);

    /*!
     * Destructor for this class.
     */
    ~LSLocateInterface();

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double time,
                              const bool initial_time);

private:
    /*!
     * Deleted default constructor.
     */
    LSLocateInterface() = delete;

    /*!
     * Deleted copy constructor.
     */
    LSLocateInterface(const LSLocateInterface& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LSLocateInterface& operator=(const LSLocateInterface& that) = delete;

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
     * Initial position of a horizontal interface.
     */
    const double d_initial_horizontal_interface_position;

    /*!
     * Initial position and radii of bubbles.
     */
    std::vector<std::pair<double, IBTK::Vector> > d_bubbles_position;
};

#endif // #ifndef included_LSLocateInterface
