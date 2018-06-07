// Filename: SetLSProperties.h
// Created on May 24, 2018 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetLSProperties
#define included_SetLSProperties

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/LSInitStrategy.h>
#include <ibtk/ibtk_utilities.h>
#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back functionS to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param ls_idx a patch data index for the current level set variable maintained by the integrator.
 * \param ctx is the pointer to SetLSProperties class object.
 */

void callSetSolidLSCallbackFunction(int ls_solid_idx,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int integrator_step,
                                    const double current_time,
                                    const bool initial_time,
                                    const bool regrid_time,
                                    void* ctx);

void callSetGasLSCallbackFunction(int ls_gas_idx,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                  const int integrator_step,
                                  const double current_time,
                                  const bool initial_time,
                                  const bool regrid_time,
                                  void* ctx);
class SetLSProperties
{
    /*!
     * \brief Class SetLSProperties is a utility class which sets level set values on the patch hierarchy
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetLSProperties(const std::string& object_name,
                    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_solid_ops,
                    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_gas_ops);

    /*!
     * Destructor for this class.
     */
    ~SetLSProperties();

    /*!
     * Set the density based on the current solid level set information
     */
    void setLSSolidPatchData(int ls_solid_idx,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int integrator_step,
                             const double current_time,
                             const bool initial_time,
                             const bool regrid_time);

    /*!
     * Set the density based on the current solid level set information
     */
    void setLSGasPatchData(int ls_gas_idx,
                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                           const int integrator_step,
                           const double current_time,
                           const bool initial_time,
                           const bool regrid_time);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetLSProperties();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetLSProperties& operator=(const SetLSProperties& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetLSProperties(const SetLSProperties& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_solid_ops;
    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_gas_ops;

}; // SetLSProperties

#endif // #ifndef included_SetLSProperties
