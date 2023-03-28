// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_INSCollocatedConvectiveOperatorManager
#define included_IBAMR_INSCollocatedConvectiveOperatorManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSCollocatedConvectiveOperatorManager is a singleton manager
 * class to provide access to generic cell-centered ConvectiveOperator
 * implementations suitable for use with class INSCollocatedHierarchyIntegrator.
 */
class INSCollocatedConvectiveOperatorManager
{
public:
    /*!
     * Default ConvectiveOperator types automatically provided by the manager
     * class.
     */
    static const std::string DEFAULT;
    static const std::string CENTERED;
    static const std::string PPM;
    static const std::string WAVE_PROP;

    /*!
     * Return a pointer to the instance of the operator manager.  Access to
     * INSCollocatedConvectiveOperatorManager objects is mediated by the getManager()
     * function.
     *
     * \return A pointer to the operator manager instance.
     */
    static INSCollocatedConvectiveOperatorManager* getManager();

    /*!
     * Deallocate the INSCollocatedConvectiveOperatorManager instance.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeManager();

    /*!
     * Allocate a new INSCollocatedConvectiveOperator object of the specified type.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocateOperator(const std::string& operator_type,
                     const std::string& operator_object_name,
                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                     ConvectiveDifferencingType difference_form,
                     const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs) const;

    /*!
     * Typedef for functions to construct cell-centered ConvectiveOperators.
     */
    using OperatorMaker = SAMRAI::tbox::Pointer<ConvectiveOperator> (*)(
        const std::string& operator_object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        ConvectiveDifferencingType difference_form,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * Register a operator factory function with the operator manager class.
     */
    void registerOperatorFactoryFunction(const std::string& operator_type, OperatorMaker operator_maker);

protected:
    /*!
     * \brief Default constructor.
     */
    INSCollocatedConvectiveOperatorManager();

    /*!
     * \brief Destructor.
     */
    ~INSCollocatedConvectiveOperatorManager() = default;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSCollocatedConvectiveOperatorManager(const INSCollocatedConvectiveOperatorManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCollocatedConvectiveOperatorManager& operator=(const INSCollocatedConvectiveOperatorManager& that) = delete;

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static INSCollocatedConvectiveOperatorManager* s_operator_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Mapping from operator type names to operator maker functions.
     */
    std::map<std::string, OperatorMaker> d_operator_maker_map;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSCollocatedConvectiveOperatorManager
