// Filename: AdvDiffConvectiveOperatorManager.h
// Created on 17 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_AdvDiffConvectiveOperatorManager
#define included_AdvDiffConvectiveOperatorManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <string>
#include <vector>

#include "CellVariable.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

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
 * \brief Class AdvDiffConvectiveOperatorManager is a singleton manager class to
 * provide access to generic cell-centered ConvectiveOperator implementations
 * suitable for use with class AdvDiffHierarchyIntegrator.
 */
class AdvDiffConvectiveOperatorManager
{
public:
    /*!
     * Default ConvectiveOperator types automatically provided by the manager
     * class.
     */
    static const std::string DEFAULT;
    static const std::string CENTERED;
    static const std::string PPM;

    /*!
     * Return a pointer to the instance of the operator manager.  Access to
     * AdvDiffConvectiveOperatorManager objects is mediated by the getManager()
     * function.
     *
     * \return A pointer to the operator manager instance.
     */
    static AdvDiffConvectiveOperatorManager* getManager();

    /*!
     * Deallocate the AdvDiffConvectiveOperatorManager instance.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeManager();

    /*!
     * Allocate a new AdvDiffConvectiveOperator object of the specified type.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocateOperator(const std::string& operator_type,
                     const std::string& operator_object_name,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                     ConvectiveDifferencingType difference_form,
                     const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs) const;

    /*!
     * Typedef for functions to construct cell-centered ConvectiveOperators.
     */
    typedef SAMRAI::tbox::Pointer<ConvectiveOperator> (*OperatorMaker)(
        const std::string& operator_object_name,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
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
    AdvDiffConvectiveOperatorManager();

    /*!
     * \brief Destructor.
     */
    ~AdvDiffConvectiveOperatorManager();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffConvectiveOperatorManager(const AdvDiffConvectiveOperatorManager& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffConvectiveOperatorManager& operator=(const AdvDiffConvectiveOperatorManager& that);

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static AdvDiffConvectiveOperatorManager* s_operator_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * Mapping from operator type names to operator maker functions.
     */
    std::map<std::string, OperatorMaker> d_operator_maker_map;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffConvectiveOperatorManager
