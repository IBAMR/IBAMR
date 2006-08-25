#ifndef included_ConvergenceMonitor
#define included_ConvergenceMonitor

// Filename: ConvergenceMonitor.h
// Last modified: <25.Aug.2006 00:58:57 boyce@bigboy.nyconnect.com>
// Created on 19 Jun 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <Variable.h>
#include <VariableContext.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <set>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Simple class to manage convergence monitoring.  Note that
 * this class should be redesigned in the near future.
 */
class ConvergenceMonitor
{
public:
    /*!
     * @brief Constructor.
     */
    ConvergenceMonitor(
        const std::string& object_name);
    
    /*!
     * @brief Empty destructor.
     */
    ~ConvergenceMonitor();

    /*!
     * \todo Write docs.
     */
    void registerMonitoredVariableAndContext(
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> var_ctx,
        SAMRAI::tbox::Pointer<SetDataStrategy> exact_soln_setter);

    /*!
     * \todo Write docs.
     */
    void monitorConvergence(
        const double data_time);
    
    /*!
     * \todo Write docs.
     */
    void initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * \todo Write docs.
     */
    void resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);
    
private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    ConvergenceMonitor();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    ConvergenceMonitor(
        const ConvergenceMonitor& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    ConvergenceMonitor& operator=(
        const ConvergenceMonitor& that);

    /*
     * The object name is used for error reporting purposes.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * SAMRAI::tbox::Pointer to the patch hierarchy this monitor object uses.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * The cell weights are used to compute norms of data defined on
     * the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_var;
    int d_wgt_idx;

    /*
     * The patch data descriptor indices and affiliated Variables,
     * VariableContexts, and exact solution data setters for the
     * variables being monitored by this object.
     */
    std::set<int> d_monitored_indices;
    std::map<int,SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_monitored_vars;
    std::map<int,SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> > d_monitored_var_ctxs;
    std::map<int,SAMRAI::tbox::Pointer<SetDataStrategy> > d_exact_soln_setters;
};
}

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "ConvergenceMonitor.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ConvergenceMonitor
