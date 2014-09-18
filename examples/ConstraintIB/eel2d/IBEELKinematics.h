//Filename : IBEELKinematics.h
//Created by Amneet Bhalla on 6/16/2011.

// Example case taken from:

//  *  Bhalla et al. A unified mathematical framework and an adaptive numerical method for 
//     fluid-structure interaction with rigid, deforming, and elastic bodies. J Comput Phys, 250:446-476 (2013). 


#ifndef included_ibeelkinematics
#define included_ibeelkinematics

/////////////////////////////////////// INCLUDES ////////////////////////////////
#include <iostream>
#include <vector>
#include <map>

#include "ibamr/ConstraintIBKinematics.h"  
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "PatchHierarchy.h"


namespace mu
{
  class Parser;
}// namespace mu


///////////////////////////////////////////////////////////////// CLASS DEFORMATIONAL KINEMATICS //////////////////

namespace IBAMR
{
  /*!
   * \brief IBEELKinematics Class.
   * 
   * IBEELKinematics is a concrete class which calculates the deformation velocity and updated shape
   * for 2D eel. It also provides routines for maneuvering and food tracking cases.
   */
  
class IBEELKinematics 
   : public ConstraintIBKinematics
   
{

public: 
      
    /*!
     * \brief ctor. This is the only ctor for this object.
     */
    IBEELKinematics(        
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        IBTK::LDataManager* l_data_manager,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
        bool register_for_restart = true);
      
    /*!
     * \brief Destructor.
     */
    virtual
    ~IBEELKinematics();
      
    /*!
     * \brief Set kinematics velocity at new time for the structure.
     */
    virtual void
    setNewKinematicsVelocity(
        const double new_time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);
    
    /*!
     * \brief Get the kinematics velocity at new time for the structure on the specified level. 
     */
    virtual const std::vector<std::vector<double> >&
    getNewKinematicsVelocity(
        const int level) const;
    
    /*!
     * \brief Get the kinematics velocity at current time for the structure on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getCurrentKinematicsVelocity(
        const int level) const;
  
    /*!
     * \brief Set the shape of the structure at new time on all levels.
     */
    virtual void
    setNewShape();
    
    /*!
     * \brief Get the shape of structure at new time  on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getNewShape(
        const int level) const;

    /*!
     * \brief Override the ConstraintIBkinematics base class method.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
  
    /*!
     * \brief The default constructor is not implemented and should not be used.
     */
    IBEELKinematics();
      
    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBEELKinematics(
        const IBEELKinematics&from);
      
    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBEELKinematics&
    operator = (
        const IBEELKinematics& that);
    
    /*!
     * \brief Set data from restart.
     */
    void
    getFromRestart();
      
    /*!
     * \brief set eel body shape related data.
     */
    void
    setImmersedBodyLayout(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);
    
    /*!
     * \brief Set deformation kinematics velocity of the eel.
     */
    void
    setEelSpecificVelocity(
        const double time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);
 
    /*!
     * \brief Rotate the maneuver axis and caluclate tangents in this orientation.
     */
    void
    transformManeuverAxisAndCalculateTangents(
        const double angleFromHorizontal);
    
    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;
    
    /*!
     * Current and new deformational velocity
     */
    std::vector<std::vector<double> >   d_current_kinematics_vel, d_new_kinematics_vel; 
    std::vector<std::vector<double> > d_new_shape;
    
    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;
    
    /*!
     * Eulerian Mesh width parameters.
     */
    std::vector<double> d_mesh_width;
    
     /*!
      * The following map is used to store eel body shape specific data.
      * The arc length 's' varies from 0 - 1. In the std::map  the arc length 's' is used as a key.
      * d_ImmersedBodyData is used to store the no of material points which represents a cross section. The 
      * width of cross section of eel varies with arc length. 
      */
     std::map<double, int> d_ImmersedBodyData;
     
     /*!
      * Initial orientation of the body axis
      */
     double d_initAngle_bodyAxis_x;
         
     /*!
      * Boolean value indicating if eel is maneuvering or not.
      * If eel is maneuvering then the traveling wave will be along a curved axis, otherwise it will be on a straight line
      */
     bool d_bodyIsManeuvering;
     
     /*!
      * Boolean to indicate if shape of maneuver axis is changing. The maneuver axis will change shape in food tracking cases.
      */
     bool d_maneuverAxisIsChangingShape;
     
     /*!
      * Vector of coordinates defining the axis of maneuvering. The reference axis will rotate with body omega.
      * The coordinates of the rotated maneuver axis is stored separately.
      */
     std::vector< std::vector<double> > d_maneuverAxisReferenceCoordinates_vec, d_maneuverAxisTransformedCoordinates_vec;
     
     /*!
      * map of tangents along the body/maneuver axis in rotated frame. The key used is arc length 's' and it stores only the abs(theta).
      * Sign of tangent is stored separately. This is done to avoid a lot of if conditions needed to determine the quadrant of the
      * angle.
      */
     std::map<double,double> d_map_transformed_tangent;
     
      /*!
      * map of tangents along the body/maneuver axis in reference/unrotated frame.The key used is arc length 's' and it stores only the abs(theta).
      * Sign of tangent is stored separately. This is done to avoid a lot of if conditions needed to determine the quadrant of the
      * angle.
      */
     std::map<double,double> d_map_reference_tangent;
     
     /*!
      * Sign of tangent vector in rotated frame. The key used is arc length 's'. 'mapped_value' is a vector which has sign of t_x and t_y 
      * respectively.
      */
     std::map<double, std::vector<int> > d_map_transformed_sign;
     
     /*!
      * Sign of tangent vector in reference/unrotated frame. The key used is arc length 's'. 'mapped_value' is a vector which has sign of t_x and t_y 
      * respectively.
      */
     std::map<double, std::vector<int> > d_map_reference_sign;
     
     /*!
      * mu::Parser object which evaluates the maneuvering axis equation.
      */
     mu::Parser* d_maneuvering_axis_parser;
     
     /*!
      * mu::Parser object which evaluates the shape of the body.
      */
     mu::Parser* d_body_shape_parser;
     
     /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
     std::vector<mu::Parser*> d_deformationvel_parsers;
     std::vector<mu::Parser*> d_all_parsers;
     
    /*!
     * Time and position variables.
     */
    double* d_parser_time;
    double* d_parser_posn;
    double* d_parser_normal;
    
    /*!
     * Array containing initial coordinates of the food location.
     */
    SAMRAI::tbox::Array<double> d_food_location;         
    
};// IBEELKinematics
  
}// IBAMR
#endif //#ifndef included_ibeelkinematics
