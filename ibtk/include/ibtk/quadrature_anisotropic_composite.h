// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_QANISOTROPIC_COMPOSITE_BASE_H
#define LIBMESH_QANISOTROPIC_COMPOSITE_BASE_H

// Local includes
#include "ibtk/quadrature_anisotropic.h"
// C++ includes

namespace libMesh
{




/**
 * This class creates anisotropic quadrature points on a non-uniform grid, and we also consider composites by integrating subintervals
 * Order = points in each subinterval 
 * only applies to Quad and Hex elements
 */

// ------------------------------------------------------------
// QAnisotropicCompositeBase class definition

class QAnisotropicCompositeBase : public QAnisotropicBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QAnisotropicCompositeBase (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QAnisotropicCompositeBase (const unsigned int _dim,
	 const std::vector<unsigned int >& vec_interval_num, bool use_anisotropic =true);
  
  QAnisotropicCompositeBase (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_interval_num, const unsigned int nqp_per_interval);  
  
  QAnisotropicCompositeBase (const unsigned int _dim,
	 const unsigned int interval_num, const unsigned int nqp_per_interval); //valid for 1D case.
  /**
   * Destructor.
   */
  ~QAnisotropicCompositeBase();

  /**
   * @returns \p QAnisotropicCompositeBase
   */
  //QuadratureType type() const { return QAnisotropicCompositeBase; }

  void tranform_1D( std::vector<Real>& vec_weights, std::vector<Point>& vec_points, unsigned int start_index, //(where to start)
		   const std::vector<Real>& standard_weights, const std::vector<Real>& standard_cods, const double x1, const double x2); 
  
   

protected:
  
   /* tranform xi, wi from standrad_1D [-1,1] to subinerval_ID [x1,x2]
    * 
    */
  bool d_use_composite;
  unsigned int d_num_qps; // per subinterval
  
  // virtual function to build standard_1D in [-1,1]
  
  //void build_standard_1D( std::vector<Real>& standard_weights,  std::vector<Real>& standard_cods){};
  
  // virtual function from QBase class

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0) {}// intentionally do nothing
	
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0) {}
	
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0){}

  
};



// ------------------------------------------------------------
// QGauss class members
inline // should not used;
QAnisotropicCompositeBase::QAnisotropicCompositeBase(const unsigned int d,
	     const unsigned int o) 
	  : QAnisotropicBase(d,o), d_use_composite(false),d_num_qps(0)
{
 
 std::cout<< " warning: Anisotropic/Composite needs 3 int-type arguments, rule is not right" << std::endl; 
}



inline
QAnisotropicCompositeBase::QAnisotropicCompositeBase(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_adaptive) 
	: QAnisotropicBase(d,o,use_adaptive),d_use_composite(false),d_num_qps(o[0]) // Composite rule always ask for 3 ints (1 int, 1 vec[int], and 1 int)
{ 
  
  
}
inline
QAnisotropicCompositeBase::QAnisotropicCompositeBase (const unsigned int d,
	 const std::vector<unsigned int> & vec_interval_num, const unsigned int nqp_per_interval)
	:QAnisotropicBase(d,vec_interval_num,true) // this will put: Order = vec_interval_num[0], which will be modified based on d_num_qps in specific rules
	,d_use_composite(true), d_num_qps(nqp_per_interval)
{
  
}
inline
QAnisotropicCompositeBase::QAnisotropicCompositeBase (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval)
	 :QAnisotropicBase(d,interval_num),d_use_composite(true),d_num_qps(nqp_per_interval)//_order = interval_num; d_vector_order(d, interval_num)
{
  if (d>1) 
  {std::cout<< " error, this is the constructor for 1D case" << std::endl;  libmesh_error();
    return ;
  }
  // else --> this will also set _order = interval_num
  d_use_anisotropic = true;

  return;
  
}


inline
void
QAnisotropicCompositeBase::tranform_1D( std::vector<Real> & vec_weights, std::vector<Point> & vec_points, 
		    const unsigned int start_index,  const std::vector<Real> & standard_weights, 
		    const std::vector<Real> & standard_cods, const double x1, const double x2)
{
  if (vec_weights.size()< standard_weights.size()) 
  { libMesh::err << "ERROR:vec_weights.size: " << vec_weights.size()
    << " < standard_weights.size :" << standard_weights.size() << std::endl;
   libmesh_error();
  return;
  }
  
 
  
  double coef_length = (x2 - x1) / 2.0;
  double coef_center = (x2 + x1) / 2.0;
  
  for (unsigned int k =0; k < standard_weights.size(); k++)
  { 
 
    vec_weights[start_index+k] = standard_weights[k] * coef_length;
    vec_points[start_index+k](0) = standard_cods[k] * coef_length + coef_center;
  }
  
  return;
}



// deconstructor

inline
QAnisotropicCompositeBase::~QAnisotropicCompositeBase()
{
}


}// namespace libMesh


#endif // LIBMESH_QUADRATURE_Composite_H
