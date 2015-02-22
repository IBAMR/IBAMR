
#ifndef LIBMESH_QANISOTROPIC_COMPOSITE_BASE_H
#define LIBMESH_QANISOTROPIC_COMPOSITE_BASE_H

// Local includes
#include "ibtk/QuadratureAnisotropicBase.h"
// C++ includes

namespace libMesh
{


// ------------------------------------------------------------
// QuadratureAnisotropicCompositeBase class definition

class QuadratureAnisotropicCompositeBase : public QuadratureAnisotropicBase
{
 public:

  /**
   * Constructor.
   */
  QuadratureAnisotropicCompositeBase (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicCompositeBase (const unsigned int _dim,
	 const std::vector<unsigned int >& vec_interval_num, bool use_anisotropic =true);
  
  QuadratureAnisotropicCompositeBase (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_interval_num, const unsigned int nqp_per_interval);  
  
  QuadratureAnisotropicCompositeBase (const unsigned int _dim,
	 const unsigned int interval_num, const unsigned int nqp_per_interval); //valid for 1D case.
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicCompositeBase();

  /**
   * @returns \p QuadratureAnisotropicCompositeBase
   */
  //QuadratureType type() const { return QuadratureAnisotropicCompositeBase; }

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
QuadratureAnisotropicCompositeBase::QuadratureAnisotropicCompositeBase(const unsigned int d,
	     const unsigned int o) 
	  : QuadratureAnisotropicBase(d,o), d_use_composite(false),d_num_qps(0)
{
 
 std::cout<< " warning: Anisotropic/Composite needs 3 int-type arguments, rule is not right" << std::endl; 
}



inline
QuadratureAnisotropicCompositeBase::QuadratureAnisotropicCompositeBase(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_adaptive) 
	: QuadratureAnisotropicBase(d,o,use_adaptive),d_use_composite(false),d_num_qps(o[0]) // Composite rule always ask for 3 ints (1 int, 1 vec[int], and 1 int)
{ 
  
  
}
inline
QuadratureAnisotropicCompositeBase::QuadratureAnisotropicCompositeBase (const unsigned int d,
	 const std::vector<unsigned int> & vec_interval_num, const unsigned int nqp_per_interval)
	:QuadratureAnisotropicBase(d,vec_interval_num,true) // this will put: Order = vec_interval_num[0], which will be modified based on d_num_qps in specific rules
	,d_use_composite(true), d_num_qps(nqp_per_interval)
{
  
}
inline
QuadratureAnisotropicCompositeBase::QuadratureAnisotropicCompositeBase (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval)
	 :QuadratureAnisotropicBase(d,interval_num),d_use_composite(true),d_num_qps(nqp_per_interval)//_order = interval_num; d_vector_order(d, interval_num)
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
QuadratureAnisotropicCompositeBase::tranform_1D( std::vector<Real> & vec_weights, std::vector<Point> & vec_points, 
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
QuadratureAnisotropicCompositeBase::~QuadratureAnisotropicCompositeBase()
{
}


}// namespace libMesh


#endif // LIBMESH_QUADRATURE_Composite_H
