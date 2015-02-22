


#ifndef included_QuadratureAnisotropicCompositeNewtonCotes
#define included_QuadratureAnisotropicCompositeNewtonCotes

// Local includes

#include "ibtk/QuadratureAnisotropicCompositeBase.h"
// C++ includes

namespace libMesh
{



// ------------------------------------------------------------
// QuadratureAnisotropicCompositeGauss class definition

class QuadratureAnisotropicCompositeNewtonCotes : public QuadratureAnisotropicCompositeBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QuadratureAnisotropicCompositeNewtonCotes (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicCompositeNewtonCotes (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_order, bool use_anisotropic =true);
  
  QuadratureAnisotropicCompositeNewtonCotes (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_interval_num, const unsigned int nqp_per_interval);
  
  QuadratureAnisotropicCompositeNewtonCotes (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval); // only for 1D
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicCompositeNewtonCotes();

  /**
   * @returns \p QuadratureAnisotropicGRID
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_COMPOSITE_NEWTON_COTES); }


 private:
   
 
  void buildStandard1D( std::vector<Real>& standard_weights,  std::vector<Real>& standard_cods);
  
  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  

};



// ------------------------------------------------------------

inline
QuadratureAnisotropicCompositeNewtonCotes::QuadratureAnisotropicCompositeNewtonCotes(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicCompositeBase(d,o)
{
}

inline
QuadratureAnisotropicCompositeNewtonCotes::QuadratureAnisotropicCompositeNewtonCotes(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicCompositeBase(d,o,use_anisotropic)
{
  
  
}

inline
QuadratureAnisotropicCompositeNewtonCotes::QuadratureAnisotropicCompositeNewtonCotes (const unsigned int d,
	 const std::vector<unsigned int> & vec_interval_num, const unsigned int nqp_per_interval)
	 :QuadratureAnisotropicCompositeBase(d, vec_interval_num, nqp_per_interval)
{
  
}

inline
 QuadratureAnisotropicCompositeNewtonCotes::QuadratureAnisotropicCompositeNewtonCotes (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval) //only for 1D case
 :QuadratureAnisotropicCompositeBase(d,interval_num,nqp_per_interval)
 
 {
 }
 
 
inline
QuadratureAnisotropicCompositeNewtonCotes::~QuadratureAnisotropicCompositeNewtonCotes()
{
}


} // namespace libMesh



#endif // included_QuadratureAnisotropicCompositeNewtonCotes 
