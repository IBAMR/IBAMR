


#ifndef included_QuadratureAnisotropicNewtonCotes
#define included_QuadratureAnisotropicNewtonCotes


#include "ibtk/QuadratureAnisotropicBase.h"
// C++ includes

namespace libMesh
{




// ------------------------------------------------------------
// QuadratureAnisotropicNewtonCotes class definition

class QuadratureAnisotropicNewtonCotes : public QuadratureAnisotropicBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QuadratureAnisotropicNewtonCotes (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicNewtonCotes (const unsigned int _dim,
	 const std::vector<unsigned int> & vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicNewtonCotes();

  /**
   * @returns \p QuadratureType
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_NEWTON_COTES); }


 private:

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

};



// ------------------------------------------------------------

inline
QuadratureAnisotropicNewtonCotes::QuadratureAnisotropicNewtonCotes(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicBase(d,o)
{
}

inline
QuadratureAnisotropicNewtonCotes::QuadratureAnisotropicNewtonCotes(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicBase(d,o,use_anisotropic)
{
  
  
}


inline
QuadratureAnisotropicNewtonCotes::~QuadratureAnisotropicNewtonCotes()
{
}


} //



#endif // included_QuadratureAnisotropicNewtonCotes
