



#ifndef included_QuadratureAnisotropicGauss
#define included_QuadratureAnisotropicGauss


#include "ibtk/QuadratureAnisotropicBase.h"
// C++ includes

namespace libMesh
{





// ------------------------------------------------------------
// QuadratureAnisotropicGauss class definition

class QuadratureAnisotropicGauss : public QuadratureAnisotropicBase

{
 public:

  /**
   * Constructor.  
   */
  QuadratureAnisotropicGauss (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicGauss (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicGauss();

  /**
   * @returns \p QuadratureType
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_GAUSS); }


 private:

  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);

};



// ------------------------------------------------------------
// QGauss class members
inline
QuadratureAnisotropicGauss::QuadratureAnisotropicGauss(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicBase(d,o)
{
}

inline
QuadratureAnisotropicGauss::QuadratureAnisotropicGauss(const unsigned int d,
	     const std::vector<unsigned>& o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicBase(d,o,use_anisotropic)
{
  
  
}


inline
QuadratureAnisotropicGauss::~QuadratureAnisotropicGauss()
{
}


} // namespace libMesh



#endif // included_QuadratureAnisotropicGauss 
