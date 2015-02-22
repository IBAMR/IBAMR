
#ifndef included_QuadratureAnisotropicCompositeGauss
#define included_QuadratureAnisotropicCompositeGauss

// Local includes

#include "ibtk/QuadratureAnisotropicCompositeBase.h"
// C++ includes

namespace libMesh
{



// ------------------------------------------------------------
// QuadratureAnisotropicCompositeGauss class definition

class QuadratureAnisotropicCompositeGauss : public QuadratureAnisotropicCompositeBase
{
 public:

  /**
   * Constructor.
   */
  QuadratureAnisotropicCompositeGauss (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicCompositeGauss (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_order, bool use_anisotropic =true);
  

   
  QuadratureAnisotropicCompositeGauss (const unsigned int d,
	 const std::vector<unsigned int>& vec_interval_num, const unsigned int nqp_per_interval);
  
  QuadratureAnisotropicCompositeGauss (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval); // only for 1D
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicCompositeGauss();

  /**
   * @returns \p QuadratureAnisotropicGRID
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_COMPOSITE_GAUSS); }


 private:
   
  void buildStandard1D( std::vector<Real>& standard_weights,  std::vector<Real>& standard_cods);
    // virtual function from Qbase class
  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0);
  

};



// ------------------------------------------------------------

inline
QuadratureAnisotropicCompositeGauss::QuadratureAnisotropicCompositeGauss(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicCompositeBase(d,o)
{
}

inline
QuadratureAnisotropicCompositeGauss::QuadratureAnisotropicCompositeGauss(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicCompositeBase(d,o,use_anisotropic)
{
 
}

inline
QuadratureAnisotropicCompositeGauss::QuadratureAnisotropicCompositeGauss (const unsigned int d,
	 const std::vector<unsigned int> &  vec_interval_num, const unsigned int nqp_per_interval)
	 :QuadratureAnisotropicCompositeBase(d, vec_interval_num, nqp_per_interval)
{

}

  inline
QuadratureAnisotropicCompositeGauss::QuadratureAnisotropicCompositeGauss (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval) // only for 1D
	:QuadratureAnisotropicCompositeBase(d,interval_num,nqp_per_interval)
{
 
  
}

inline
QuadratureAnisotropicCompositeGauss::~QuadratureAnisotropicCompositeGauss()
{
}


} // 



#endif //  included_QuadratureAnisotropicCompositeGauss
