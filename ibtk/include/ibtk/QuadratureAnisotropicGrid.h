


#ifndef included_QuadratureAnisotropicGrid
#define included_QuadratureAnisotropicGrid


#include "ibtk/QuadratureAnisotropicBase.h"
// C++ includes

namespace libMesh
{



// ------------------------------------------------------------
// QuadratureAnisotropicGrid class definition

class QuadratureAnisotropicGrid : public QuadratureAnisotropicBase
{
 public:

  /**
   * Constructor.
   */
  QuadratureAnisotropicGrid (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicGrid (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicGrid();

  /**
   * @returns \p QuadratureAnisotropicGRID
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_GRID); }


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
QuadratureAnisotropicGrid::QuadratureAnisotropicGrid(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicBase(d,o)
{
}

inline
QuadratureAnisotropicGrid::QuadratureAnisotropicGrid(const unsigned int d,
	     const std::vector<unsigned int>& o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicBase(d,o,use_anisotropic)
{
  
  
}


inline
QuadratureAnisotropicGrid::~QuadratureAnisotropicGrid()
{
}


} // 



#endif // 
