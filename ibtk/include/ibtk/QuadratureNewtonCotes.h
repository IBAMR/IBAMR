


#ifndef LIBMESH_QUADRATURE_NEWTONCOTES_H
#define LIBMESH_QUADRATURE_NEWTONCOTES_H

// Local includes
#include "libmesh/quadrature.h"
#include "ibtk/EnumQuadratureAnisotropicType.h"
// C++ includes

namespace libMesh
{




/**
 * This class creates quadrature points based on Newton-Cotes method
 * 
 */

// ------------------------------------------------------------
// QuadratureNewtonCotes class definition

class QuadratureNewtonCotes : public QBase
{
 public:

  /**
   * Constructor.  
   */
  QuadratureNewtonCotes (const unsigned int _dim,
	 const Order _order=INVALID_ORDER);

  /**
   * Destructor.
   */
  ~QuadratureNewtonCotes();

  /**
   * @returns \p QuadratureNewtonCotes
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QNEWTON_COTES); }


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
QuadratureNewtonCotes::QuadratureNewtonCotes(const unsigned int d,
	     const Order o) : QBase(d,o)
{
}




inline
QuadratureNewtonCotes::~QuadratureNewtonCotes()
{
}


} 



#endif 
