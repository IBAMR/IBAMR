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



#ifndef LIBMESH_ANISOTROPIC_COMPOSITE_QUADRATURE_GAUSS_H
#define LIBMESH_ANISOTROPIC_COMPOSITE_QUADRATURE_GAUSS_H

// Local includes
//#include "libmesh/quadrature.h"
#include "ibtk/quadrature_anisotropic_composite.h"
// C++ includes

namespace libMesh
{




/**
 * This class creates quadrature points on a uniform grid, Order
 * points on a side.
 */

// ------------------------------------------------------------
// QAnisotropicCompositeGauss class definition

class QAnisotropicCompositeGauss : public QAnisotropicCompositeBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QAnisotropicCompositeGauss (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QAnisotropicCompositeGauss (const unsigned int _dim,
	 const std::vector<unsigned int>& vec_order, bool use_anisotropic =true);
  

   
  QAnisotropicCompositeGauss (const unsigned int d,
	 const std::vector<unsigned int>& vec_interval_num, const unsigned int nqp_per_interval);
  
  QAnisotropicCompositeGauss (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval); // only for 1D
  /**
   * Destructor.
   */
  ~QAnisotropicCompositeGauss();

  /**
   * @returns \p QAnisotropicGRID
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_COMPOSITE_GAUSS); }


 private:
   
  void build_standard_1D( std::vector<Real>& standard_weights,  std::vector<Real>& standard_cods);
    
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
QAnisotropicCompositeGauss::QAnisotropicCompositeGauss(const unsigned int d,
	     const unsigned int o) : QAnisotropicCompositeBase(d,o)
{
}

inline
QAnisotropicCompositeGauss::QAnisotropicCompositeGauss(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QAnisotropicCompositeBase(d,o,use_anisotropic)
{
 
}

inline
QAnisotropicCompositeGauss::QAnisotropicCompositeGauss (const unsigned int d,
	 const std::vector<unsigned int> &  vec_interval_num, const unsigned int nqp_per_interval)
	 :QAnisotropicCompositeBase(d, vec_interval_num, nqp_per_interval)
{

}

  inline
QAnisotropicCompositeGauss::QAnisotropicCompositeGauss (const unsigned int d,
	 const unsigned int interval_num, const unsigned int nqp_per_interval) // only for 1D
	:QAnisotropicCompositeBase(d,interval_num,nqp_per_interval)
{
 
  
}

inline
QAnisotropicCompositeGauss::~QAnisotropicCompositeGauss()
{
}


} // namespace libMesh



#endif //  LIBMESH_ANISOTROPIC_QUADRATURE_GRID_H
