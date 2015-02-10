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



#ifndef LIBMESH_QAnisotropicBase_H
#define LIBMESH_QAnisotropicBase_H

// Local includes
#include "libmesh/quadrature.h"
#include "ibtk/enum_quadrature_anisotropic_type.h"
// C++ includes

namespace libMesh
{




/**
 * This class creates anisotropic quadrature points on a non-uniform grid, Order
 * points on a side.
 */

// ------------------------------------------------------------
// QAnisotropicBase class definition

class QAnisotropicBase : public QBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QAnisotropicBase (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QAnisotropicBase (const unsigned int _dim,
	 const std::vector<unsigned int> & vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QAnisotropicBase();

  /**
   * @returns \p QAnisotropicBase
   */
  //QuadratureType type() const { return QAnisotropicBase; }

  void tensor_2product_quad (const QAnisotropicBase& q1D1, const QAnisotropicBase& q1D2);

  /**
   * Computes the tensor product quadrature rule
   * [q1D x q1D x q1D] from the 1D rule q1D.
   * Used in the init_3D routines for
   * hexahedral element types.
   */
  void tensor_3product_hex (const QAnisotropicBase& q1D1, const QAnisotropicBase& q1D2, const QAnisotropicBase& q1D3);
  
  //virtual // to check_vec_nps; will be overriden by composite rule
  std::vector<unsigned int> check_vec_nps( void) {return d_vec_order;}
protected:
  std::vector<unsigned int> d_vec_order;
  bool d_use_anisotropic;
  
  
  // virtual function from QBase class
  void init_1D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0){} // intentionally do nothing
  void init_2D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0){}
  void init_3D (const ElemType _type=INVALID_ELEM,
		unsigned int p_level=0){}

  
};



// ------------------------------------------------------------
// QGauss class members
inline
QAnisotropicBase::QAnisotropicBase(const unsigned int d,
	     const unsigned int o) 
	  : QBase(d,static_cast<Order>(o)),d_vec_order(d,o),d_use_anisotropic(false)
{
  if (d==1) d_use_anisotropic = true; // this is called when init_2D, if d>1, only one int means uniform qp rule
}

inline
QAnisotropicBase::QAnisotropicBase(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QBase(d, static_cast<Order>(o[0])), d_vec_order(o), d_use_anisotropic(use_anisotropic)
{
  
  
}


//q1d1 x q1d2
inline 
void QAnisotropicBase::tensor_2product_quad(const QAnisotropicBase& q1D1, const QAnisotropicBase& q1D2)
{

  const unsigned int np1 = q1D1.n_points();
  
  const unsigned int np2 = q1D2.n_points();
 
  _points.resize(np1 * np2);

  _weights.resize(np1 * np2);
  
  unsigned int q=0;

  for (unsigned int j=0; j<np2; j++)
    for (unsigned int i=0; i<np1; i++)
      {
	_points[q](0) = q1D1.qp(i)(0);
	_points[q](1) = q1D2.qp(j)(0);

	_weights[q] = q1D1.w(i)*q1D2.w(j);

	q++;
      }
}




inline
void QAnisotropicBase::tensor_3product_hex(const QAnisotropicBase& q1D1, const QAnisotropicBase& q1D2,const QAnisotropicBase& q1D3 )
{
  const unsigned int np3 = q1D3.n_points();
  const unsigned int np1 = q1D1.n_points();
  const unsigned int np2 = q1D2.n_points();
  _points.resize(np1 * np2 * np3);

  _weights.resize(np1 * np2 * np3);

  unsigned int q=0;

  for (unsigned int k=0; k<np3; k++)
    for (unsigned int j=0; j<np2; j++)
      for (unsigned int i=0; i<np1; i++)
	{
	  _points[q](0) = q1D1.qp(i)(0);
	  _points[q](1) = q1D2.qp(j)(0);
	  _points[q](2) = q1D3.qp(k)(0);

	  _weights[q] = q1D1.w(i) * q1D2.w(j) * q1D3.w(k);

	  q++;
	}
}





inline
QAnisotropicBase::~QAnisotropicBase()
{
}


} // namespace libMesh



#endif // LIBMESH_ADAPTIVE_QUADRATURE_GRID_H
