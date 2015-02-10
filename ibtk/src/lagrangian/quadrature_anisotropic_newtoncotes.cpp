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



// C++ includes

// Local includes
#include "ibtk/quadrature_anisotropic_newtoncotes.h"


namespace libMesh
{



void QAnisotropicNewtonCotes::init_1D(const ElemType,
                    unsigned int)
{
  
  //----------------------------------------------------------------------
  // 1D quadrature rules --> only support order = 1,2,3 and 4.

  // We ignore p - the grid rule is just for experimentation
  // 0.7071067811865475 (2.0/3.0) for Degree 4
  _points.resize(_order + 1);
  _weights.resize(_order + 1); // npts = order +1;
  
  double L=2.0; // [-1, 1]
  // first to do the points
  for (int i = 0; i != _order + 1; ++i)
    {
      _points[i](0) = 2.0 * (double)(i+1) / (double)(_order +1) - 1.0;
      
    }
  // then to do the weight
  if (_order == 1)
  { _weights[0]=L/2.0; _weights[1]=L/2.0;return;}
  if (_order == 2) // default weight can be negative, --> so we change the points
  { _weights[0]=L/3.0;_weights[1]=L/2.0;_weights[2]=L/2.0; _points[0](0)=-0.7071067811865475; _points[2](0)=0.7071067811865475;return;}
  if (_order == 3)
  {_weights[0] = L*11.0/24.0; _weights[1] = L*1/24.0;_weights[2] = L*1/24.0;_weights[3] = L*11.0/24.0;return;}
  if (_order >= 4) // there is negative coefficient --> not recommended
  { 
    libMesh::err << "ERROR: order > 3 is not supported, such as here: " << _order<< std::endl;
    libmesh_error();
    return;
  }
  // _order <= 0
   libMesh::err << "ERROR: _order should be at least 1: not " << _order<< std::endl;
   libmesh_error();
  return;
}


void QAnisotropicNewtonCotes::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules

  // We ignore p - the grid rule is just for experimentation
if (!QAnisotropicNewtonCotes::d_use_anisotropic) // this qp rule is not suggested, if anisotropic rule is not used
{
  switch (type_in)
    {


      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	// We compute the 2D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QAnisotropicNewtonCotes q1D(1,_order);
	q1D.init(EDGE2);
	tensor_product_quad( q1D );
	return;
      }


      //---------------------------------------------
      // Triangle quadrature rules
    case TRI3:
    case TRI6:
      {
        _points.resize((_order + 1)*(_order + 2)/2);
        _weights.resize((_order + 1)*(_order + 2)/2);

        unsigned int pt = 0;
        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                _points[pt](0) = (double)i / (double)_order;
                _points[pt](1) = (double)j / (double)_order;
                _weights[pt] = 1.0 / (double)(_order+1) /
                  (double)(_order+2);
                pt++;
              }
          }
	return;
      }

      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "Element type not supported!:" << type_in << std::endl;
	libmesh_error();
      }
    }
}
else // use the order
{
  switch (type_in)
    {


      //---------------------------------------------
      // Quadrilateral quadrature rules
    case QUAD4:
    case QUAD8:
    case QUAD9:
      {
	// We compute the 2D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QAnisotropicNewtonCotes q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QAnisotropicNewtonCotes q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);	
	
	tensor_2product_quad( q1D1, q1D2 );
	
	
	return;
      }
      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "Only suppoer Quad_element: Element type not supported!:" << type_in << std::endl;
	libmesh_error();
      }
    }  
      
}
  libmesh_error();

  return;

#endif
}



void QAnisotropicNewtonCotes::init_3D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules

  // We ignore p - the grid rule is just for experimentation
  switch (type_in)
    {
      //---------------------------------------------
      // Hex quadrature rules
    case HEX8:
    case HEX20:
    case HEX27:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule.
	QAnisotropicNewtonCotes q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QAnisotropicNewtonCotes q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);
	QAnisotropicNewtonCotes q1D3(1,d_vec_order[2]);
	q1D3.init(EDGE2);
	tensor_3product_hex( q1D1,q1D2,q1D3 );

	return;
      }



      

      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "ERROR:Only suppoer Hex_element: Unsupported type: " << type_in << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

} // namespace libMesh
