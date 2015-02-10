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



// Local includes
#include "ibtk/quadrature_newtoncotes.h"

namespace libMesh
{
  
void QNewtonCotes::init_1D(const ElemType,
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
  { _weights[0]=L/3.0; _weights[1]=L/2.0; _weights[2]=L/2.0; 
    _points[0](0) =-0.7071067811865475; _points[2](0)=0.7071067811865475;
    return;}
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

void QNewtonCotes::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules

  // We ignore p - the grid rule is just for experimentation

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
	QNewtonCotes q1D(1,_order);
	q1D.init(EDGE2);
	tensor_product_quad( q1D );
	return;
      }


      //---------------------------------------------
      // Triangle quadrature rules
      // >> equal weights: open netwon-cotes; 
      // reference region: pt(0,0) -- pt(1,0) -- pt(0,1)
    case TRI3:
    case TRI6:
      {
        _points.resize((_order + 1)*(_order + 2)/2);
        _weights.resize((_order + 1)*(_order + 2)/2);
	unsigned n_interval = _order +1;
        unsigned int pt = 0;
        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                _points[pt](0) = 1.0 / (double) n_interval * ((double)i + 0.5); //(double)i / (double)_order;
                _points[pt](1) = 1.0 / (double) n_interval * ((double)i + 0.5); //(double)j / (double)_order;
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

  libmesh_error();

  return;

#endif
}


void QNewtonCotes::init_3D(const ElemType type_in,
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
	QNewtonCotes q1D(1,_order);
	q1D.init(EDGE2);

	tensor_product_hex( q1D );

	return;
      }



      //---------------------------------------------
      // Tetrahedral quadrature rules
    case TET4:
    case TET10:      
      // >> equal weights: open netwon-cotes; 
      // reference region: pt(0,0,0) -- pt(1,0,0) -- pt(0,0,1) -- pt(0,1,0)
      {
        _points.resize((_order+1)*(_order+2)*(_order+3)/6);
        _weights.resize((_order+1)*(_order+2)*(_order+3)/6);
	unsigned n_interval = _order +1;
        unsigned int pt = 0;
        for (int i = 0; i != _order + 1; ++i)
          {
            for (int j = 0; j != _order + 1 - i; ++j)
              {
                for (int k = 0; k != _order + 1 - i - j; ++k)
                  {
                    _points[pt](0) = 1.0 / (double) n_interval * ((double)i + 0.5); //(double)i / (double)_order;
                    _points[pt](1) = 1.0 / (double) n_interval * ((double)i + 0.5); // (double)j / (double)_order;
                    _points[pt](2) = 1.0 / (double) n_interval * ((double)i + 0.5); //(double)k / (double)_order;
                    _weights[pt] = 1.0 / (double)(_order+1) /
                      (double)(_order+2) / (double)(_order+3);
                    pt++;
                  }
              }
          }
	return;
      }


      // Prism quadrature rules
    case PRISM6:
    case PRISM15:
    case PRISM18:
      {
	// We compute the 3D quadrature rule as a tensor
	// product of the 1D quadrature rule and a 2D
	// triangle quadrature rule

	QNewtonCotes q1D(1,_order);
	QNewtonCotes q2D(2,_order);

	// Initialize
	q1D.init(EDGE2);
	q2D.init(TRI3);

	tensor_product_prism(q1D, q2D);

	return;
      }



      //---------------------------------------------
      // Pyramid
    case PYRAMID5:
    case PYRAMID14:
      { /*
        _points.resize((_order+1)*(_order+2)*(_order+3)/6);
        _weights.resize((_order+1)*(_order+2)*(_order+3)/6);
	
        unsigned int pt = 0;
        for (int k = 0; k != _order + 1; ++k)
          {
            for (int i = 0; i != _order + 1 - k; ++i)
              {
                for (int j = 0; j != _order + 1 - k; ++j)
                  {
                    _points[pt](0) = 2.0 * (double)i / (double)_order
                      - 1.0 + (double)k / (double)_order;
                    _points[pt](1) = 2.0 * (double)j / (double)_order
                      - 1.0 + (double)k / (double)_order;
                    _points[pt](2) = (double)k / (double)_order;
                    _weights[pt] = 1.0 / (double)(_order+1) /
                      (double)(_order+2) / (double)(_order+3);
                    pt++;
                  }
              }
          }
	  */	
	libMesh::err << "Element type not supported!:" << type_in << std::endl;
	libmesh_error();
	return;
      }



      //---------------------------------------------
      // Unsupported type
    default:
      {
	libMesh::err << "ERROR: Unsupported type: " << type_in << std::endl;
	libmesh_error();
      }
    }

  libmesh_error();

  return;

#endif
}

} // namespace libMesh
