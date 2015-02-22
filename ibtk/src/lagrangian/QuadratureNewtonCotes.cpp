// Filename: QuadratureNewtonCotes.cpp
// Created on 22 Feb 2015 by Wenjun Kou
//
// Copyright (c) 2002-2015, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "ibtk/QuadratureNewtonCotes.h"

namespace libMesh
{
  
void QuadratureNewtonCotes::init_1D(const ElemType,
                    unsigned int)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules --> only support order = 1, 2 and 3.

 
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

void QuadratureNewtonCotes::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules



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
	QuadratureNewtonCotes q1D(1,_order);
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
                _points[pt](0) = 1.0 / (double) n_interval * ((double)i + 0.5); 
                _points[pt](1) = 1.0 / (double) n_interval * ((double)i + 0.5); 
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


void QuadratureNewtonCotes::init_3D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules

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
	QuadratureNewtonCotes q1D(1,_order);
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
                    _points[pt](0) = 1.0 / (double) n_interval * ((double)i + 0.5); 
                    _points[pt](1) = 1.0 / (double) n_interval * ((double)i + 0.5); 
                    _points[pt](2) = 1.0 / (double) n_interval * ((double)i + 0.5); 
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

	QuadratureNewtonCotes q1D(1,_order);
	QuadratureNewtonCotes q2D(2,_order);

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
      {
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
