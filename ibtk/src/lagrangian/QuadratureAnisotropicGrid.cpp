// Filename: QuadratureAnisotropicGrid.cpp
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
#include "ibtk/QuadratureAnisotropicGrid.h"


namespace libMesh
{



void QuadratureAnisotropicGrid::init_1D(const ElemType,
                    unsigned int)
{



  _points.resize(_order + 1);
  _weights.resize(_order + 1);
  for (int i = 0; i != _order + 1; ++i)
    {
      _points[i](0) = 2.0 * (double)i / (double)_order - 1.0;
      _weights[i] = 2.0 / (double)(_order + 1);
    }
  return;
}


void QuadratureAnisotropicGrid::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1



if (!QuadratureAnisotropicGrid::d_use_anisotropic)
{
  switch (type_in)
    {


    case QUAD4:
    case QUAD8:
    case QUAD9:
      {

	QuadratureAnisotropicGrid q1D(1,_order);
	q1D.init(EDGE2);
	tensor_product_quad( q1D );
	return;
      }


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


    case QUAD4:
    case QUAD8:
    case QUAD9:
      {

	QuadratureAnisotropicGrid q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QuadratureAnisotropicGrid q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);	
	
	tensorProductForQuad( q1D1, q1D2 );
	
	
	return;
      }

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



void QuadratureAnisotropicGrid::init_3D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM == 3


  switch (type_in)
    {

    case HEX8:
    case HEX20:
    case HEX27:
      {

	QuadratureAnisotropicGrid q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QuadratureAnisotropicGrid q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);
	QuadratureAnisotropicGrid q1D3(1,d_vec_order[2]);
	q1D3.init(EDGE2);
	tensorProductForHex( q1D1,q1D2,q1D3 );

	return;
      }


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
