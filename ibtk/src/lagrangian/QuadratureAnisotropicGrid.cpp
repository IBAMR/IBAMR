
#include "ibtk/QuadratureAnisotropicGrid.h"


namespace libMesh
{



void QuadratureAnisotropicGrid::init_1D(const ElemType,
                    unsigned int)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules


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

  //-----------------------------------------------------------------------
  // 2D quadrature rules

if (!QuadratureAnisotropicGrid::d_use_anisotropic)
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
	QuadratureAnisotropicGrid q1D(1,_order);
	q1D.init(EDGE2);
	tensor_product_quad( q1D );
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
	QuadratureAnisotropicGrid q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QuadratureAnisotropicGrid q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);	
	
	tensorProductForQuad( q1D1, q1D2 );
	
	
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



void QuadratureAnisotropicGrid::init_3D(const ElemType type_in,
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
	QuadratureAnisotropicGrid q1D1(1,d_vec_order[0]);
	q1D1.init(EDGE2);
	QuadratureAnisotropicGrid q1D2(1,d_vec_order[1]);
	q1D2.init(EDGE2);
	QuadratureAnisotropicGrid q1D3(1,d_vec_order[2]);
	q1D3.init(EDGE2);
	tensorProductForHex( q1D1,q1D2,q1D3 );

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
