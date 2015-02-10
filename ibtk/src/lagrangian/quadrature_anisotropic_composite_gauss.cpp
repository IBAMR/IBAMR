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
#include "ibtk/quadrature_anisotropic_composite_gauss.h"


namespace libMesh
{

void QAnisotropicCompositeGauss::build_standard_1D( std::vector<Real> & standard_weights,  std::vector<Real> & standard_cods)
{
 
   switch(d_num_qps) // number of means points
  {
  
    case 1:
      {
	
	standard_cods.resize (1);
	standard_weights.resize(1);

	standard_cods[0]  = 0.;

	standard_weights[0]    = 2.;

	return;
      }
    case 2:
    //case 3:
      {
	standard_cods.resize (2);
	standard_weights.resize(2);

	standard_cods[0] = -5.7735026918962576450914878050196e-01L; // -sqrt(3)/3
	standard_cods[1]    = -standard_cods[0];

	standard_weights[0]   = 1.;
	standard_weights[1]   = standard_weights[0];

	return;
      }
    case 3:
    //case 5:
      {
	standard_cods.resize (3);
	standard_weights.resize(3);

	standard_cods[ 0]= -7.7459666924148337703585307995648e-01L;
	standard_cods[ 1]= 0.;
	standard_cods[ 2]    = -standard_cods[0];

	standard_weights[ 0]   = 5.5555555555555555555555555555556e-01L;
	standard_weights[ 1]   = 8.8888888888888888888888888888889e-01L;
	standard_weights[ 2]   = standard_weights[0];

	return;
      }
    case 4:
    //case 7:
      {
	standard_cods.resize (4);
	standard_weights.resize(4);

	standard_cods[ 0]= -8.6113631159405257522394648889281e-01L;
	standard_cods[ 1]= -3.3998104358485626480266575910324e-01L;
	standard_cods[ 2]    = -standard_cods[1];
	standard_cods[ 3]    = -standard_cods[0];

	standard_weights[ 0]   = 3.4785484513745385737306394922200e-01L;
	standard_weights[ 1]   = 6.5214515486254614262693605077800e-01L;
	standard_weights[ 2]   = standard_weights[1];
	standard_weights[ 3]   = standard_weights[0];

	return;
      }
    case 5:
    //case 9:
      {
	standard_cods.resize (5);
	standard_weights.resize(5);

	standard_cods[ 0]= -9.0617984593866399279762687829939e-01L;
	standard_cods[ 1]= -5.3846931010568309103631442070021e-01L;
	standard_cods[ 2]= 0.;
	standard_cods[ 3]    = -standard_cods[1];
	standard_cods[ 4]    = -standard_cods[0];

	standard_weights[ 0]   = 2.3692688505618908751426404071992e-01L;
	standard_weights[ 1]   = 4.7862867049936646804129151483564e-01L;
	standard_weights[ 2]   = 5.6888888888888888888888888888889e-01L;
	standard_weights[ 3]   = standard_weights[1];
	standard_weights[ 4]   = standard_weights[0];

	return;
      }
    default:
      {
	libMesh::err << "ERROR:only support n_pts in (1,5), not this one: " << d_num_qps << std::endl;
	libmesh_error();
      
      }
      
  } // switch

  libmesh_error();  
  return;
}


void QAnisotropicCompositeGauss::init_1D(const ElemType,
                    unsigned int)
{
  //----------------------------------------------------------------------
  // 1D quadrature rules
  // here _order means number of points
 if ((!d_use_composite) || (!d_use_anisotropic) || (d_num_qps <1))
 {	
   libMesh::err << "ERROR:input is invalid for one anisotropic composite qp rule " << std::endl;
   libmesh_error();
   return ;
 }
 
 // step 1: build_standard_1D
 std::vector<Real> standard_weights;
 std::vector<Real> standard_cods;
 
 build_standard_1D( standard_weights,  standard_cods);

  unsigned int cur_num_intervals = _order;
  unsigned int total_qps = cur_num_intervals * d_num_qps;
  double h = 2.0 / cur_num_intervals; // sub-interval size

  _weights.resize(total_qps);
    
  _points.resize(total_qps);

  unsigned int id_interval;
  for (id_interval =0; id_interval < cur_num_intervals; id_interval++)
  {
    unsigned int start_index = id_interval * d_num_qps;
    double x1 = -1.0 + h *id_interval;
    double x2 = x1 + h;
    
    tranform_1D( _weights, _points, start_index, 
		   standard_weights, standard_cods, x1, x2);
  }
 
  return; 
 
}


void QAnisotropicCompositeGauss::init_2D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM > 1

  //-----------------------------------------------------------------------
  // 2D quadrature rules

  // We ignore p - the grid rule is just for experimentation
  if ((!d_use_composite) || (!d_use_anisotropic) || (d_num_qps <1))//(!d_use_anisotropic)
    {
  

      //---------------------------------------------
      // Unsupported type
   
	libMesh::err << "Only support anisotropic composite qp rule!:" << type_in << std::endl;
	libmesh_error();
	return;
    }

//libmesh_error();
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

	QAnisotropicCompositeGauss q1D1(1,d_vec_order[0], d_num_qps); // number of intervals
	q1D1.init(EDGE2);
	
	QAnisotropicCompositeGauss q1D2(1,d_vec_order[1], d_num_qps); //d_vec_order[1]);
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
 libmesh_error();     
 return;

#endif
}



void QAnisotropicCompositeGauss::init_3D(const ElemType type_in,
                    unsigned int)
{
#if LIBMESH_DIM == 3

  //-----------------------------------------------------------------------
  // 3D quadrature rules
    if ((!d_use_composite) || (!d_use_anisotropic) || (d_num_qps <1))//(!QAnisotropicCompositeGauss::d_use_anisotropic)
    {
  

      //---------------------------------------------
      // Unsupported type
   
	libMesh::err << "Only support anisotropic composite qp rule!:" << type_in << std::endl;
	libmesh_error();
	return;
    }
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

	QAnisotropicCompositeGauss q1D1(1,d_vec_order[0], d_num_qps); // number of intervals
	q1D1.init(EDGE2);
	
	QAnisotropicCompositeGauss q1D2(1,d_vec_order[1], d_num_qps); //d_vec_order[1]);
	q1D2.init(EDGE2);
	
	QAnisotropicCompositeGauss q1D3(1,d_vec_order[2], d_num_qps); //d_vec_order[2]);
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
