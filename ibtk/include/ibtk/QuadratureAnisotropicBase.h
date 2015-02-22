



#ifndef included_QuadratureAnisotropicBase
#define included_QuadratureAnisotropicBase

// Local includes
#include "libmesh/quadrature.h"
#include "ibtk/EnumQuadratureAnisotropicType.h"
// C++ includes

namespace libMesh
{




/**
 * This class is the base class for anisotropic quadrature rule
 *
 */

// ------------------------------------------------------------
// QuadratureAnisotropicBase class definition

class QuadratureAnisotropicBase : public QBase
{
 public:

  /**
   * Constructor.  
   */
  QuadratureAnisotropicBase (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicBase (const unsigned int _dim,
	 const std::vector<unsigned int> & vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicBase();

  /**
   * @returns \p QuadratureAnisotropicBase
   */


  void tensorProductForQuad (const QuadratureAnisotropicBase& q1D1, const QuadratureAnisotropicBase& q1D2);

  /**
   * Computes the tensor product quadrature rule
   * from the 1D rule q1D.
   * This is used in the init_3D routines for
   * hexahedral element types.
   */
  void tensorProductForHex (const QuadratureAnisotropicBase& q1D1, const QuadratureAnisotropicBase& q1D2, const QuadratureAnisotropicBase& q1D3);
  
  //virtual // to check_vec_nps; will be overriden by composite rule
  std::vector<unsigned int> checkVecNumQPs( void) {return d_vec_order;}
  
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

inline
QuadratureAnisotropicBase::QuadratureAnisotropicBase(const unsigned int d,
	     const unsigned int o) 
	  : QBase(d,static_cast<Order>(o)),d_vec_order(d,o),d_use_anisotropic(false)
{
  if (d==1) d_use_anisotropic = true; // this is called when init_2D, if d>1, only one int means uniform qp rule
}

inline
QuadratureAnisotropicBase::QuadratureAnisotropicBase(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QBase(d, static_cast<Order>(o[0])), d_vec_order(o), d_use_anisotropic(use_anisotropic)
{
  
  
}


//q1d1 x q1d2
inline 
void QuadratureAnisotropicBase::tensorProductForQuad(const QuadratureAnisotropicBase& q1D1, const QuadratureAnisotropicBase& q1D2)
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
void QuadratureAnisotropicBase::tensorProductForHex(
    const QuadratureAnisotropicBase& q1D1, const QuadratureAnisotropicBase& q1D2,const QuadratureAnisotropicBase& q1D3 )
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
QuadratureAnisotropicBase::~QuadratureAnisotropicBase()
{
}


} // namespace libMesh



#endif // included_QuadratureAnisotropicBase
