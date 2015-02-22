// Filename: QuadratureAnisotropicNewtonCotes.h
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


#ifndef included_QuadratureAnisotropicNewtonCotes
#define included_QuadratureAnisotropicNewtonCotes


#include "ibtk/QuadratureAnisotropicBase.h"
// C++ includes

namespace libMesh
{




// ------------------------------------------------------------
// QuadratureAnisotropicNewtonCotes class definition

class QuadratureAnisotropicNewtonCotes : public QuadratureAnisotropicBase
{
 public:

  /**
   * Constructor.  Declares the order of the quadrature rule.
   */
  QuadratureAnisotropicNewtonCotes (const unsigned int _dim,
	 const unsigned int _order=INVALID_ORDER);
  
  QuadratureAnisotropicNewtonCotes (const unsigned int _dim,
	 const std::vector<unsigned int> & vec_order,bool use_anisotropic =true);
  /**
   * Destructor.
   */
  ~QuadratureAnisotropicNewtonCotes();

  /**
   * @returns \p QuadratureType
   */
  QuadratureType type() const { return static_cast<QuadratureType>(QANISOTROPIC_NEWTON_COTES); }


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
QuadratureAnisotropicNewtonCotes::QuadratureAnisotropicNewtonCotes(const unsigned int d,
	     const unsigned int o) : QuadratureAnisotropicBase(d,o)
{
}

inline
QuadratureAnisotropicNewtonCotes::QuadratureAnisotropicNewtonCotes(const unsigned int d,
	     const std::vector<unsigned int> & o,
	     bool use_anisotropic) 
	: QuadratureAnisotropicBase(d,o,use_anisotropic)
{
  
  
}


inline
QuadratureAnisotropicNewtonCotes::~QuadratureAnisotropicNewtonCotes()
{
}


} //



#endif // included_QuadratureAnisotropicNewtonCotes
