// Filename: EnumQuadratureAnisotropicType.h
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


#ifndef included_EnumQuadratureAnisotropicType
#define included_EnumQuadratureAnisotropicType

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_quadrature_type.h"
// ------------------------------------------------------------
// enum QuadratureType definition
namespace libMeshEnums {


  enum MoreQuadratureType{    QNEWTON_COTES = 10, //  OPEN TYPE WITH NO ANISOTROPIC PROPERTY
			      
			      // Anisotropic type
			      QANISOTROPIC_NEWTON_COTES =30,
			      QANISOTROPIC_GRID =31, // not recommeded
			      QANISOTROPIC_GAUSS =32,
			      
			      // Anisotropic_Composite_type
			      QANISOTROPIC_COMPOSITE_NEWTON_COTES = 40,
			      QANISOTROPIC_COMPOSITE_GAUSS = 41,
			      QANISOTROPIC_COMPOSITE_GRID = 42 //NOT RECOMMENDED
  };
			      
}

using namespace libMeshEnums;


// >> for preprocess in IBFEMethod, we want consistency: return QuadratureType
inline
QuadratureType
string_to_quadrature_type(const std::string& s)
{ 
  if (s.compare("QNEWTON_COTES") == 0)
    return static_cast<QuadratureType>(QNEWTON_COTES);
  
  if (s.compare("QANISOTROPIC_NEWTON_COTES") == 0)
    return static_cast<QuadratureType>(QANISOTROPIC_NEWTON_COTES);  
  if (s.compare("QANISOTROPIC_GRID") == 0)
    return static_cast<QuadratureType>(QANISOTROPIC_GRID);
  if (s.compare("QANISOTROPIC_GAUSS") == 0)
    return static_cast<QuadratureType>(QANISOTROPIC_GAUSS);
  
  if (s.compare("QANISOTROPIC_COMPOSITE_NEWTON_COTES") == 0)
    return static_cast<QuadratureType>(QANISOTROPIC_COMPOSITE_NEWTON_COTES);
  if (s.compare("QANISOTROPIC_COMPOSITE_GAUSS") == 0)
    return static_cast<QuadratureType>(QANISOTROPIC_COMPOSITE_GAUSS);
  
  // else
   return libMesh::Utility::string_to_enum<QuadratureType> (s); 
}

inline
std::string
quadrature_type_to_string(QuadratureType type)
{
  switch (static_cast<MoreQuadratureType>(type)) 
      {
	case QNEWTON_COTES:
	  return "QNEWTON_COTES";
	case QANISOTROPIC_GAUSS:
	  return "QANISOTROPIC_GAUSS";
	case QANISOTROPIC_GRID:
	  return "QANISOTROPIC_GRID";
	case QANISOTROPIC_NEWTON_COTES:
	  return "QANISOTROPIC_NEWTON_COTES";
	case QANISOTROPIC_COMPOSITE_GAUSS:
	  return "QANISOTROPIC_COMPOSITE_GAUSS";
	case QANISOTROPIC_COMPOSITE_NEWTON_COTES:
	  return "QANISOTROPIC_COMPOSITE_NEWTON_COTE";
	
	default:
          return  libMesh::Utility::enum_to_string<QuadratureType> (type);      
      }
  
}
#endif // included_EnumQuadratureAnisotropicType
