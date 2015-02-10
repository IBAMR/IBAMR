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



#ifndef LIBMESH_ENUM_QUADRATURE_ANISOTROPIC_TYPE_H
#define LIBMESH_ENUM_QUADRATURE_ANISOTROPIC_TYPE_H

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_quadrature_type.h"
// ------------------------------------------------------------
// enum QuadratureType definition
namespace libMeshEnums {

  /**
   * Defines an \p enum for currently available quadrature rules.
   */
  /*
  enum QuadratureType {QGAUSS            = 0,

		       QJACOBI_1_0       = 1,
		       QJACOBI_2_0       = 2,

		       QSIMPSON          = 3,
		       QTRAP             = 4,
		       QGRID             = 5,
		       QGRUNDMANN_MOLLER = 6,
		       QMONOMIAL         = 7,
		       QCONICAL          = 8,

		       QCLOUGH           = 21,

		       INVALID_Q_RULE    = 127};
		       */
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
  switch (static_cast<MoreQuadratureType>(type)) // need to modify the rool
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
#endif // LIBMESH_ENUM_QUADRATURE_ANISOTROPIC_TYPE_H
