


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
