// Filename: FESystemDataCache.h
// Created on 10 Mar 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_FESystemDataCache
#define included_FESystemDataCache

/////////////////////////////// INCLUDES /////////////////////////////////////

// BLITZ++ INCLUDES
#include <blitz/array.h>

// LIBMESH INCLUDES
#define LIBMESH_REQUIRE_SEPARATE_NAMESPACE
#include <../base/variable.h>
#include <point.h>
#include <quadrature.h>
#include <system.h>
#include <vector_value.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FESystemDataCache provides a mechanism for computing and storing
 * FE system-, mesh-, and quadrature-rule dependent data (e.g., the values of
 * shape functions and their derivatives) for reuse.
 */
class FESystemDataCache
{
public:
    /*!
     * \brief Constructor.
     */
    FESystemDataCache(
        const libMesh::System& system,
        const libMesh::MeshBase& mesh);

    /*!
     * \brief Destructor.
     */
    ~FESystemDataCache();

    /*!
     * \brief Compute the cached data items over the specified range of
     * elements.
     */
    template<typename ElementIterator>
    void
    computeCachedData(
        const ElementIterator& elem_begin,
        const ElementIterator& elem_end,
        libMesh::QBase* qrule=NULL,
        libMesh::QBase* qrule_face=NULL);

    /*!
     * \brief Clear the cached data items.
     */
    void
    clearCachedData();

    /*!
     * \brief Configuration options.
     */
    inline bool& compute_q_point()  { return d_compute_q_point;  }
    inline bool& comute_JxW()       { return d_compute_JxW;      }
    inline bool& compute_phi()      { return d_compute_phi;      }
    inline bool& compute_phi_JxW()  { return d_compute_phi_JxW;  }
    inline bool& compute_dphi()     { return d_compute_dphi;     }
    inline bool& compute_dphi_JxW() { return d_compute_dphi_JxW; }

    inline bool& compute_q_point_face()    { return d_compute_q_point_face;    }
    inline bool& comute_JxW_face()         { return d_compute_JxW_face;        }
    inline bool& compute_normal_face()     { return d_compute_normal_face;     }
    inline bool& compute_normal_JxW_face() { return d_compute_normal_JxW_face; }
    inline bool& compute_phi_face()        { return d_compute_phi_face;        }
    inline bool& compute_phi_JxW_face()    { return d_compute_phi_JxW_face;    }
    inline bool& compute_dphi_face()       { return d_compute_dphi_face;       }
    inline bool& compute_dphi_JxW_face()   { return d_compute_dphi_JxW_face;   }

    /*!
     * \brief Cached data items.
     */
    inline int num_elems() const { return d_elems.size(); }

    inline const blitz::Array<libMesh::Elem*,1>& elems() const { return d_elems; }
    inline const blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1>& dof_indices() const { return d_dof_indices; }

    inline const blitz::Array<blitz::Array<libMesh::Point,1>,1>& q_point() const { return d_q_point; }
    inline const blitz::Array<blitz::Array<double,1>,1>& JxW() const { return d_JxW; }
    inline const blitz::Array<blitz::Array<double,2>,1>& phi() const { return d_phi; }
    inline const blitz::Array<blitz::Array<double,2>,1>& phi_JxW() const { return d_phi_JxW; }
    inline const blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>& dphi() const { return d_dphi; }
    inline const blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>& dphi_JxW() const { return d_dphi_JxW; }

    inline const blitz::Array<bool,1>& elem_touches_physical_bdry() const { return d_elem_touches_physical_bdry; }
    inline const blitz::Array<blitz::Array<bool,1>,1>& elem_side_at_physical_bdry() const { return d_elem_side_at_physical_bdry; }
    inline const blitz::Array<blitz::Array<bool,1>,1>& elem_side_at_dirichlet_bdry() const { return d_elem_side_at_dirichlet_bdry; }

    inline const blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1>& q_point_face() const { return d_q_point_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<double,1>,1>,1>& JxW_face() const { return d_JxW_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1>& normal_face() const { return d_normal_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1>& normal_JxW_face() const { return d_normal_JxW_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>& phi_face() const { return d_phi_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1>& phi_JxW_face() const { return d_phi_JxW_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1>& dphi_face() const { return d_dphi_face; }
    inline const blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1>& dphi_JxW_face() const { return d_dphi_JxW_face; }

    inline libMesh::Elem* const elems(const int e) const { return d_elems(e); }
    inline const blitz::Array<std::vector<unsigned int>,1>& dof_indices(const int e) const { return d_dof_indices(e); }

    inline const blitz::Array<libMesh::Point,1>& q_point(const int e) const { return d_q_point(e); }
    inline const blitz::Array<double,1>& JxW(const int e) const { return d_JxW(e); }
    inline const blitz::Array<double,2>& phi(const int e) const { return d_phi(e); }
    inline const blitz::Array<double,2>& phi_JxW(const int e) const { return d_phi_JxW(e); }
    inline const blitz::Array<libMesh::VectorValue<double>,2>& dphi(const int e) const { return d_dphi(e); }
    inline const blitz::Array<libMesh::VectorValue<double>,2>& dphi_JxW(const int e) const { return d_dphi_JxW(e); }

    inline const bool& elem_touches_physical_bdry(const int e) const { return d_elem_touches_physical_bdry(e); }
    inline const blitz::Array<bool,1>& elem_side_at_physical_bdry(const int e) const { return d_elem_side_at_physical_bdry(e); }
    inline const blitz::Array<bool,1>& elem_side_at_dirichlet_bdry(const int e) const { return d_elem_side_at_dirichlet_bdry(e); }

    inline const blitz::Array<blitz::Array<libMesh::Point,1>,1>& q_point_face(const int e) const { return d_q_point_face(e); }
    inline const blitz::Array<blitz::Array<double,1>,1>& JxW_face(const int e) const { return d_JxW_face(e); }
    inline const blitz::Array<blitz::Array<libMesh::Point,1>,1>& normal_face(const int e) const { return d_normal_face(e); }
    inline const blitz::Array<blitz::Array<libMesh::Point,1>,1>& normal_JxW_face(const int e) const { return d_normal_JxW_face(e); }
    inline const blitz::Array<blitz::Array<double,2>,1>& phi_face(const int e) const { return d_phi_face(e); }
    inline const blitz::Array<blitz::Array<double,2>,1>& phi_JxW_face(const int e) const { return d_phi_JxW_face(e); }
    inline const blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>& dphi_face(const int e) const { return d_dphi_face(e); }
    inline const blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>& dphi_JxW_face(const int e) const { return d_dphi_JxW_face(e); }

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FESystemDataCache();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FESystemDataCache(
        const FESystemDataCache& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FESystemDataCache&
    operator=(
        const FESystemDataCache& that);

    /*!
     * System and mesh configuration.
     */
    const libMesh::System& d_system;
    const libMesh::MeshBase& d_mesh;

    /*!
     * \brief Configuration options.
     */
    bool d_compute_q_point;
    bool d_compute_JxW;
    bool d_compute_phi;
    bool d_compute_phi_JxW;
    bool d_compute_dphi;
    bool d_compute_dphi_JxW;

    bool d_compute_q_point_face;
    bool d_compute_JxW_face;
    bool d_compute_normal_face;
    bool d_compute_normal_JxW_face;
    bool d_compute_phi_face;
    bool d_compute_phi_JxW_face;
    bool d_compute_dphi_face;
    bool d_compute_dphi_JxW_face;

    /*!
     * \brief Cached data items.
     */
    blitz::Array<libMesh::Elem*,1> d_elems;
    blitz::Array<blitz::Array<std::vector<unsigned int>,1>,1> d_dof_indices;

    blitz::Array<blitz::Array<libMesh::Point,1>,1> d_q_point;
    blitz::Array<blitz::Array<double,1>,1> d_JxW;
    blitz::Array<blitz::Array<double,2>,1> d_phi;
    blitz::Array<blitz::Array<double,2>,1> d_phi_JxW;
    blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1> d_dphi;
    blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1> d_dphi_JxW;

    blitz::Array<bool,1> d_elem_touches_physical_bdry;
    blitz::Array<blitz::Array<bool,1>,1> d_elem_side_at_physical_bdry;
    blitz::Array<blitz::Array<bool,1>,1> d_elem_side_at_dirichlet_bdry;

    blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1> d_q_point_face;
    blitz::Array<blitz::Array<blitz::Array<double,1>,1>,1> d_JxW_face;
    blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1> d_normal_face;
    blitz::Array<blitz::Array<blitz::Array<libMesh::Point,1>,1>,1> d_normal_JxW_face;
    blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1> d_phi_face;
    blitz::Array<blitz::Array<blitz::Array<double,2>,1>,1> d_phi_JxW_face;
    blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1> d_dphi_face;
    blitz::Array<blitz::Array<blitz::Array<libMesh::VectorValue<double>,2>,1>,1> d_dphi_JxW_face;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/FESystemDataCache.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FESystemDataCache
