//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// This class
#include "grins/immersed_boundary_coupling_functor.h"


// GRINS
#include "grins/common.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

namespace GRINS
{
  void ImmersedBoundaryCouplingFunctor::operator()
    ( const libMesh::MeshBase::const_element_iterator & range_begin,
      const libMesh::MeshBase::const_element_iterator & range_end,
      libMesh::processor_id_type p,
      std::unordered_map<const libMesh::Elem *,const libMesh::CouplingMatrix*> & coupled_elements )
  {
    const auto & solid_map = _overlapping_map.solid_map();
    const auto & fluid_map = _overlapping_map.fluid_map();

    for( const auto & elem : libMesh::as_range(range_begin,range_end) )
      {
        libMesh::dof_id_type elem_id = elem->id();

        auto solid_map_it = solid_map.find(elem_id);
        auto fluid_map_it = fluid_map.find(elem_id);

        // If this element is a solid element, then we need to populate
        // the coupled_elements with all the fluid elements associated with
        // this solid element. We use the same coupling matrix for all of them.
        if( solid_map_it != solid_map.end() )
          {
            const auto & fluid_group = solid_map_it->second;

            for( const auto & fluid_it : fluid_group )
              {
                const libMesh::Elem * fluid_elem = _mesh.elem_ptr(fluid_it.first);

                if(!fluid_elem)
                  libmesh_error_msg("ERROR: fluid_elem is NULL!");

                if( fluid_elem->processor_id() != p )
                  coupled_elements.insert( std::make_pair(fluid_elem,&_coupling_matrix) );
              }
          }

        // If this element is a fluid element, then we need to populate
        // the coupled_elements with all the solid elements associated with
        // this fluid element. While we don't need the algebraic coupling
        // this will generate, it is needed to get the correct sparsity
        // pattern.
        if( fluid_map_it != fluid_map.end() )
          {
            const auto & solid_group = fluid_map_it->second;

            for( const auto & solid_it : solid_group )
              {
                const libMesh::Elem * solid_elem = _mesh.elem_ptr(solid_it.first);

                if(!solid_elem)
                  libmesh_error_msg("ERROR: fluid_elem is NULL!");

                if( solid_elem->processor_id() != p )
                  coupled_elements.insert( std::make_pair(solid_elem,&_coupling_matrix) );
              }
          }
      }
  }

} // end namespace GRINS
