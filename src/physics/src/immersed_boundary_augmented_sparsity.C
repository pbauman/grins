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
#include "grins/immersed_boundary_augmented_sparsity.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/overlapping_fluid_solid_map.h"

// libMesh
#include "libmesh/mesh_base.h"

namespace GRINS
{
  void ImmersedBoundaryAugmentedSparsity::augment_sparsity_pattern( libMesh::SparsityPattern::Graph &,
                                                                    std::vector<libMesh::dof_id_type> & n_nz,
                                                                    std::vector<libMesh::dof_id_type> & n_oz )
  {
    const libMesh::MeshBase & mesh = _system.get_mesh();

    const std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > > &
      fluid_map = _fluid_solid_overlap.fluid_map();

    // Loop over all the fluid elements that have overlapping solid elements
    for( std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > >::const_iterator
           fluid_elem_map_it = fluid_map.begin();
         fluid_elem_map_it != fluid_map.end();
         ++fluid_elem_map_it )
      {
        libMesh::dof_id_type fluid_elem_id = fluid_elem_map_it->first;

        const libMesh::Elem & fluid_elem = mesh.elem_ref(fluid_elem_id);

        const std::map<libMesh::dof_id_type,std::vector<unsigned int> > & solid_elem_map =
          fluid_elem_map_it->second;

        // Extract solid element ids for this fluid element
        std::vector<libMesh::dof_id_type> solid_elem_ids;
        solid_elem_ids.reserve(solid_elem_map.size());

        for( std::map<libMesh::dof_id_type,std::vector<unsigned int> >::const_iterator
               solid_map_it = solid_elem_map.begin();
             solid_map_it != solid_elem_map.end();
             ++solid_map_it )
          solid_elem_ids.push_back(solid_map_it->first);

        // Now set the sparsity for this set of overlappling fluid/solid elements
        this->set_sparsity(_system,fluid_elem,solid_elem_ids,n_nz,n_oz);
      }
  }

  void ImmersedBoundaryAugmentedSparsity::set_sparsity
  ( MultiphysicsSystem & system,
    const libMesh::Elem & fluid_elem, const std::vector<libMesh::dof_id_type> & solid_elem_ids,
    std::vector<libMesh::dof_id_type> & n_nz, std::vector<libMesh::dof_id_type> & n_oz )
  {
    const libMesh::DofMap & dof_map = system.get_dof_map();
    const libMesh::MeshBase & mesh = _system.get_mesh();

    std::vector<libMesh::dof_id_type> u_fluid_dofs, v_fluid_dofs, w_fluid_dofs;
    dof_map.dof_indices( &fluid_elem, u_fluid_dofs, _flow_vars.u() );
    dof_map.dof_indices( &fluid_elem, v_fluid_dofs, _flow_vars.v() );
    if( _flow_vars.dim() == 3)
      dof_map.dof_indices( &fluid_elem, w_fluid_dofs, _flow_vars.w() );

    std::vector<libMesh::dof_id_type> fluid_dofs;
    fluid_dofs.reserve( u_fluid_dofs.size()+v_fluid_dofs.size()+w_fluid_dofs.size() );

    fluid_dofs.insert( fluid_dofs.end(), u_fluid_dofs.begin(), u_fluid_dofs.end() );
    fluid_dofs.insert( fluid_dofs.end(), v_fluid_dofs.begin(), v_fluid_dofs.end() );
    if( _flow_vars.dim() == 3)
      fluid_dofs.insert( fluid_dofs.end(), w_fluid_dofs.begin(), w_fluid_dofs.end() );

    std::vector<libMesh::dof_id_type> u_solid_dofs, v_solid_dofs, w_solid_dofs, solid_dofs;

    for( std::vector<libMesh::dof_id_type>::const_iterator id = solid_elem_ids.begin();
         id < solid_elem_ids.end(); ++id )
      {
        const libMesh::Elem & solid_elem = mesh.elem_ref(*id);

        u_solid_dofs.clear();
        v_solid_dofs.clear();
        w_solid_dofs.clear();

        dof_map.dof_indices( &solid_elem, u_solid_dofs, _disp_vars.u() );
        dof_map.dof_indices( &solid_elem, v_solid_dofs, _disp_vars.v() );
        if( _flow_vars.dim() == 3)
          dof_map.dof_indices( &solid_elem, w_solid_dofs, _disp_vars.w() );

        solid_dofs.clear();
        solid_dofs.reserve( u_solid_dofs.size()+v_solid_dofs.size()+w_solid_dofs.size() );

        solid_dofs.insert( solid_dofs.end(), u_solid_dofs.begin(), u_solid_dofs.end() );
        solid_dofs.insert( solid_dofs.end(), v_solid_dofs.begin(), v_solid_dofs.end() );
        if( _flow_vars.dim() == 3)
          solid_dofs.insert( solid_dofs.end(), w_solid_dofs.begin(), w_solid_dofs.end() );

        // Add the fluid->solid coupling
        for( unsigned int i = 0; i < fluid_dofs.size(); i++ )
          {
            // Only deal with augmenting the sparsity pattern if the fluid dof is local
            if( this->dof_is_local(dof_map,fluid_dofs[i]) )
              {
                unsigned int dof_offset = fluid_dofs[i] - dof_map.first_dof();

                libmesh_assert_less(dof_offset, n_nz.size());
                libmesh_assert_less(dof_offset, n_oz.size());

                unsigned int n_local_coupled_dofs = 0;
                unsigned int n_remote_coupled_dofs = 0;

                for( unsigned int j = 0; j < solid_dofs.size(); j++ )
                  {
                    if( this->dof_is_local(dof_map,solid_dofs[j]) )
                      n_local_coupled_dofs++;
                    else
                      n_remote_coupled_dofs++;
                  }

                n_nz[dof_offset] += n_local_coupled_dofs;
                n_oz[dof_offset] += n_remote_coupled_dofs;
              }
          }

        // Now add the solid -> fluid coupling
        for( unsigned int i = 0; i < solid_dofs.size(); i++ )
          {
            // Only deal with augmenting the sparsity pattern if the fluid dof is local
            if( this->dof_is_local(dof_map,solid_dofs[i]) )
              {
                unsigned int dof_offset = solid_dofs[i] - dof_map.first_dof();

                libmesh_assert_less(dof_offset, n_nz.size());
                libmesh_assert_less(dof_offset, n_oz.size());

                unsigned int n_local_coupled_dofs = 0;
                unsigned int n_remote_coupled_dofs = 0;

                for( unsigned int j = 0; j < fluid_dofs.size(); j++ )
                  {
                    if( this->dof_is_local(dof_map,fluid_dofs[j]) )
                      n_local_coupled_dofs++;
                    else
                      n_remote_coupled_dofs++;
                  }

                n_nz[dof_offset] += n_local_coupled_dofs;
                n_oz[dof_offset] += n_remote_coupled_dofs;
              }
          }

      } // end loop over solid elems
  }

} // end namespace GRINS
