//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// this class
#include "grins/immersed_boundary.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/elasticity_tensor.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"

// includes for IBM instantiation
#include "grins/elastic_membrane.h"
#include "grins/elastic_cable.h"
#include "grins/hookes_law_1d.h"
#include "grins/hookes_law.h"
#include "grins/mooney_rivlin.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"


// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

namespace GRINS
{

  template<typename SolidMech>
  ImmersedBoundary<SolidMech>::ImmersedBoundary(const std::string& physics_name, libMesh::UniquePtr<SolidMech> & solid_mech_ptr, const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::physics_velocity_variable_name(input,physics_name))),
      _disp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<DisplacementVariable>(VariablesParsing::physics_disp_variable_name(input,physics_name))),
      _solid_mech(solid_mech_ptr.release()),
      _fluid_mechanics(input("Physics/ImmersedBoundary/fluid_mechanics","DIE!")),
      _solid_mechanics(input("Physics/ImmersedBoundary/solid_mechanics","DIE!"))
  {
    // Get the fluid mechanics from the input, this is needed to query the fluid subdomain ids
    if( !input.have_variable( "Physics/ImmersedBoundary/fluid_mechanics" ) )
      libmesh_error_msg("Error: Must specify fluid_mechanics to identify the fluid in the IBM physics.\n");

    // Get the solid mechanics from the input, this is needed to prerequest data in init_context
    if( !input.have_variable( "Physics/ImmersedBoundary/solid_mechanics" ) )
      libmesh_error_msg("Error: Must specify solid_mechanics to identify the solid in the IBM physics.\n");

    // Get the solid subdomain ids from the input
    std::string solid_id_str = "Physics/"+_solid_mechanics+"/enabled_subdomains";
    unsigned int n_solid_subdomains = input.vector_variable_size(solid_id_str);
    if( n_solid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify "+solid_id_str+" to identify the solid in the IBM physics.\n");

    for( unsigned int i = 0; i < n_solid_subdomains; i++)
      _solid_subdomain_set.insert( input(solid_id_str, -1, i) );

    // Get the subdomain ids for the fluid
    std::string fluid_id_str = "Physics/"+_fluid_mechanics+"/enabled_subdomains";
    unsigned int n_fluid_subdomains = input.vector_variable_size(fluid_id_str);
    if( n_fluid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify "+fluid_id_str+" to identify the fluid in the IBM physics.\n");

    for( unsigned int i = 0; i < n_fluid_subdomains; i++)
      _fluid_subdomain_set.insert( input(fluid_id_str, -1, i) );

    // TODO: Need to check that Mesh has all the fluid and solid subdomain ids
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::auxiliary_init( MultiphysicsSystem & system )
  {
    // Get the point locator object that will find the right fluid element
    _point_locator = system.get_mesh().sub_point_locator();

    // Build helper FEMContexts. We'll use this to handle
    // supplemental finite element data for variables that
    // we need, but are not defined on the "current" subdomain.
    {
      libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      _solid_context.reset(context);
    }

    {
      libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      _fluid_context.reset(context);
    }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Velocity are first order in time
    system->time_evolving(_flow_vars.u(),1);
    system->time_evolving(_flow_vars.v(),1);

    if ( _flow_vars.dim() == 3 )
      system->time_evolving(_flow_vars.w(),1);


    // Displacements are second order in time
    system->time_evolving(_disp_vars.u(),2);

    if( _disp_vars.dim() >= 2 )
      system->time_evolving(_disp_vars.v(),2);

    if ( _disp_vars.dim() == 3 )
      system->time_evolving(_disp_vars.w(),2);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.

    if (_solid_mechanics == "ElasticMembrane")
      {
        //membrane context inits
        context.get_element_fe(_disp_vars.u(),2)->get_JxW();
        context.get_element_fe(_disp_vars.u(),2)->get_phi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphidxi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphideta();
        //for metric tensors of the membrane
        context.get_element_fe(_disp_vars.u(),2)->get_dxyzdxi();
        context.get_element_fe(_disp_vars.u(),2)->get_dxyzdeta();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidx();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidy();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidz();
        context.get_element_fe(_disp_vars.u(),2)->get_detadx();
        context.get_element_fe(_disp_vars.u(),2)->get_detady();
        context.get_element_fe(_disp_vars.u(),2)->get_detadz();
      }
    else if (_solid_mechanics == "ElasticCable")
      {
        //cable context inits
        context.get_element_fe(_disp_vars.u(),1)->get_JxW();
        context.get_element_fe(_disp_vars.u(),1)->get_phi();
        context.get_element_fe(_disp_vars.u(),1)->get_dphidxi();
        //for metric tensors of cable
        context.get_element_fe(_disp_vars.u(),1)->get_dxyzdxi();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidx();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidy();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidz();
      }
    else
      {
        std::string err = "ERROR: solid_mechanics not properly specified";
        libmesh_error_msg(err);
      }

    context.get_element_fe(_flow_vars.u())->get_JxW();

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::preassembly( MultiphysicsSystem & system )
  {
    const libMesh::MeshBase & mesh = system.get_mesh();

    for( std::set<libMesh::subdomain_id_type>::const_iterator solid_id_it = _solid_subdomain_set.begin();
         solid_id_it != _solid_subdomain_set.end(); ++solid_id_it )
      for( libMesh::MeshBase::const_element_iterator e = mesh.active_local_subdomain_elements_begin(*solid_id_it);
           e != mesh.active_local_subdomain_elements_end(*solid_id_it);
           ++e )
        {
          // Convenience
          const libMesh::Elem * elem = *e;

          libmesh_assert( is_solid_elem(elem->subdomain_id()) );

          const std::vector<libMesh::Point>& qpoints =
            _solid_context->get_element_fe(_disp_vars.u(),2)->get_xyz();

          _solid_context->get_element_fe(_disp_vars.u(),2)->get_phi();

          _solid_context->pre_fe_reinit(system,elem);
          _solid_context->elem_fe_reinit();

          // Find what fluid element contains each of the quadrature points and cache
          for( unsigned int qp = 0; qp < qpoints.size(); qp++ )
            {
              libMesh::Real u_disp = 0;
              libMesh::Real v_disp = 0;
              libMesh::Real w_disp = 0;

              _solid_context->interior_value(this->_disp_vars.u(), qp, u_disp);
              if( this->_disp_vars.dim() >= 2 )
                _solid_context->interior_value(this->_disp_vars.v(), qp, v_disp);
              if( this->_disp_vars.dim() == 3 )
                _solid_context->interior_value(this->_disp_vars.w(), qp, w_disp);

              libMesh::Point U( u_disp, v_disp, w_disp );

              // We need to look for overlap with *displaced* solid point
              libMesh::Point x = qpoints[qp]+U;

              const libMesh::Elem * fluid_elem =
                (*_point_locator)( x, &_fluid_subdomain_set );

              if( !fluid_elem )
                libmesh_error_msg("ERROR: Could not find fluid element for given displacement! Likely solid displaced off the fluid mesh!");

              SolidElemToQpMap & solid_elem_map = _fluid_id_to_solid_ids_qps[fluid_elem->id()];

              std::vector<unsigned int>& solid_qps = solid_elem_map[elem->id()];
              solid_qps.push_back(qp);
            }
        }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext & context,
                                                             CachedValues & /*cache*/ )
  {
    if( this->is_fluid_elem( context.get_elem().subdomain_id() ) )
      this->assemble_fluid_var_residual_contributions(compute_jacobian,context);

    if( this->is_solid_elem( context.get_elem().subdomain_id() ) )
      this->assemble_solid_var_residual_contributions(compute_jacobian,context);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::assemble_fluid_var_residual_contributions( bool compute_jacobian,
                                                                               AssemblyContext & context )
  {
    // For clarity
    AssemblyContext & fluid_context = context;

    FluidToSolidMap::const_iterator fluid_elem_map_it =
      _fluid_id_to_solid_ids_qps.find(fluid_context.get_elem().id());

    // Check if this fluid element had any solid elements overlapping
    if( fluid_elem_map_it != _fluid_id_to_solid_ids_qps.end() )
      {
        MultiphysicsSystem & system = context.get_multiphysics_system();

        const SolidElemToQpMap & solid_elem_map = fluid_elem_map_it->second;

        for( SolidElemToQpMap::const_iterator solid_elem_map_it = solid_elem_map.begin();
             solid_elem_map_it != solid_elem_map.end();
             ++solid_elem_map_it )
          {
            libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
            AssemblyContext & solid_context = libMesh::cast_ref<AssemblyContext &>(*raw_context.get());

            // Prepare and reinit helper solid FEMContext for the solid element.
            libMesh::dof_id_type solid_elem_id = solid_elem_map_it->first;
            const libMesh::Elem* solid_elem = system.get_mesh().elem(solid_elem_id);

            const std::vector<unsigned int> & solid_qp_indices = solid_elem_map_it->second;

            unsigned int us_var = this->_disp_vars.u();

            const std::vector<libMesh::Point>& all_solid_qpoints =
              solid_context.get_element_fe(us_var,2)->get_xyz();

            const std::vector<std::vector<libMesh::Real> > & solid_phi =
              solid_context.get_element_fe(us_var,2)->get_phi();

             const std::vector<libMesh::Real> & solid_JxW =
               solid_context.get_element_fe(us_var,2)->get_JxW();

            solid_context.pre_fe_reinit(system,solid_elem);
            solid_context.elem_fe_reinit();

            const unsigned int n_solid_dofs =
              solid_context.get_dof_indices(us_var).size();

            // Extract the subset of solid quadrature points that overlap with
            // this fluid element
            std::vector<libMesh::Point> solid_qpoints(solid_qp_indices.size());
            for( unsigned int qp = 0; qp < solid_qp_indices.size(); qp++ )
              solid_qpoints[qp] = all_solid_qpoints[ solid_qp_indices[qp] ];

            // Now construct a fluid finite element using the
            // solid element quadrature points and
            libMesh::FEType fluid_fe_type = fluid_context.get_element_fe(_flow_vars.u())->get_fe_type();

            const libMesh::Elem & fluid_elem = fluid_context.get_elem();

            libMesh::UniquePtr<libMesh::FEGenericBase<libMesh::Real> > fluid_fe =
              libMesh::FEGenericBase<libMesh::Real>::build(2,fluid_fe_type);

            const std::vector<std::vector<libMesh::Real> > & fluid_phi = fluid_fe->get_phi();

            fluid_fe->reinit(&fluid_elem,&solid_qpoints);

            unsigned int u_dot_var = system.get_second_order_dot_var(us_var);

            // Now assemble fluid part of velocity matching term
            // into solid residual.
            libMesh::DenseSubVector<libMesh::Number> & Fus =
              solid_context.get_elem_residual(u_dot_var);

            libMesh::DenseSubVector<libMesh::Number> * Fvs = NULL;
            libMesh::DenseSubVector<libMesh::Number> * Fws = NULL;

            unsigned int v_var = libMesh::invalid_uint;
            unsigned int v_dot_var = libMesh::invalid_uint;
            if ( this->_disp_vars.dim() >= 2 )
              {
                v_var = this->_disp_vars.v();
                v_dot_var = system.get_second_order_dot_var(v_var);

                Fvs = &solid_context.get_elem_residual(v_dot_var);
              }

            unsigned int w_var = libMesh::invalid_uint;
            unsigned int w_dot_var = libMesh::invalid_uint;
            if ( this->_disp_vars.dim() == 3 )
              {
                w_var = this->_disp_vars.w();
                w_dot_var = system.get_second_order_dot_var(w_var);

                Fws = &solid_context.get_elem_residual(w_dot_var);
              }

            // Build up solid dof indices
            unsigned int sshift = 1;
            if( us_var != u_dot_var )
              sshift = 2;

            std::vector<libMesh::dof_id_type> solid_dof_indices(_disp_vars.dim()*sshift*n_solid_dofs);
            std::vector<libMesh::dof_id_type>::iterator sdof_start = solid_dof_indices.begin();
            solid_dof_indices.insert( sdof_start,
                                      solid_context.get_dof_indices(us_var).begin(),
                                      solid_context.get_dof_indices(us_var).end() );

            solid_dof_indices.insert( sdof_start+n_solid_dofs,
                                      solid_context.get_dof_indices(v_var).begin(),
                                      solid_context.get_dof_indices(v_var).end() );

            if( us_var != u_dot_var )
              {
                solid_dof_indices.insert( sdof_start+2*n_solid_dofs,
                                          solid_context.get_dof_indices(u_dot_var).begin(),
                                          solid_context.get_dof_indices(u_dot_var).end() );

                solid_dof_indices.insert( sdof_start+3*n_solid_dofs,
                                          solid_context.get_dof_indices(v_dot_var).begin(),
                                          solid_context.get_dof_indices(v_dot_var).end() );
              }

            unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

            std::vector<libMesh::dof_id_type> velocity_dof_indices(_flow_vars.dim()*n_fluid_dofs);
            std::vector<libMesh::dof_id_type>::iterator vdof_start = velocity_dof_indices.begin();
            velocity_dof_indices.insert( vdof_start,
                                         fluid_context.get_dof_indices(this->_flow_vars.u()).begin(),
                                         fluid_context.get_dof_indices(this->_flow_vars.u()).end() );

            velocity_dof_indices.insert( vdof_start+n_fluid_dofs,
                                         fluid_context.get_dof_indices(this->_flow_vars.v()).begin(),
                                         fluid_context.get_dof_indices(this->_flow_vars.v()).end() );

            libMesh::DenseMatrix<libMesh::Number> K;
            libMesh::DenseSubMatrix<libMesh::Number> Kus_uf(K), Kvs_vf(K), Kws_wf(K);

            if( compute_jacobian)
              {
                K.resize( this->_disp_vars.dim()*sshift*n_solid_dofs, this->_flow_vars.dim()*n_fluid_dofs );

                if( us_var == u_dot_var )
                  Kus_uf.reposition( 0, 0, n_solid_dofs, n_fluid_dofs );
                else
                  Kus_uf.reposition( 2*n_solid_dofs, 0, n_solid_dofs, n_fluid_dofs );

                if( this->_disp_vars.dim() >= 2 )
                  {
                    if( us_var == u_dot_var )
                      Kvs_vf.reposition( n_solid_dofs, n_fluid_dofs, n_solid_dofs, n_fluid_dofs );
                    else
                      Kvs_vf.reposition( 3*n_solid_dofs, n_fluid_dofs, n_solid_dofs, n_fluid_dofs );

                  }
                if( this->_disp_vars.dim() == 3 )
                  {
                    if( us_var == u_dot_var )
                      Kws_wf.reposition( 2*n_solid_dofs, 2*n_fluid_dofs, n_solid_dofs, n_fluid_dofs );
                    else
                      Kws_wf.reposition( 5*n_solid_dofs, 2*n_fluid_dofs, n_solid_dofs, n_fluid_dofs );
                  }
              }

            unsigned int n_qpoints = solid_qpoints.size();

            std::cout << "computing fluid part" << std::endl;

            for( unsigned int qp = 0; qp < n_qpoints; qp++ )
              {
                libMesh::Real Vx = 0.0;
                libMesh::Real Vy = 0.0;
                libMesh::Real Vz = 0.0;

                // Compute the fluid velocity at the solid element quadrature points.
                // The fluid solution dofs come from the incoming context. This way,
                // we can still have correct numerical Jacobians.
                for( unsigned int i = 0; i < n_fluid_dofs; i++ )
                  {
                    Vx += fluid_phi[i][qp]*(fluid_context.get_elem_solution(this->_flow_vars.u()))(i);

                    if ( this->_disp_vars.dim() >= 2 )
                      Vy += fluid_phi[i][qp]*(fluid_context.get_elem_solution(this->_flow_vars.v()))(i);

                    if ( this->_disp_vars.dim() == 3 )
                      Vz += fluid_phi[i][qp]*(fluid_context.get_elem_solution(this->_flow_vars.w()))(i);
                  }

                // We're working with a subset of the solid quadrature points
                // so grab the index to the current one
                unsigned int solid_qp_idx = solid_qp_indices[qp];

                libMesh::Real jac = solid_JxW[solid_qp_idx];

                for( unsigned int i = 0; i < n_solid_dofs; i++ )
                  {
                    /*! \todo [PB]: For manifolds, I don't think these are the correct shape
                                    functions because there are missing terms. */
                    Fus(i) += Vx*solid_phi[i][solid_qp_idx]*jac;

                    if ( this->_disp_vars.dim() >= 2 )
                      (*Fvs)(i) += Vy*solid_phi[i][solid_qp_idx]*jac;

                    if ( this->_disp_vars.dim() ==3 )
                      (*Fws)(i) += Vz*solid_phi[i][solid_qp_idx]*jac;

                    if( compute_jacobian )
                      {
                        for( unsigned int j = 0; j < n_fluid_dofs; j++ )
                          {
                            libMesh::Real diag_value =
                              fluid_phi[j][qp]*solid_phi[i][solid_qp_idx]*jac*
                              fluid_context.get_elem_solution_derivative();

                            Kus_uf(i,j) += diag_value;

                            if ( this->_disp_vars.dim() >= 2 )
                              Kvs_vf(i,j) += diag_value;

                            if ( this->_disp_vars.dim() == 3 )
                              Kws_wf(i,j) += diag_value;
                          }
                      }
                  }

              } // end loop over active solid quadrature points


            // Since we manually build the solid context, we have to manually
            // constrain and add the residuals and Jacobians. This
            /*! \todo  We're hardcoding to the case that the residual is always
              assembled and homogeneous constraints. */
            if( compute_jacobian )
              {
                system.get_dof_map().constrain_element_matrix_and_vector
                  ( K,
                    solid_context.get_elem_residual(),
                    solid_context.get_dof_indices(), false );

                system.matrix->add_matrix( K,
                                           solid_dof_indices,
                                           velocity_dof_indices );
              }
            else
              {
                system.get_dof_map().constrain_element_vector
                  ( solid_context.get_elem_residual(),
                    solid_context.get_dof_indices(), false );
              }

            system.rhs->add_vector( solid_context.get_elem_residual(),
                                    solid_dof_indices );

            if( solid_context.get_elem_residual().l2_norm() > 0.0 )
              std::cout << "elem " << context.get_elem().id()
                        << ", solid residual = "
                        << solid_context.get_elem_residual()
                        << std::endl;

        } // end loop over solid elem map
      } // end if fluid element has overlapping solid elem
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::assemble_solid_var_residual_contributions( bool compute_jacobian,
                                                                               AssemblyContext & context )
  {
    // For clarity
    AssemblyContext & solid_context = context;

    MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_var = this->_disp_vars.u();

    const unsigned int n_solid_dofs = solid_context.get_dof_indices(u_var).size();

    std::cout << "n_solid_dofs = " << n_solid_dofs << std::endl;

    // Global coordinates of the solid qp points
    const std::vector<libMesh::Point> & solid_qpoints =
      solid_context.get_element_fe(u_var,2)->get_xyz();

    const std::vector<libMesh::Real> & solid_JxW =
      solid_context.get_element_fe(u_var,2)->get_JxW();

    const std::vector<std::vector<libMesh::Real> > & solid_phi =
      solid_context.get_element_fe(u_var,2)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & solid_dphi =
      solid_context.get_element_fe(u_var,2)->get_dphi();

    libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
    libMesh::FEMContext & fluid_context = libMesh::cast_ref<libMesh::FEMContext &>(*raw_context.get());

    const std::vector<std::vector<libMesh::RealGradient> > & fluid_dphi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);

    // Solid residuals
    libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> * Fvs = NULL;
    libMesh::DenseSubVector<libMesh::Number> * Fws = NULL;

    // Solid Jacobians
    libMesh::DenseSubMatrix<libMesh::Number>& Kus_us =
      solid_context.get_elem_jacobian(u_dot_var,u_var);

    libMesh::DenseSubMatrix<libMesh::Number>* Kvs_vs = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kws_ws = NULL;

    unsigned int v_var = libMesh::invalid_uint;
    unsigned int v_dot_var = libMesh::invalid_uint;
    if ( this->_disp_vars.dim() >= 2 )
      {
        v_var = this->_disp_vars.v();
        v_dot_var = system.get_second_order_dot_var(v_var);

        Fvs = &solid_context.get_elem_residual(v_dot_var);
        Kvs_vs = &solid_context.get_elem_jacobian(v_dot_var,v_var);
      }

    unsigned int w_var = libMesh::invalid_uint;
    unsigned int w_dot_var = libMesh::invalid_uint;
    if ( this->_disp_vars.dim() == 3 )
      {
        w_var = this->_disp_vars.w();
        w_dot_var = system.get_second_order_dot_var(w_var);

        Fws = &solid_context.get_elem_residual(w_var);
        Kws_ws = &solid_context.get_elem_jacobian(w_dot_var,w_var);
      }

    const unsigned int n_qpoints = solid_context.get_element_qrule().n_points();

    libMesh::DenseMatrix<libMesh::Number> Kmat;

    libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kmat), Kuf_vs(Kmat), Kuf_ws(Kmat);
    libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kmat), Kvf_vs(Kmat), Kvf_ws(Kmat);
    libMesh::DenseSubMatrix<libMesh::Number> Kwf_us(Kmat), Kwf_vs(Kmat), Kwf_ws(Kmat);

    std::vector<libMesh::dof_id_type> velocity_dof_indices;

    std::cout << "computing solid  part" << std::endl;

    for(unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // We only do this if the time solver is unsteady. Otherwise,
        // the solution_rate is not allocated.
        if( !system.get_time_solver().is_steady() )
          {
            // Velocity of solid at quadrature point.
            libMesh::Point Udot;
            solid_context.interior_rate(this->_disp_vars.u(), qp, Udot(0));
            if ( this->_disp_vars.dim() >= 2 )
              solid_context.interior_rate(this->_disp_vars.v(), qp, Udot(1));
            if ( this->_disp_vars.dim() == 3 )
              solid_context.interior_rate(this->_disp_vars.w(), qp, Udot(2));

            // Velocity matching term, loop over solid dofs. These
            // are minus since the fluid velocity part was positive.
            for( unsigned int i = 0; i < n_solid_dofs; i++ )
              {
                Fus(i) -= Udot(0)*solid_phi[i][qp]*solid_JxW[qp];

                if( this->_disp_vars.dim() >= 2 )
                  (*Fvs)(i) -= Udot(1)*solid_phi[i][qp]*solid_JxW[qp];

                if( this->_disp_vars.dim() == 3 )
                  (*Fws)(i) -= Udot(2)*solid_phi[i][qp]*solid_JxW[qp];

                if(compute_jacobian)
                  {
                    for( unsigned int j = 0; j < n_solid_dofs; j++ )
                      {
                        libMesh::Real diag_value =
                          solid_phi[j][qp]*solid_phi[i][qp]*solid_JxW[qp]*
                          solid_context.get_elem_solution_rate_derivative();

                        Kus_us(i,j) -= diag_value;

                        if ( this->_disp_vars.dim() >= 2 )
                          (*Kvs_vs)(i,j) -= diag_value;

                        if ( this->_disp_vars.dim() == 3 )
                          (*Kws_ws)(i,j) -= diag_value;
                      }
                  }
              }
          } // end is time_solver steady

        // Quadrature point in the initial configuration
        const libMesh::Point & X = solid_qpoints[qp];

        // Displacement at quadrature point
        libMesh::Point U;
        solid_context.interior_value(this->_disp_vars.u(), qp, U(0));
        if ( this->_disp_vars.dim() >= 2 )
          solid_context.interior_value(this->_disp_vars.v(), qp, U(1));
        if ( this->_disp_vars.dim() == 3 )
          solid_context.interior_value(this->_disp_vars.w(), qp, U(2));

        // The fluid shape functions need to be evaluated at the *displaced*
        // location of the solid.
        const libMesh::Elem * fluid_elem = (*_point_locator)(X+U,&_fluid_subdomain_set);

        if( !fluid_elem )
          libmesh_error_msg("ERROR: Could not find fluid element for given displacement! Likely solid displaced off the fluid mesh!");

        // Need fluid shape function gradients at the current quadrature point of the deformed solid
        std::vector<libMesh::Point> this_qp(1,solid_qpoints[qp]);

        fluid_context.pre_fe_reinit(system,fluid_elem);
        fluid_context.elem_fe_reinit(&this_qp);

        unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

        std::cout << "fluid_dim = " << _flow_vars.dim() << std::endl;
        std::cout << "n_fluid_dofs = " << n_fluid_dofs << std::endl;

        std::cout << "n_ufluid_dofs = " << fluid_context.get_dof_indices(this->_flow_vars.u()).size() << std::endl;
        std::cout << "n_vfluid_dofs = " << fluid_context.get_dof_indices(this->_flow_vars.v()).size() << std::endl;

        velocity_dof_indices.clear();
        velocity_dof_indices.resize(_flow_vars.dim()*n_fluid_dofs);

        std::cout << "vdof size = " << velocity_dof_indices.size() << std::endl;

        std::vector<libMesh::dof_id_type>::iterator vdof_start = velocity_dof_indices.begin();
        const std::vector<libMesh::dof_id_type>& u_dof_indices = fluid_context.get_dof_indices(this->_flow_vars.u());

        for( unsigned int i = 0; i < u_dof_indices.size(); i++ )
          velocity_dof_indices[i] = u_dof_indices[i];

        const std::vector<libMesh::dof_id_type>& v_dof_indices = fluid_context.get_dof_indices(this->_flow_vars.u());

        for( unsigned int i = 0; i < v_dof_indices.size(); i++ )
          velocity_dof_indices[i+n_fluid_dofs] = v_dof_indices[i];

        libMesh::DenseSubVector<libMesh::Number> & Fu = fluid_context.get_elem_residual(this->_flow_vars.u());
        libMesh::DenseSubVector<libMesh::Number> & Fv = fluid_context.get_elem_residual(this->_flow_vars.v());
        libMesh::DenseSubVector<libMesh::Number> * Fw = NULL;

        if( _flow_vars.dim() == 3 )
          Fw = &fluid_context.get_elem_residual(this->_flow_vars.w());


        if( compute_jacobian)
          {
            std::cout << "flow dim = " << this->_flow_vars.dim() << ", solid_dim = " << this->_disp_vars.dim() << std::endl;

            Kmat.resize( this->_flow_vars.dim()*n_fluid_dofs, this->_disp_vars.dim()*n_solid_dofs );

            // We need to manually manage the indexing since we're working only on this particular subblock
            Kuf_us.reposition( 0, 0, n_fluid_dofs, n_solid_dofs );

            Kvf_us.reposition( n_fluid_dofs, 0, n_fluid_dofs, n_solid_dofs );

            if( _flow_vars.dim() == 3 )
              Kwf_us.reposition( 2*n_fluid_dofs, 0, n_fluid_dofs, n_solid_dofs );

            if( this->_disp_vars.dim() >= 2 )
              {
                Kuf_vs.reposition( 0, n_solid_dofs, n_fluid_dofs, n_solid_dofs );

                Kvf_vs.reposition( n_fluid_dofs, n_solid_dofs, n_fluid_dofs, n_solid_dofs );

                if( _flow_vars.dim() == 3 )
                  Kwf_vs.reposition( 2*n_fluid_dofs, n_solid_dofs, n_fluid_dofs, n_solid_dofs );
              }

            if( this->_disp_vars.dim() == 3 )
              {
                Kuf_ws.reposition( 0, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );

                Kvf_ws.reposition( n_fluid_dofs, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );

                if( _flow_vars.dim() == 3 )
                  Kwf_ws.reposition( 2*n_fluid_dofs, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
              }
          }

        libMesh::Gradient grad_u, grad_v, grad_w;
        solid_context.interior_gradient(this->_disp_vars.u(), qp, grad_u);
        solid_context.interior_gradient(this->_disp_vars.v(), qp, grad_v);

        libMesh::TensorValue<libMesh::Real> grad_U;
        grad_U(0,0) = grad_u(0);
        grad_U(0,1) = grad_u(1);
        grad_U(0,2) = grad_u(2);
        grad_U(1,0) = grad_v(0);
        grad_U(1,1) = grad_v(1);
        grad_U(1,2) = grad_v(2);

        libMesh::TensorValue<libMesh::Real> F(grad_U);
        F(0,0) += 1;
        F(1,1) += 1;

        // We need to use F^T a few times so just cache it.
        libMesh::TensorValue<libMesh::Real> Ftrans = F.transpose();

        libMesh::TensorValue<libMesh::Real> E(Ftrans*F);
        E(0,0) -= 1;
        E(1,1) -= 1;
        E *= 0.5;

        if( E.norm() > 1.0e-5 )
          std::cout << "E = " << E << std::endl;

        libMesh::Real Em = 1000000000;
        libMesh::Real nu = 0.3;
        libMesh::Real lambda = Em*nu/(1+nu)*(1-2*nu);
        libMesh::Real mu = Em/(2*(1+nu));

        libMesh::Real trE = E.tr();
        libMesh::TensorValue<libMesh::Real> S(2.0*E*mu);
        S(0,0) += lambda*trE;
        S(1,1) += lambda*trE;

        libMesh::TensorValue<libMesh::Real> P(F*S);

        // The F^T comes from needing the derivative of the fluid
        // shape function w.r.t. solid coordinates
        libMesh::TensorValue<libMesh::Real> tau(P*Ftrans);

        // Gradients w.r.t. the master element coordinates
        _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

        // Piola-kirchoff stress tensor in the reference configuration
        // TODO: tau needs to be scaled basd on mesh dimension
        libMesh::TensorValue<libMesh::Real> blah;
        ElasticityTensor C;
        _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,blah,C);

        /*
        if( tau.norm() > 1.0e-3 )
          std::cout << "grad_u = " << grad_u << std::endl
                    << "grad_v = " << grad_v << std::endl
                    << "tau = " << tau << std::endl;
        */

        for (unsigned int i=0; i != n_fluid_dofs; i++)
          {
            if( fluid_dphi[i].size() != 1 )
              libmesh_error_msg("ERROR: Unexpected fluid_dphi size!");

            // Zero index for fluid dphi/JxW since we only requested one quad. point.
            for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
              {
                Fu(i) -= tau(alpha,0)*fluid_dphi[i][0](alpha)*solid_JxW[qp];
                Fv(i) -= tau(alpha,1)*fluid_dphi[i][0](alpha)*solid_JxW[qp];

                if (this->_flow_vars.dim() == 3)
                  (*Fw)(i) -= tau(alpha,2)*fluid_dphi[i][qp](alpha)*solid_JxW[qp];
              }

            if( compute_jacobian )
              {
                libMesh::Real factor = solid_JxW[qp]*solid_context.get_elem_solution_derivative();

                for (unsigned int j=0; j != n_solid_dofs; j++)
                  {
                    libMesh::Real term1 = (((S*Ftrans).transpose())*solid_dphi[j][qp])*fluid_dphi[i][0];

                    libMesh::Real term3 = (P*solid_dphi[j][qp])*fluid_dphi[i][0];

                    Kuf_us(i,j) -= (term1+term3)*factor;
                    Kvf_us(i,j) -= (term1+term3)*factor;

                    if(this->_flow_vars.dim() == 3)
                      Kwf_us(i,j) -= (term1+term3)*factor;

                    if( this->_disp_vars.dim() >= 2 )
                      {
                        Kuf_vs(i,j) -= (term1+term3)*factor;
                        Kvf_vs(i,j) -= (term1+term3)*factor;

                        if(this->_flow_vars.dim() == 3)
                          Kwf_vs(i,j) -= (term1+term3)*factor;
                      }

                    if( this->_disp_vars.dim() == 3 )
                      {
                        Kuf_ws(i,j) -= (term1+term3)*factor;
                        Kvf_ws(i,j) -= (term1+term3)*factor;

                        if(this->_flow_vars.dim() == 3)
                          Kwf_ws(i,j) -= (term1+term3)*factor;
                      }

                    for( unsigned int I = 0; I < _disp_vars.dim(); I++ )
                      for( unsigned int J = 0; J < _disp_vars.dim(); J++ )
                        for( unsigned int K = 0; K < _disp_vars.dim(); K++ )
                          for( unsigned int L = 0; L < _disp_vars.dim(); L++ )
                            {
                              Kuf_us(i,j) -= F(0,I)*C(I,J,K,L)*F(0,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);
                              Kvf_us(i,j) -= F(1,I)*C(I,J,K,L)*F(0,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);

                              if(this->_flow_vars.dim() == 3)
                                Kwf_us(i,j) -= F(2,I)*C(I,J,K,L)*F(0,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);

                              if( this->_disp_vars.dim() >= 2 )
                                {
                                  Kuf_vs(i,j) -= F(0,I)*C(I,J,K,L)*F(1,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);
                                  Kvf_vs(i,j) -= F(1,I)*C(I,J,K,L)*F(1,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);

                                  if(this->_flow_vars.dim() == 3)
                                    Kwf_vs(i,j) -= F(2,I)*C(I,J,K,L)*F(1,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);
                                }

                              if( this->_disp_vars.dim() == 3 )
                                {
                                  Kuf_ws(i,j) -= F(0,I)*C(I,J,K,L)*F(2,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);
                                  Kvf_ws(i,j) -= F(1,I)*C(I,J,K,L)*F(2,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);

                                  if(this->_flow_vars.dim() == 3)
                                    Kwf_ws(i,j) -= F(2,I)*C(I,J,K,L)*F(2,K)*(Ftrans.row(J)*fluid_dphi[i][0])*solid_dphi[j][qp](L);
                                }
                            }

                  }
              }// compute_jacobian

          } //fluid dof loop

        // Since we manually built the fluid context, we have to manually
        // constrain and add the residuals and Jacobians.
        /*! \todo  We're hardcoding to the case that the residual is always
          assembled and homogeneous constraints. */
        if( compute_jacobian )
          {
            system.get_dof_map().constrain_element_matrix
              ( Kmat,
                velocity_dof_indices,
                solid_context.get_dof_indices(), false );

            system.matrix->add_matrix( Kmat,
                                       velocity_dof_indices,
                                       solid_context.get_dof_indices() );
          }

        system.get_dof_map().constrain_element_vector
          ( fluid_context.get_elem_residual(),
            fluid_context.get_dof_indices(), false );

        system.rhs->add_vector( fluid_context.get_elem_residual(),
                                fluid_context.get_dof_indices() );

        if( fluid_context.get_elem_residual().l2_norm() > 0.0 )
          std::cout << "elem " << solid_context.get_elem().id()
                << ", fluid residual = "
                    << fluid_context.get_elem_residual()
                    << std::endl;

      } // end quadrature point loop

  }

    //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

} // namespace GRINS
