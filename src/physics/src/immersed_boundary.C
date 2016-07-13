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
      _solid_mech(solid_mech_ptr.release())
  {
    // Get the fluid mechanics from the input, this is needed to query the fluid subdomain ids
    bool have_fluid_mech = input.have_variable( "Physics/ImmersedBoundary/fluid_mechanics" );
    if( !have_fluid_mech )
      libmesh_error_msg("Error: Must specify fluid_mechanics to identify the fluid in the IBM physics.\n");

    else
      this->_fluid_mechanics = input( "Physics/ImmersedBoundary/fluid_mechanics", "DIE!" );

    // Get the solid mechanics from the input, this is needed to prerequest data in init_context
    bool have_solid_mech = input.have_variable( "Physics/ImmersedBoundary/solid_mechanics" );
    if( !have_solid_mech )
      libmesh_error_msg("Error: Must specify solid_mechanics to identify the solid in the IBM physics.\n");

    else
      this->_solid_mechanics = input( "Physics/ImmersedBoundary/solid_mechanics", "DIE!" );

    // Get the solid subdomain ids from the input
    unsigned int n_solid_subdomains = input.vector_variable_size( "Physics/ImmersedBoundary/enabled_subdomains" );
    if( n_solid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify enabled_subdomain ids to identify the solid in the IBM physics.\n");

    else
      {
        std::string solid_id_loc = "Physics/ImmersedBoundary/enabled_subdomains";

        for( unsigned int solid_dom_idx = 0; solid_dom_idx < n_solid_subdomains; solid_dom_idx++)
          _solid_subdomain_set.insert( input(solid_id_loc, -1, solid_dom_idx) );
      }

    // Get the subdomain ids for the fluid
    unsigned int n_fluid_subdomains = input.vector_variable_size( "Physics/"+ _fluid_mechanics  +"/enabled_subdomains" );
    if( n_fluid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify enabled_subdomains for the fluid. \n");

    else
      {
        std::string fluid_id_loc = "Physics/"+ _fluid_mechanics + "/enabled_subdomains";
        for( unsigned int fluid_dom_idx = 0; fluid_dom_idx < n_fluid_subdomains; fluid_dom_idx++)
          _fluid_subdomain_set.insert( input(fluid_id_loc, -1, fluid_dom_idx) );
      }

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_variables( libMesh::FEMSystem* system )
  {
    // Get the point locator object that will find the right fluid element
    // TODO: move this to somewhere not so hackish
    _point_locator = system->get_mesh().sub_point_locator();

    // Build helper FEMContexts. We'll use this to handle
    // supplemental finite element data for variables that
    // we need, but are not defined on the "current" subdomain.
    {
      libMesh::UniquePtr<libMesh::DiffContext> raw_context = system->build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      _solid_context.reset(context);
    }

    {
      libMesh::UniquePtr<libMesh::DiffContext> raw_context = system->build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      _fluid_context.reset(context);
    }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march velocity and displacements forward in time
    system->time_evolving(_flow_vars.u());
    system->time_evolving(_flow_vars.v());

    if ( _flow_vars.dim() == 3 )
      {
        system->time_evolving(_flow_vars.w());
      }
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
        // Get reference to system
        libMesh::System & base_system = const_cast<libMesh::System &>(fluid_context.get_system());
        libMesh::FEMSystem & system = libMesh::cast_ref<libMesh::FEMSystem &>(base_system);

        const SolidElemToQpMap & solid_elem_map = fluid_elem_map_it->second;

        for( SolidElemToQpMap::const_iterator solid_elem_map_it = solid_elem_map.begin();
             solid_elem_map_it != solid_elem_map.end();
             ++solid_elem_map_it )
          {
            libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
            libMesh::FEMContext & solid_context = libMesh::cast_ref<libMesh::FEMContext &>(*raw_context.get());

            // Prepare and reinit helper solid FEMContext for the solid element.
            libMesh::dof_id_type solid_elem_id = solid_elem_map_it->first;
            const libMesh::Elem* solid_elem = system.get_mesh().elem(solid_elem_id);

            const std::vector<unsigned int> & solid_qp_indices = solid_elem_map_it->second;


            const std::vector<libMesh::Point>& all_solid_qpoints =
              solid_context.get_element_fe(_disp_vars.u(),2)->get_xyz();

            const std::vector<std::vector<libMesh::Real> > & solid_phi =
              solid_context.get_element_fe(_disp_vars.u(),2)->get_phi();

             const std::vector<libMesh::Real> & solid_JxW =
               solid_context.get_element_fe(_disp_vars.u(),2)->get_JxW();

            solid_context.pre_fe_reinit(system,solid_elem);
            solid_context.elem_fe_reinit();

            const unsigned int n_solid_dofs =
              solid_context.get_dof_indices(this->_disp_vars.u()).size();

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


            // Now assemble fluid part of velocity matching term
            // into solid residual.
            libMesh::DenseSubVector<libMesh::Number> & Fus =
              solid_context.get_elem_residual(this->_disp_vars.u());

            libMesh::DenseSubVector<libMesh::Number> * Fvs = NULL;
            libMesh::DenseSubVector<libMesh::Number> * Fws = NULL;

            if ( this->_disp_vars.dim() >= 2 )
              Fvs = &solid_context.get_elem_residual(this->_disp_vars.v());

            if ( this->_disp_vars.dim() == 3 )
              Fws = &solid_context.get_elem_residual(this->_disp_vars.w());

            unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

            libMesh::DenseMatrix<libMesh::Number> K;
            libMesh::DenseSubMatrix<libMesh::Number> Kus_uf(K), Kvs_vf(K), Kws_wf(K);

            if( compute_jacobian)
              {
                K.resize( this->_disp_vars.dim()*n_solid_dofs, this->_flow_vars.dim()*n_fluid_dofs );

                Kus_uf.reposition( this->_disp_vars.u()*n_solid_dofs, this->_flow_vars.u()*n_fluid_dofs,
                                   n_solid_dofs, n_fluid_dofs );

                if( this->_disp_vars.dim() >= 2 )
                  {
                    Kvs_vf.reposition( this->_disp_vars.v()*n_solid_dofs, this->_flow_vars.v()*n_fluid_dofs,
                                       n_solid_dofs, n_fluid_dofs );
                  }
                if( this->_disp_vars.dim() == 3 )
                  {
                    Kws_wf.reposition( this->_disp_vars.w()*n_solid_dofs, this->_flow_vars.w()*n_fluid_dofs,
                                       n_solid_dofs, n_fluid_dofs );
                  }
              }

            unsigned int n_qpoints = solid_qpoints.size();

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
                                           solid_context.get_dof_indices(),
                                           fluid_context.get_dof_indices() );
              }
            else
              {
                system.get_dof_map().constrain_element_vector
                  ( solid_context.get_elem_residual(),
                    solid_context.get_dof_indices(), false );
              }

            system.rhs->add_vector( solid_context.get_elem_residual(),
                                    solid_context.get_dof_indices() );

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

    // Get reference to system
    libMesh::System & base_system = const_cast<libMesh::System &>(context.get_system());
    libMesh::FEMSystem & system = libMesh::cast_ref<libMesh::FEMSystem &>(base_system);

    const unsigned int n_solid_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    // Global coordinates of the solid qp points
    const std::vector<libMesh::Point> & solid_qpoints =
      context.get_element_fe(this->_disp_vars.u(),2)->get_xyz();

    const std::vector<libMesh::Real> & solid_JxW =
      context.get_element_fe(this->_disp_vars.u(),2)->get_JxW();

    const std::vector<std::vector<libMesh::Real> > & solid_phi =
      context.get_element_fe(this->_disp_vars.u(),2)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & solid_dphi =
      context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

    libMesh::UniquePtr<libMesh::DiffContext> raw_context = system.build_context();
    libMesh::FEMContext & fluid_context = libMesh::cast_ref<libMesh::FEMContext &>(*raw_context.get());

    const std::vector<std::vector<libMesh::RealGradient> > & fluid_dphi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // Solid residuals
    libMesh::DenseSubVector<libMesh::Number> & Fus = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> * Fvs = NULL;
    libMesh::DenseSubVector<libMesh::Number> * Fws = NULL;

    // Solid Jacobians
    libMesh::DenseSubMatrix<libMesh::Number>& Kus_us =
      context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());

    libMesh::DenseSubMatrix<libMesh::Number>* Kvs_vs = NULL;
    libMesh::DenseSubMatrix<libMesh::Number>* Kws_ws = NULL;

    if ( this->_disp_vars.dim() >= 2 )
      {
        Fvs = &context.get_elem_residual(this->_disp_vars.v());
        Kvs_vs = &context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());
      }
    if ( this->_disp_vars.dim() == 3 )
      {
        Fws = &context.get_elem_residual(this->_disp_vars.w());
        Kws_ws = &context.get_elem_jacobian(this->_disp_vars.w(),this->_disp_vars.w());
      }

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

    libMesh::DenseMatrix<libMesh::Number> Kmat;

    for(unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // We only do this if the time solver is unsteady. Otherwise,
        // the solution_rate is not allocated.
        if( !system.get_time_solver().is_steady() )
          {
            // Velocity of solid at quadrature point.
            libMesh::Point Udot;
            context.interior_rate(this->_disp_vars.u(), qp, Udot(0));
            if ( this->_disp_vars.dim() >= 2 )
              context.interior_rate(this->_disp_vars.v(), qp, Udot(1));
            if ( this->_disp_vars.dim() == 3 )
              context.interior_rate(this->_disp_vars.w(), qp, Udot(2));

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
                          context.get_elem_solution_rate_derivative();

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
        context.interior_value(this->_disp_vars.u(), qp, U(0));
        if ( this->_disp_vars.dim() >= 2 )
          context.interior_value(this->_disp_vars.v(), qp, U(1));
        if ( this->_disp_vars.dim() == 3 )
          context.interior_value(this->_disp_vars.w(), qp, U(2));

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

        libMesh::DenseSubVector<libMesh::Number> & Fu = fluid_context.get_elem_residual(this->_flow_vars.u());
        libMesh::DenseSubVector<libMesh::Number> & Fv = fluid_context.get_elem_residual(this->_flow_vars.v());
        libMesh::DenseSubVector<libMesh::Number> * Fw = NULL;

        if( _flow_vars.dim() == 3 )
          Fw = &fluid_context.get_elem_residual(this->_flow_vars.w());

        libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kmat), Kuf_vs(Kmat), Kuf_ws(Kmat);
        libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kmat), Kvf_vs(Kmat), Kvf_ws(Kmat);
        libMesh::DenseSubMatrix<libMesh::Number> Kwf_us(Kmat), Kwf_vs(Kmat), Kwf_ws(Kmat);

        if( compute_jacobian)
          {
            Kmat.resize( this->_flow_vars.dim()*n_fluid_dofs, this->_disp_vars.dim()*n_solid_dofs );

            Kuf_us.reposition( this->_flow_vars.u()*n_fluid_dofs, this->_disp_vars.u()*n_solid_dofs,
                               n_fluid_dofs, n_solid_dofs );

            Kvf_us.reposition( this->_flow_vars.v()*n_fluid_dofs, this->_disp_vars.u()*n_solid_dofs,
                               n_fluid_dofs, n_solid_dofs );

            if( _flow_vars.dim() == 3 )
              Kwf_us.reposition( this->_flow_vars.w()*n_fluid_dofs, this->_disp_vars.u()*n_solid_dofs,
                                 n_fluid_dofs, n_solid_dofs );

            if( this->_disp_vars.dim() >= 2 )
              {
                Kuf_vs.reposition( this->_flow_vars.u()*n_fluid_dofs, this->_disp_vars.v()*n_solid_dofs,
                                   n_fluid_dofs, n_solid_dofs );

                Kvf_vs.reposition( this->_flow_vars.v()*n_fluid_dofs, this->_disp_vars.v()*n_solid_dofs,
                                   n_fluid_dofs, n_solid_dofs );

                if( _flow_vars.dim() == 3 )
                  Kwf_vs.reposition( this->_flow_vars.w()*n_fluid_dofs, this->_disp_vars.v()*n_solid_dofs,
                                     n_fluid_dofs, n_solid_dofs );
              }

            if( this->_disp_vars.dim() == 3 )
              {
                Kuf_ws.reposition( this->_flow_vars.u()*n_fluid_dofs, this->_flow_vars.w()*n_solid_dofs,
                                   n_fluid_dofs, n_solid_dofs );

                Kvf_ws.reposition( this->_flow_vars.v()*n_fluid_dofs, this->_flow_vars.w()*n_solid_dofs,
                                   n_fluid_dofs, n_solid_dofs );

                if( _flow_vars.dim() == 3 )
                  Kwf_ws.reposition( this->_flow_vars.w()*n_fluid_dofs, this->_flow_vars.w()*n_solid_dofs,
                                     n_fluid_dofs, n_solid_dofs );
              }
          }

        libMesh::Gradient grad_u, grad_v, grad_w;
        context.interior_gradient(this->_disp_vars.u(), qp, grad_u);
        context.interior_gradient(this->_disp_vars.v(), qp, grad_v);

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
      } // end quadrature point loop

    // Since we manually built the fluid context, we have to manually
    // constrain and add the residuals and Jacobians.
    /*! \todo  We're hardcoding to the case that the residual is always
      assembled and homogeneous constraints. */
    if( compute_jacobian )
      {
        system.get_dof_map().constrain_element_matrix_and_vector
          ( Kmat,
            fluid_context.get_elem_residual(),
            fluid_context.get_dof_indices(), false );

        system.matrix->add_matrix( Kmat,
                                   fluid_context.get_dof_indices(),
                                   solid_context.get_dof_indices() );
      }
    else
      {
        system.get_dof_map().constrain_element_vector
          ( fluid_context.get_elem_residual(),
            fluid_context.get_dof_indices(), false );
      }

    system.rhs->add_vector( fluid_context.get_elem_residual(),
                            fluid_context.get_dof_indices() );

    if( fluid_context.get_elem_residual().l2_norm() > 0.0 )
      std::cout << "elem " << context.get_elem().id()
                << ", fluid residual = "
                << fluid_context.get_elem_residual()
                << std::endl;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative_fluid( AssemblyContext & context )
  {
    // Since IBM only acts as a source on the fluid, need to inverse map solid quad points onto the
    // fluid elements and reinit the fluid fe with the mapped quad points
    libMesh::PointLocatorBase &  lctr = *_point_locator; //locator object

    // Global coordinates of the solid qp points
    const std::vector<libMesh::Point> & solid_qp_pts = context.get_element_fe( this->_disp_vars.u() )->get_xyz();

    const unsigned int n_qpoints = context.get_element_qrule().n_points();

    // The main qp loop which calculates residual contributions
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        // The shape functions need to be evaluated at the displaced body.
        // first get the displacement of the solid variables
        libMesh::Number x_disp = 0;
        libMesh::Number y_disp = 0;
        libMesh::Number z_disp = 0;
        context.interior_value(this->_disp_vars.u(), qp, x_disp);
        if ( this->_disp_vars.dim() >= 2 )
          context.interior_value(this->_disp_vars.v(), qp, y_disp);
        if ( this->_disp_vars.dim() == 3 )
          context.interior_value(this->_disp_vars.w(), qp, z_disp);
        libMesh::Point qp_disp(x_disp,y_disp,z_disp);

        // Get a pointer to the fluid element that contains the displaced solid qp
        const libMesh::Elem * fluid_elem_qp = lctr( solid_qp_pts[qp] + qp_disp , &_fluid_subdomain_set );
        libmesh_assert(fluid_elem_qp); //ensure we find an element

        // TODO: the dimension shouldnt be hardcoded as 2 for the fluid element
        libMesh::UniquePtr<libMesh::FEGenericBase<libMesh::Real> > fluid_element_fe(libMesh::FEGenericBase<libMesh::Real>::build(2, context.get_element_fe(_flow_vars.u())->get_fe_type()));
        const std::vector<std::vector<libMesh::RealGradient> > & fluid_dphi = fluid_element_fe->get_dphi();

        // Reinit the fluid element with the mapped solid quadrature points
        fluid_element_fe->reinit(fluid_elem_qp, &solid_qp_pts);

        // Jacobian values vector at the quadrature points
        // TODO: hardcoded dim 2
        const std::vector<libMesh::Real> & JxW = context.get_element_fe(_disp_vars.u(),2)->get_JxW();

        const unsigned int fl_dof_size = libMesh::FEInterface::n_dofs(2, fluid_element_fe->get_fe_type(), fluid_elem_qp->type());

        /*
        // Build up info for fluid residual vectors
        std::vector<libMesh::dof_id_type> fluid_u_dof_indices;
        std::vector<libMesh::dof_id_type> fluid_v_dof_indices;
        std::vector<libMesh::dof_id_type> fluid_w_dof_indices;
        context.get_system().get_dof_map().dof_indices(fluid_elem_qp, fluid_u_dof_indices, _flow_vars.u());
        context.get_system().get_dof_map().dof_indices(fluid_elem_qp, fluid_v_dof_indices, _flow_vars.v());
        if (this->_flow_vars.dim() == 3)
          context.get_system().get_dof_map().dof_indices(fluid_elem_qp, fluid_w_dof_indices, _flow_vars.w());

        libMesh::DenseVector<libMesh::Number> Fu_new(fl_dof_size);
        libMesh::DenseVector<libMesh::Number> Fv_new(fl_dof_size);
        libMesh::DenseVector<libMesh::Number> Fw_new(fl_dof_size);
        */

        // Fluid Residuals that we're populating
        libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(this->_flow_vars.u());
        libMesh::DenseSubVector<libMesh::Number> & Fv = context.get_elem_residual(this->_flow_vars.v());
        libMesh::DenseSubVector<libMesh::Number> * Fw = NULL;

        // Gradients w.r.t. the master element coordinates
        libMesh::Gradient grad_u,grad_v,grad_w;
        _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

        // Piola-kirchoff stress tensor in the reference configuration
        // TODO: tau needs to be scaled basd on mesh dimension
        libMesh::TensorValue<libMesh::Real> tau;
        ElasticityTensor C;
        _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);

        for (unsigned int i=0; i != fl_dof_size; i++)
          {
            // Tau:grad_phi source term contributions
            for (unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++)
              {
                for (unsigned int beta = 0; beta < this->_disp_vars.dim(); beta++)
                  {
                    //Fu_new(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);
                    //Fv_new(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);

                    Fu(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);
                    Fv(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);

                    if (this->_flow_vars.dim() == 3)
                      {
                        //Fw_new(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);
                        (*Fw)(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);

                      }
                  }
              }
          } //fluid dof loop

        /*
        // Now  constrain the fluid residual vectors and then add them to the system
        libMesh::FEMSystem & femsys = libMesh::cast_ref<libMesh::FEMSystem&>(const_cast<libMesh::System&>(context.get_system()));

        context.get_system().get_dof_map().constrain_element_vector(Fu_new, fluid_u_dof_indices, false);
        context.get_system().get_dof_map().constrain_element_vector(Fv_new, fluid_v_dof_indices, false);
        if ( this->_flow_vars.dim() == 3 )
          context.get_system().get_dof_map().constrain_element_vector(Fw_new, fluid_w_dof_indices, false);

        femsys.rhs->add_vector(Fu_new, fluid_u_dof_indices);
        femsys.rhs->add_vector(Fv_new, fluid_v_dof_indices);
        if ( this->_flow_vars.dim() == 3 )
          femsys.rhs->add_vector(Fw_new, fluid_w_dof_indices);
        */
      }// qp loop
  } //element_time_derivative_fluid


  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative_solid( AssemblyContext& context )
  {
    // Solid variables evolve via the underlying velocity of the displaced solid dofs
    // Contributions depend on the dimensionality of the solid body

    libMesh::PointLocatorBase &  lctr = *_point_locator; //locator object

    libMesh::Elem * solid_elem = &context.get_elem();

    // Solid residuals
    libMesh::DenseSubVector<libMesh::Number> & u_disp = context.get_elem_residual(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> * v_disp = NULL;
    libMesh::DenseSubVector<libMesh::Number> * w_disp = NULL;

    // Only get the solid residuals for the dimension we need
    if ( this->_disp_vars.dim() >= 2 )
      v_disp = &context.get_elem_residual(this->_disp_vars.v());
    if ( this->_disp_vars.dim() == 3 )
      w_disp = &context.get_elem_residual(this->_disp_vars.w());

    const unsigned int sl_dof_size = context.get_dof_indices(_disp_vars.u()).size();
    std::vector<libMesh::dof_id_type> solid_dof_indices;
    context.get_system().get_dof_map().dof_indices(solid_elem, solid_dof_indices, _disp_vars.u());

    //for (unsigned int i=0; i != solid_dof_indices.size(); i++)
    for (unsigned int i=0; i != sl_dof_size; i++)
      {
        //std::cout << "the solid elem id is " << solid_elem->id() << std::endl;

        libMesh::Point solid_dof_loc = solid_elem->point(i);

        // How much each dof moved
        libMesh::Number dof_x_disp = 0.;
        libMesh::Number dof_y_disp = 0.;
        libMesh::Number dof_z_disp = 0.;

        // Find how much the solid dofs have displaced
        // Note: constrained dofs dont move
        if (! context.get_system().get_dof_map().is_constrained_dof( solid_dof_indices[i] ))
          {
            dof_x_disp = context.point_value(this->_disp_vars.u(), solid_dof_loc);
            if (this->_disp_vars.dim() >= 2)
              dof_y_disp = context.point_value(this->_disp_vars.v(), solid_dof_loc);
            if (this->_disp_vars.dim() == 3)
              dof_z_disp = context.point_value(this->_disp_vars.w(), solid_dof_loc);
          }

        libMesh::Point new_dof_loc(solid_dof_loc(0) + dof_x_disp, solid_dof_loc(1) + dof_y_disp, solid_dof_loc(2) + dof_z_disp);
        std::vector<libMesh::Point> new_dof_loc_vec;
        new_dof_loc_vec.push_back(new_dof_loc);

        // Get a pointer to the fluid element that contains the displaced solid dof
        const libMesh::Elem * fluid_elem_dof = lctr( new_dof_loc , &_fluid_subdomain_set );
        libmesh_assert(fluid_elem_dof); // ensure we found an element

        // TODO: the dimension shouldnt be hardcoded as 2 for the fluid element
        libMesh::UniquePtr<libMesh::FEGenericBase<libMesh::Real>> fluid_element_fe_dof(libMesh::FEGenericBase<libMesh::Real>::build(2, context.get_element_fe(_flow_vars.u())->get_fe_type()));
        const std::vector<std::vector<libMesh::Real> > & fl_phi = fluid_element_fe_dof->get_phi();

        // Reinit the fluid element to obtain its velocities to supply the solid displacements
        fluid_element_fe_dof->reinit(fluid_elem_dof, &new_dof_loc_vec);

        // Fluid velocities (to be computed) at displaced solid dofs
        libMesh::Number fl_vel_x=0.;
        libMesh::Number fl_vel_y=0.;
        libMesh::Number fl_vel_z=0.;

        // TODO: hard coded dim 2 for the  fluid
        const unsigned int n_fl_dofs = libMesh::FEInterface::n_dofs(2, fluid_element_fe_dof->get_fe_type(), fluid_elem_dof->type());

        std::vector<libMesh::dof_id_type> fluid_dof_indices_u;
        std::vector<libMesh::dof_id_type> fluid_dof_indices_v;
        std::vector<libMesh::dof_id_type> fluid_dof_indices_w;

        std::vector<libMesh::Number> ucoef;
        context.get_system().get_dof_map().dof_indices(fluid_elem_dof,fluid_dof_indices_u, _flow_vars.u());
        context.get_system().current_local_solution->get(fluid_dof_indices_u, ucoef);

        for (unsigned int l = 0; l != n_fl_dofs; l++)
          {
            fl_vel_x += fl_phi[l][0] * ucoef[l];
          }
        if (this->_flow_vars.dim() >= 2 )
          {
            std::vector<libMesh::Number> vcoef;
            context.get_system().get_dof_map().dof_indices(fluid_elem_dof,fluid_dof_indices_v, _flow_vars.v());
            context.get_system().current_local_solution->get(fluid_dof_indices_v, vcoef);
            for (unsigned int l = 0; l != n_fl_dofs; l++)
              {
                fl_vel_y += fl_phi[l][0] * vcoef[l];
              }
          }
        if (this->_flow_vars.dim() == 3 )
          {
            std::vector<libMesh::Number> wcoef;
            context.get_system().get_dof_map().dof_indices(fluid_elem_dof,fluid_dof_indices_w, _flow_vars.w());
            context.get_system().current_local_solution->get(fluid_dof_indices_w, wcoef);
            for (unsigned int l = 0; l != n_fl_dofs; l++)
              {
                fl_vel_z += fl_phi[l][0] * wcoef[l];
              }
          }

        // Finally, set the computed velocities as the displacements of the solid
        //u_disp(i) = (context.get_elem_solution_rate(_disp_vars.u())(i) - fl_vel_x);
        u_disp(i) = fl_vel_x;
        if (this->_disp_vars.dim() >= 2 )
          //(*v_disp)(i) = (context.get_elem_solution_rate(_disp_vars.v())(i) - fl_vel_y);
          (*v_disp)(i) = fl_vel_y;
        if (this->_disp_vars.dim() == 3 )
          //(*w_disp)(i) = (context.get_elem_solution_rate(_disp_vars.w())(i) - fl_vel_z);
          (*w_disp)(i) = fl_vel_z;

      } //solid dof loop
  } //element_time_derivative_solid

    //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

} // namespace GRINS
