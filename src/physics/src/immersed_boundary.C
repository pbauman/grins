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
      {
        std::string err = "Error: Must specify solid_mechanics to identify the solid in the IBM physics.\n";
        libmesh_error_msg(err);
      }
    else
      {
        this->_solid_mechanics = input( "Physics/ImmersedBoundary/solid_mechanics", "DIE!" );
      }

    // Get the solid subdomain ids from the input
    unsigned int n_solid_subdomains = input.vector_variable_size( "Physics/ImmersedBoundary/enabled_subdomains" );
    if( n_solid_subdomains == 0 )
      {
        std::string err = "Error: Must specify enabled_subdomain ids to identify the solid in the IBM physics.\n";
        libmesh_error_msg(err);
      }
    else
      {
        std::string solid_id_loc = "Physics/ImmersedBoundary/enabled_subdomains";
        for( unsigned int solid_dom_idx = 0; solid_dom_idx < n_solid_subdomains; solid_dom_idx++)
          {
            _solid_subdomain_set.insert( input(solid_id_loc, -1, solid_dom_idx) );
          }
      }

    // Get the subdomain ids for the fluid
    unsigned int n_fluid_subdomains = input.vector_variable_size( "Physics/"+ _fluid_mechanics  +"/enabled_subdomains" );
    if( n_fluid_subdomains == 0 )
      {
        std::string err = "Error: Must specify enabled_subdomains for the fluid. \n";
        libmesh_error_msg(err);
      }
    else
      {
        std::string fluid_id_loc = "Physics/"+ _fluid_mechanics + "/enabled_subdomains";
        for( unsigned int fluid_dom_idx = 0; fluid_dom_idx < n_fluid_subdomains; fluid_dom_idx++)
          {
            _fluid_subdomain_set.insert( input(fluid_id_loc, -1, fluid_dom_idx) );
          }
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

    for( libMesh::MeshBase::const_element_iterator e = mesh.active_local_elements_begin();
         e != mesh.active_local_elements_end();
         ++e )
      {
        // Convenience
        const libMesh::Elem * elem = *e;

        if( is_solid_elem(elem->subdomain_id()) )
          {
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

                SolidElemToQpMap & solid_elem_map = _fluid_id_to_solid_ids_qps[fluid_elem->id()];

                std::vector<unsigned int>& solid_qps = solid_elem_map[elem->id()];
                solid_qps.push_back(qp);
              }
          }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext & context,
                                                             CachedValues & /*cache*/ )
  {
    // is the element a fluid elem?
    if ( this->is_fluid_elem( context.get_elem().subdomain_id() ) )
        element_time_derivative_fluid(context);

    // or is it a solid?
    if ( this->is_solid_elem( context.get_elem().subdomain_id() ) )
        element_time_derivative_solid(context);

    //TODO: immersed_boundary solid jacobian not implemented
    if (compute_jacobian)
        libmesh_not_implemented();

  } // element_time_derivative


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
