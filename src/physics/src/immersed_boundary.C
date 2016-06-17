//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

    std::cout << "Constructing IBM physics" << std::endl;

    // Get the solid subdomain id from the input
    unsigned int n_solid_domains = input.vector_variable_size( "Physics/ImmersedBoundary/enabled_subdomains" );
    if( n_solid_domains == 0 )
      {
        std::string err = "Error: Must specify enabled_subdomains id to identify the solid in the IBM physics.\n";
        libmesh_error_msg(err);
      }
    else
      {
        this->set_parameter(_solid_subdomain_id, input, "Physics/ImmersedBoundary/enabled_subdomains", 1 );
      }

    // Get the fluid mechanics from the input, this is needed to query the fluid subdomain ids
    bool have_fluid_mech = input.have_variable( "Physics/ImmersedBoundary/fluid_mechanics" );
    if( !have_fluid_mech )
      {
        std::string err = "Error: Must specify fluid_mechanics to identify the fluid in the IBM physics.\n";
        libmesh_error_msg(err);
      }
    else
      {
        this->_fluid_mechanics = input( "Physics/ImmersedBoundary/fluid_mechanics", "DIE!" );
      }

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
    _pnt_lctr = system->get_mesh().sub_point_locator();
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
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext& context,
                                                             CachedValues& /*cache*/ )
  {

    //IBM acts as a source term only under solid elements
    if (context.get_elem().subdomain_id() == this->_solid_subdomain_id )
      {
        // Since IBM only acts as a source on the fluid, need to inverse map solid quad points onto the
        // fluid elements and reinit the fluid fe with the mapped quad points
        libMesh::PointLocatorBase &  lctr = *_pnt_lctr; //locator object

        // Global coordinates of the solid qp points
        const std::vector<libMesh::Point> solid_qp_pts = context.get_element_fe( this->_disp_vars.u() )->get_xyz();

        libMesh::Elem * solid_elem = &context.get_elem();

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
            const std::vector<std::vector<libMesh::RealGradient> >& fluid_dphi = fluid_element_fe->get_dphi();

            // Reinit the fluid element with the mapped solid quadrature points
            fluid_element_fe->reinit(fluid_elem_qp, &solid_qp_pts);

            // Jacobian values vector at the quadrature points
            // TODO: hardcoded dim 2
            const std::vector<libMesh::Real> &JxW = context.get_element_fe(_disp_vars.u(),2)->get_JxW();

            // Fluid Residuals that we're populating
            libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(this->_flow_vars.u()); //R_{u}
            libMesh::DenseSubVector<libMesh::Number> & Fv = context.get_elem_residual(this->_flow_vars.v()); //R_{v}
            libMesh::DenseSubVector<libMesh::Number> * Fw = NULL;

            // Solid residuals
            libMesh::DenseSubVector<libMesh::Number> & u_disp = context.get_elem_residual(this->_disp_vars.u());
            libMesh::DenseSubVector<libMesh::Number> * v_disp = NULL;
            libMesh::DenseSubVector<libMesh::Number> * w_disp = NULL;

            // Only get the residuals for the dimension we need
            if ( this->_disp_vars.dim() >= 2 )
              v_disp = &context.get_elem_residual(this->_disp_vars.v());
            if ( this->_disp_vars.dim() == 3 )
              w_disp = &context.get_elem_residual(this->_disp_vars.w());
            if ( this->_flow_vars.dim() == 3 )
              Fw = &context.get_elem_residual(this->_flow_vars.w()); //R_{w}

            // Gradients w.r.t. the master element coordinates
            libMesh::Gradient grad_u,grad_v,grad_w;
            _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

            // Piola-kirchoff stress tensor in the reference configuration
            libMesh::TensorValue<libMesh::Real> tau;
            ElasticityTensor C;
            _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);


            const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                // Tau:grad_phi source term contributions
                for (int alpha = 0; alpha < this->_disp_vars.dim(); alpha++)
                  {
                    for (int beta = 0; beta < this->_disp_vars.dim(); beta++)
                      {
                        Fu(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);
                        Fv(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);

                        if (this->_flow_vars.dim() == 3)
                          {
                            (*Fw)(i) += tau(alpha,beta)*JxW[qp]*fluid_dphi[i][qp](beta);
                          }
                      }
                  }

                // Solid variable contributions
                // Solid variables evolve via the underlying velocity of the displaced solid dofs
                // Contributions depend on the dimensionality of the solid body
                libMesh::Point solid_dof_loc = solid_elem->point(i);

                std::vector<libMesh::dof_id_type> solid_dof_indices;
                context.get_system().get_dof_map().dof_indices(solid_elem,solid_dof_indices);

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
                  libMesh::UniquePtr<libMesh::FEGenericBase<libMesh::Real> > fluid_element_fe_dof(libMesh::FEGenericBase<libMesh::Real>::build(2, context.get_element_fe(_flow_vars.u())->get_fe_type()));
                  const std::vector<std::vector<libMesh::Real> >& fl_phi = fluid_element_fe_dof->get_phi();

                  // Reinit the fluid element to obtain its velocities to supply the solid displacements
                  fluid_element_fe_dof->reinit(fluid_elem_dof, &new_dof_loc_vec);

                  // Fluid velocities (to be computed) at displaced solid dofs
                  libMesh::Number fl_vel_xval=0.;
                  libMesh::Number fl_vel_yval=0.;
                  libMesh::Number fl_vel_zval=0.;

                  // TODO: hard coded dim 2 for the  fluid
                  const unsigned int n_fl_dofs = libMesh::FEInterface::n_dofs(2, fluid_element_fe_dof->get_fe_type(), fluid_elem_dof->type());

                  std::vector<libMesh::dof_id_type> fluid_dof_indices;
                  context.get_system().get_dof_map().dof_indices(fluid_elem_dof,fluid_dof_indices);

                  std::vector<libMesh::Number> ucoef;
                  context.get_system().current_local_solution->get(fluid_dof_indices, ucoef);

                  for (unsigned int l = 0; l != n_fl_dofs; l++)
                    {
                      fl_vel_xval += fl_phi[l][0] * ucoef[l];
                    }

                  if (this->_flow_vars.dim() >= 2 )
                    {
                      std::vector<libMesh::Number> vcoef;
                      context.get_system().current_local_solution->get(fluid_dof_indices, vcoef);
                      for (unsigned int l = 0; l != n_fl_dofs; l++)
                        {
                          fl_vel_yval += fl_phi[l][0] * vcoef[l];
                        }
                    }

                  if (this->_flow_vars.dim() == 3 )
                    {
                      std::vector<libMesh::Number> wcoef;
                      context.get_system().current_local_solution->get(fluid_dof_indices, wcoef);
                      for (unsigned int l = 0; l != n_fl_dofs; l++)
                        {
                          fl_vel_zval += fl_phi[l][0] * wcoef[l];
                        }

                    }

                  // Finally, set the computed velocities as the displacements of the solid
                  u_disp(i) = fl_vel_xval;
                  if (this->_disp_vars.dim() >= 2 )
                    (*v_disp)(i) = fl_vel_yval;
                  if (this->_disp_vars.dim() == 3 )
                    (*w_disp)(i) = fl_vel_zval;

                  } //dof loop

                //TODO : Factor 2: acceleration term: this will be implemented in the mass_residual :
                // how to get velocity values at a current vs previous point in time

                if (compute_jacobian)
                  {
                    //TODO: immersed_boundary jacobian not implemented
                    libmesh_not_implemented();
                  }

              }// qp loop
          } // if elem subdomain id == solid
      } //end elem time derivative


    //instantiate IBM classes
    template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
    template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
    template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

  } // namespace GRINS
