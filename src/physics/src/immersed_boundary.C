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
            _fluid_subdomain_set.insert( input("fluid_id_loc", -1, fluid_dom_idx) );
          }
      }
    
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_variables( libMesh::FEMSystem* system )
  {
    // Get the point locator object that will find the right fluid element
    _pnt_lctr = system->get_mesh().sub_point_locator();
  }
  
  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march velocity and displacements forward in time
    system->time_evolving(_flow_vars.u());
    system->time_evolving(_flow_vars.v());

    system->time_evolving(_disp_vars.u());
    system->time_evolving(_disp_vars.v());

    if (_dim == 3)
      {
        system->time_evolving(_flow_vars.w());
        system->time_evolving(_disp_vars.w());
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_flow_vars.u())->get_JxW();
    context.get_element_fe(_flow_vars.u())->get_phi();
    context.get_element_fe(_flow_vars.u())->get_dphi();
    context.get_element_fe(_flow_vars.u())->get_xyz();
  }
  
  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext& context,
                                                             CachedValues& /*cache*/ )
  {   

    //IBM acts as a source term only under solid elements 
    if (context.get_elem().subdomain_id() == this->_solid_subdomain_id )
      {
        //Gather general info for the solid element 
        const unsigned int n_qpoints = context.get_element_qrule().n_points();
        const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();            
        libMesh::PointLocatorBase &  lctr = *_pnt_lctr;
        
        // Since IBM only acts as a source on the fluid, need to inverse map solid quad points onto the
        // fluid elements and reinit the fluid fe with the mapped quad points

        // Global coordinates of the solid qp points
        const  std::vector<libMesh::Point> solid_qp_pts = context.get_element_fe(this->_flow_vars.u())->get_xyz();

        // First loop over solid qp to precompute needed info
        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            //we should precompute gradphi and such here (needs to be inside the init_context )
          }        

        // Now do the real work over the solid qps
        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            //get a pointer to the fluid element that contains the solid qp
            const libMesh::Elem * fluid_elem = lctr( solid_qp_pts[qp], &_fluid_subdomain_set ); 

            libMesh::FEBase * fluid_element_fe;
            context.get_element_fe<libMesh::Real>(this->_flow_vars.u(), fluid_element_fe);

            // Reinit the fluid element with the mapped solid quadrature points
            fluid_element_fe->reinit(fluid_elem, &solid_qp_pts);

            // The velocity shape function gradients (in global coords.) at interior quadrature points of the fluid element
            const std::vector<std::vector<libMesh::Number> >& dphidxi = fluid_element_fe->get_dphidxi();
            const std::vector<std::vector<libMesh::Number> >& dphideta = fluid_element_fe->get_dphideta();
            const std::vector<std::vector<libMesh::Number> >& dphidzeta = fluid_element_fe->get_dphidzeta();

            // Jacobian vector at the quadrature points
            const std::vector<libMesh::Real> &JxW = fluid_element_fe->get_JxW();
                    
            // Residuals that we're populating
            libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u());
            libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v());
            //libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_flow_vars.w());
                      
            //Gradients w.r.t. the master element coordinates
            libMesh::Gradient grad_u,grad_v,grad_w;
            _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

            // Piola-kirchoff stress tensor in the reference configuration
            libMesh::TensorValue<libMesh::Real> tau;
            ElasticityTensor C;
            _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);

            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                // tau:grad_phi source term
                for (int alpha = 0; alpha < _dim; alpha++)
                  {
                    for (int beta = 0;beta < _dim; beta++)
                      {
                        Fu(i) += tau(alpha,beta)*JxW[qp]*(dphidxi[i][qp] + dphideta[i][qp] + dphidzeta[i][qp] );
                        Fv(i) += tau(alpha,beta)*JxW[qp]*(dphidxi[i][qp] + dphideta[i][qp] + dphidzeta[i][qp] );
                      
                        if (_dim == 3){
                          //            Fw(i) += tau(alpha,beta)*JxW[qp]*(dphidxi[i][qp] + dphideta[i][qp] + dphidzeta[i][qp] );
                        }
                      }
                  }
              } //dof loop

            //TODO : Factor 2: acceleration term: how to get velocity values at a current vs previous point in time
          
            if (compute_jacobian)
              {
                //TODO: not implemented
              }
                       
          }// qp loop
      } // if elem subdomain id is right
  } //end elem time derivative


  //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;
  
  
  
} // namespace GRINS
