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

    std::cout << "constructing ibm physics" << std::endl;
                   
    //get the subdomain id for the solid from the input
    this->_subdomain_id = input.vector_variable_size( "Physics/ImmersedBoundary/enabled_subdomains" );
    if( this->_subdomain_id == 0 )
      {
        std::string err = "Error: Must specify enabled_subdomains id to identify the solid in the IBM physics.\n";
        libmesh_error_msg(err);
      }
    
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

    //ibm acts as a source term only on solid elements 
    if (context.get_elem().subdomain_id() == this->_subdomain_id)
      {

        //Gather general info we might need for ibm 
        const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();
        const unsigned int n_qpoints = context.get_element_qrule().n_points();
        const std::vector<libMesh::Real> &JxW = context.get_element_fe(this->_flow_vars.u())->get_JxW();
        
        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            // since ibm only acts as a source on the fluid, need to inverse map quad points on the
            // fluid elements and reinit the fluid fe with the mapped quad points
            // the regular context.get_element_fe shape functions wont exist in the solid: need to locate the fluid element on which the solid element lies

            // get the point locator object
            libMesh::UniquePtr < libMesh::PointLocatorBase > pnt_lctr = context.get_system().get_mesh().sub_point_locator();
            libMesh::PointLocatorBase & lctr = *pnt_lctr;            

            //get a pointer to the fluid element that contains the solid qp
            //TODO: this call should also provide the fluid subdomain id in order to speed up the search
            const libMesh::Elem * fluid_elem = lctr( context.get_elem().point(qp) ); 
            
            libMesh::FEBase * fluid_element_fe;
            context.get_element_fe<libMesh::Real>(this->_flow_vars.u(), fluid_element_fe);

            //reinit the fluid element with the mapped quadrature points
            fluid_element_fe->reinit(fluid_elem);
              
            
            // Residuals that we're populating
            libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u());
            libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v());
            libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_flow_vars.w());
            
            // The velocity shape functions at interior quadrature points. (we only need their gradients prolly in elem_time_deriv)
            const std::vector<std::vector<libMesh::Real> >& u_phi = context.get_element_fe(this->_flow_vars.u())->get_phi();



            
            // The velocity shape function gradients (in global coords.) at interior quadrature points.
            const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi = context.get_element_fe(this->_flow_vars.u())->get_dphi();

            libMesh::Real jac = JxW[qp];

            //these are gradients w.r.t. the master element coordinates
            libMesh::Gradient grad_u,grad_v,grad_w;
            _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

            // the piola kirchoff stress tensor in the reference configuration
            libMesh::TensorValue<libMesh::Real> tau;
            ElasticityTensor C;
            _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);
                        
            for (unsigned int i=0; i != n_u_dofs; i++)
              {
                
                for (int alpha = 0; alpha < _dim; alpha++)
                  {
                    for (int beta = 0;beta < _dim; beta++)
                      {
                        Fu(i) += tau(alpha,beta)*jac*u_gradphi[i][qp](alpha);
                        Fv(i) += tau(alpha,beta)*jac*u_gradphi[i][qp](alpha);
                      
                        if (_dim == 3){
                          Fw(i) += tau(alpha,beta)*jac*u_gradphi[i][qp](alpha);
                        }
                      }
                  }
              } //dof loop

            //TODO : Factor 2: acceleration term: how to get velocity values at a current vs previous point in time
            
          }// qp loop
      } // if elem subdomain id is right
  } //end elem time derivative


  //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;
  
  
  
} // namespace GRINS
