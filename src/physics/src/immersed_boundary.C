
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

// for the instantiation
#include "grins/elastic_cable.h"
#include "grins/elastic_membrane.h"
#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/hyperelasticity.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"

namespace GRINS
{

  template<typename SolidMech>
  ImmersedBoundary<SolidMech>::ImmersedBoundary(const GRINS::PhysicsName& physics_name, const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(input,physics_name),
      _disp_vars(input,physics_name,false,true,false) /* is2d, is3d, is constraint*/
  {
    this->_solid_mech = new SolidMech(physics_name, input, false); /*is_compressible*/

    //get the subdomain id for the solid from the input
    this->_subdomain_id = input.vector_variable_size( "ImmersedBoundary/Solid/enabled_subdomains" );
    if( this->_subdomain_id == 0 )
      {
        std::cerr << "Error: Must specify at least one subdomain id to identify the solid." << std::endl;
        libmesh_error();
      }
    
  }
  
  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_variables( libMesh::FEMSystem* system )
  {
    this->_dim = system->get_mesh().mesh_dimension();
    this->_flow_vars.init(system);
    this->_disp_vars.init(system);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity and displacements forward in time
    system->time_evolving(_flow_vars.u());
    system->time_evolving(_flow_vars.v());

    system->time_evolving(_disp_vars.u());
    system->time_evolving(_disp_vars.v());

    if (dim == 3)
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

    //these arent needed?
    context.get_side_fe(_flow_vars.u())->get_JxW();
    context.get_side_fe(_flow_vars.u())->get_phi();
    context.get_side_fe(_flow_vars.u())->get_dphi();
    context.get_side_fe(_flow_vars.u())->get_xyz();
  }
  
  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext& context,
                                                             CachedValues& /*cache*/ )
  {
    //Gather info we might need for ibm 
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_qpoints = context.get_element_qrule().n_points();
    const std::vector<libMesh::Real> &JxW = context.get_element_fe(this->_flow_vars.u())->get_JxW();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_flow_vars.w());

    // The velocity shape functions at interior quadrature points. (we only need their gradients prol)
    const std::vector<std::vector<libMesh::Real> >& u_phi = context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The velocity shape function gradients (in global coords.) at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi = context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // these shape functions wont exist in the solid. I will need to locate the fluid element on which the solid element lies

    if (context.get_elem().subdomain_id() == this->_subdomain_id)
      {
        for (unsigned int qp=0; qp != n_qpoints; qp++)
          {
            libMesh::Real jac = JxW[qp];

            //these are gradients w.r.t. the master element coordinates
            libMesh::Gradient grad_u,grad_v,grad_w;
            _solid_mech->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

            libMesh::TensorValue<libMesh::Real> tau;// piola kirchoff stress tensor in the reference configuration
            ElasticityTensor C;
            // get the stress tensor for the solid physics
            _solid_mech->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);
                        
            //how to do a double dot product? by hand!
            for (unsigned int i=0; i != n_u_dofs; i++)
              {

                
                Fu(i) += jac*u_phi[i][qp];

                //tau._coords[0]*u_gradphi[i][qp] + tau._coords[1]*u_gradphi[0][1][qp] + tau._coords[2]*u_gradphi[0][2][qp] + \
                //tau._coords[3]*u_gradphi[1][0][qp] + tau._coords[4]*u_gradphi[1][1][qp] + tau._coords[5]*u_gradphi[1][2][qp] + \
                //tau._coords[6]*u_gradphi[2][0][qp] + tau._coords[7]*u_gradphi[2][1][qp] + tau._coords[8]*u_gradphi[2][2][qp] ;
                Fu(i) += jac*u_phi[i][qp];
                Fv(i) += jac*u_phi[i][qp];

                if (_dim == 3){
                  Fw(i) += jac*u_phi[i][qp];
                }

              }
            //Factor 2: acceleration term: how to get velocity values at a current vs previous point in time
            
          }//qp loop
      }
  } //end elem time derivative
} // namespace GRINS


//Instantiate
template class GRINS::ImmersedBoundary< GRINS::ElasticCable<GRINS::HookesLaw1D> >;
template class GRINS::ImmersedBoundary< GRINS::ElasticMembrane<GRINS::HookesLaw> >;

