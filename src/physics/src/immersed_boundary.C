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


#include "grins/immersed_boundary.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/elasticity_tensor.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// for the instantiation
#include "grins/elastic_cable.h"
#include "grins/elastic_membrane.h"
#include "grins/hookes_law.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"
#include "libmesh/quadrature.h"

namespace GRINS
{

  template<typename SolidMech>
  ImmersedBoundary<SolidMech>::ImmersedBoundary(const GRINS::PhysicsName& physics_name, const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(input,physics_name),
      _disp_vars(input,physics_name,false,true,false) /* is2d, is3d, is constraint (how do we know this?)*/
  {
    this->_solid_mech = new SolidMech(physics_name, input, false); /*is_compressible*/
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
        system->time_evolving(_disp_vars.v());
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
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const std::vector<libMesh::Real> &JxW = this->get_fe(context)->get_JxW();

    // Residuals that we're populating
    libMesh::DenseSubVector<libMesh::Number> &Fu = context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Number> &Fv = context.get_elem_residual(this->_flow_vars.v());
    libMesh::DenseSubVector<libMesh::Number> &Fw = context.get_elem_residual(this->_flow_vars.w());

    // The velocity shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& u_phi =
      context.get_element_fe(this->_flow_vars.u())->get_phi();

    // The velocity shape function gradients (in global coords.)
    // at interior quadrature points.
    const std::vector<std::vector<libMesh::RealGradient> >& u_gradphi =
      context.get_element_fe(this->_flow_vars.u())->get_dphi();

    // these shape functions wont exist in the solid. I will need to locate the fluid element on which the solid element lies  (solid vs fluid is updated using the elem->subdomain_id())
    

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        //check to see if we are part of the solid. Do we check qp dofs the element? I think should be dof
        bool is_solid = 1;
        if (is_solid)
          {

            //Factor 1: fsi term
            // get the stress tensor for the solid physics

            //these are gradients w.r.t. the master element coordinates
            libMesh::Gradient grad_u,grad_v,grad_w;
            this->get_grad_uvw(context, qp, grad_u,grad_v,grad_w);

            libMesh::TensorValue<libMesh::Real> tau;
            ElasticityTensor C;
            this->get_stress_and_elasticity(context,qp,grad_u,grad_v,grad_w,tau,C);
            //this tau is the piola kirch stress tensor in the reference configuration
                        
            //how to do a double dot product? by hand!


            //Factor 2: acceleration term
            //drho density factor how to get dimension of the solid mesh? needed for density (eq 17)
            
            //how to get velocity values at a current vs previous point in time
            
          }//is solid
      }//qp loop
  } //end elem time derivative

} // namespace GRINS


//Instantiate
template class GRINS::ImmersedBoundary<GRINS::ElasticCable<GRINS::HookesLaw> >;

