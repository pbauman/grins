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

#ifndef GRINS_ELASTIC_CABLE_H
#define GRINS_ELASTIC_CABLE_H


//GRINS
#include "grins/elastic_cable_base.h"
#include "grins/elasticity_tensor.h"


namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticCable : public ElasticCableBase<StressStrainLaw>
  {
  public:

    ElasticCable( const PhysicsName& physics_name, const GetPot& input,
                  bool is_compressible );

    virtual ~ElasticCable(){};

    //! Register postprocessing variables for ElasticCable
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context,
                                          CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext& context,
                                CachedValues& /*cache*/ )
    { this->mass_residual_impl(compute_jacobian,
                               context,
                               &libMesh::FEMContext::interior_accel,
                               &libMesh::DiffContext::get_elem_solution_accel_derivative); }

    //! Compute the registered postprocessed quantities
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

    //! Precompute data needed for residual inline function
    void get_grad_uvw(const AssemblyContext& context, unsigned int qp,
                                   libMesh::Gradient &gradu,
                                   libMesh::Gradient &gradv,
                                   libMesh::Gradient &gradw);


    //! Precompute tau, needed for residual
    void get_stress_and_elasticity(const AssemblyContext& context, unsigned int qp,
                               const libMesh::Gradient &gradu,
                               const libMesh::Gradient &gradv,
                               const libMesh::Gradient &gradw,
                               libMesh::TensorValue<libMesh::Real> & t, /*stress (tau)*/
                               ElasticityTensor & C );
   
  private:

    ElasticCable();

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _stress_indices;

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _strain_indices;

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _force_indices;
    
  }; //end class ElasticCable


  /* ------------- Inline Functions ---------------*/

  template<typename StressStrainLaw> inline
  void ElasticCable<StressStrainLaw>::get_grad_uvw(const AssemblyContext& context, unsigned int qp,
                                                         libMesh::Gradient &gradu,
                                                         libMesh::Gradient &gradv,
                                                         libMesh::Gradient &gradw)
  { 
    const unsigned int n_u_dofs = context.get_dof_indices(this->_disp_vars.u()).size();
    
    // All shape function gradients are w.r.t. master element coordinates
    const std::vector<std::vector<libMesh::Real> >& dphi_dxi = this->get_fe(context)->get_dphidxi();

    const libMesh::DenseSubVector<libMesh::Number>& u_coeffs = context.get_elem_solution( this->_disp_vars.u() );
    const libMesh::DenseSubVector<libMesh::Number>& v_coeffs = context.get_elem_solution( this->_disp_vars.v() );
    const libMesh::DenseSubVector<libMesh::Number>& w_coeffs = context.get_elem_solution( this->_disp_vars.w() );

    // Compute gradients  w.r.t. master element coordinates
    for( unsigned int d = 0; d < n_u_dofs; d++ )
      {
        libMesh::RealGradient u_gradphi( dphi_dxi[d][qp] );
        gradu += u_coeffs(d)*u_gradphi;
        gradv += v_coeffs(d)*u_gradphi;
        gradw += w_coeffs(d)*u_gradphi;
      }
  }

  
  template<typename StressStrainLaw> inline
  void ElasticCable<StressStrainLaw>::get_stress_and_elasticity(const AssemblyContext& context, unsigned int qp,
                                                     const libMesh::Gradient &gradu,
                                                     const libMesh::Gradient &gradv,
                                                     const libMesh::Gradient &gradw,
                                                     libMesh::TensorValue<libMesh::Real> & tau,
                                                     ElasticityTensor & C)
  {
    // Need these to build up the covariant and contravariant metric tensors
    const std::vector<libMesh::RealGradient>& dxdxi  = this->get_fe(context)->get_dxyzdxi();

    const unsigned int dim = 1; // The cable dimension is always 1 for this physics

    // Compute & store gradients  w.r.t. actual element coordinates
    libMesh::RealGradient grad_x( dxdxi[qp](0) );
    libMesh::RealGradient grad_y( dxdxi[qp](1) );
    libMesh::RealGradient grad_z( dxdxi[qp](2) );
    
    libMesh::TensorValue<libMesh::Real> a_cov, a_contra, A_cov, A_contra;
    libMesh::Real lambda_sq = 0;
    
    this->compute_metric_tensors( qp, *(this->get_fe(context)), context,
                                  gradu, gradv, gradw,
                                  a_cov, a_contra, A_cov, A_contra,
                                  lambda_sq );

    // Compute stress tensor
    this->_stress_strain_law.compute_stress_and_elasticity(dim,a_contra,a_cov,A_contra,A_cov,tau,C);
  }



} // end namespace GRINS


#endif /* GRINS_ELASTIC_CABLE_H */
