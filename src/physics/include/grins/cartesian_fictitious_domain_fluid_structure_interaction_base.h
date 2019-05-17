//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H
#define GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H

// GRINS
#include "grins/fictitious_domain_fluid_structure_interaction_abstract.h"
#include "grins/hyperelastic_strain_energy.h"
#include "grins/mooney_rivlin.h"
#include "grins/cartesian_hyperelasticity.h"
#include "grins/incompressible_hyperelasticity_weak_form.h"

namespace GRINS
{
  template<unsigned int Dim, bool UseOldDisplacement>
  class CartesianFictitiousDomainFluidStructureInteractionBase :
    public FictitiousDomainFluidStructureInteractionAbstract
  {
  public:

    CartesianFictitiousDomainFluidStructureInteractionBase( const PhysicsName & physics_name,
                                                            const GetPot & input );

    CartesianFictitiousDomainFluidStructureInteractionBase() = delete;

    virtual ~CartesianFictitiousDomainFluidStructureInteractionBase() = default;

    //! Context initialization
    virtual void init_context( AssemblyContext & context ) override;

    virtual void element_time_derivative( bool compute_jacobian, AssemblyContext& context ) override;

  protected:

    std::unique_ptr<HyperelasticStrainEnergy<MooneyRivlin>> _strain_energy;

    void check_variable_dim_consistency() const;

    std::string dim_error_msg(unsigned int var_dim) const;

    //! 2D deformation gradient for plane strain
    /* F(2,2) = 1 will be consistent with plane strain for everything that consumes F */
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v ) const;

    //! 3D deformation gradient
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v,
                                       const libMesh::Gradient & grad_w ) const;

    libMesh::Tensor form_fluid_def_gradient(AssemblyContext & solid_context,
                                            const libMesh::Tensor & F,
                                            const unsigned int qp) const;



    void compute_ibm_residuals(MultiphysicsSystem & system,
                               AssemblyContext & solid_context,
                               AssemblyContext & fluid_context,
                               const std::vector<unsigned int> & qp_indices);


    void compute_residuals(AssemblyContext & solid_context,const AssemblyContext & fluid_context,unsigned int qp,
                           IncompressibleHyperelasticityWeakForm<MooneyRivlin> & weak_form,
                           libMesh::DenseSubVector<libMesh::Number> & Fuf,
                           libMesh::DenseSubVector<libMesh::Number> & Fvf,
                           libMesh::DenseSubVector<libMesh::Number> & Fus,
                           libMesh::DenseSubVector<libMesh::Number> & Fvs,
                           libMesh::DenseSubVector<libMesh::Number> & Fulm,
                           libMesh::DenseSubVector<libMesh::Number> & Fvlm,
                           libMesh::DenseSubVector<libMesh::Number> & Fps);

    void compute_numerical_jacobians(MultiphysicsSystem & system,
                                     AssemblyContext & solid_context,
                                     AssemblyContext & fluid_context,
                                     const std::vector<unsigned int> & quad_points,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vlm,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvf_ulm,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vf,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_uf,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kuf_us,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vs,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvf_us,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kus_pf,
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvs_pf);

    void finite_difference_residuals(MultiphysicsSystem & system,
                                     const std::vector<unsigned int> & quad_points,
                                     AssemblyContext & solid_context,
                                     AssemblyContext & fluid_context,
                                     const libMesh::Real delta,
                                     libMesh::Number & coeff,
                                     libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                                     libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual);

    void compute_lambda_derivs(MultiphysicsSystem & system,
                               const std::vector<unsigned int> & quad_points,
                               AssemblyContext & solid_context,
                               AssemblyContext & fluid_context,
                               const libMesh::Real delta,
                               libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                               libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                               libMesh::DenseSubVector<libMesh::Number> & lambda_coeff,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kuf,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kvf,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kus,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kvs);

    void compute_fluid_derivs(MultiphysicsSystem & system,
                              const std::vector<unsigned int> & quad_points,
                              AssemblyContext & solid_context,
                              AssemblyContext & fluid_context,
                              const libMesh::Real delta,
                              libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                              libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                              libMesh::DenseSubVector<libMesh::Number> & fluid_coeff,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kulm,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kvlm);

    void compute_solid_derivs(MultiphysicsSystem & system,
                              const std::vector<unsigned int> & quad_points,
                              AssemblyContext & solid_context,
                              AssemblyContext & fluid_context,
                              const libMesh::Real delta,
                              libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                              libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                              libMesh::DenseSubVector<libMesh::Number> & solid_coeff,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kuf,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kvf,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kus,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kvs,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kulm,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kvlm,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kp);

    void compute_press_derivs(MultiphysicsSystem & system,
			      const std::vector<unsigned int> & quad_points,
			      AssemblyContext & solid_context,
                              AssemblyContext & fluid_context,
			      const libMesh::Real delta,
			      libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
			      libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
			      libMesh::DenseSubVector<libMesh::Number> & press_coeff,
			      libMesh::DenseSubMatrix<libMesh::Number> & Kus,
			      libMesh::DenseSubMatrix<libMesh::Number> & Kvs);

    void compute_fluid_press_derivs(MultiphysicsSystem & system,
                                    const std::vector<unsigned int> & quad_points,
                                    AssemblyContext & solid_context,
                                    AssemblyContext & fluid_context,
                                    const libMesh::Real delta,
                                    libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                                    libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                                    libMesh::DenseSubVector<libMesh::Number> & press_coeff,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kus,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kvs);

  };

} // end namespace GRINS

#endif // GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H
