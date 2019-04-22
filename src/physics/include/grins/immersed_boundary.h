//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_IMMERSED_BOUNDARY_H
#define GRINS_IMMERSED_BOUNDARY_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/single_variable.h"
#include "grins/common.h"
#include "grins/solid_mechanics_abstract.h"
#include "grins/overlapping_fluid_solid_map.h"
#include "grins/immersed_boundary_coupling_functor.h"

//libMesh
#include "libmesh/fem_context.h"

namespace GRINS
{

  //! Physics class for the Immersed Boundary Method
  /*!
    This physics class implements the classical FEM Immersed  Boundary Method.
    This is a templated class, the class SolidMech can be instantiated as a specific type
    (right now as: ElasticCable or ElasticMembrane)
   */

  template<typename SolidMech>
  class ImmersedBoundary : public Physics
  {
  public:

    ImmersedBoundary( const std::string & my_physics_name,
                      std::unique_ptr<SolidMech> & solid_mech_ptr,
                      const GetPot& input );

    ImmersedBoundary() = delete;

    virtual ~ImmersedBoundary() = default;

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Init point locator and local FEMContext
    virtual void auxiliary_init( MultiphysicsSystem & system );

    //! Context initialization
    virtual void init_context( AssemblyContext& context );

    //! Residual contributions from the solid in the flow
    virtual void element_time_derivative( bool compute_jacobian, AssemblyContext& context );

    //! Cache mesh information needed for residual computation
    virtual void presolve( MultiphysicsSystem & system );

    //! We need to reinit the point locator when libMesh::System::reinit is called.
    virtual void reinit( MultiphysicsSystem & system );

    //! Override to point to solid Physics for displacement initial conditions
    /*! It would make sense for the user, under the current schema, to put the
        displacement initial conditions in the solid Physics part of the input
        file. But, since that Physics doesn't actually get added, only this
        ImmersedBoundary, then we need to internally point to the solids ics function.*/
    virtual void init_ics( libMesh::FEMSystem* system,
                           libMesh::CompositeFunction<libMesh::Number>& all_ics );

  private:

    //! FE variables for the flow
    VelocityVariable & _flow_vars;

    //! FE variables for the solid
    DisplacementVariable & _disp_vars;

    //! FE variables for the lagrange multiplier
    LagrangeMultVectorVariable & _lambda_var;

    PressureFEVariable & _press_var;

    PressureFEVariable & _solid_press_var;

    //! Solid Mechanics from the ibm factory
    std::unique_ptr<SolidMech> _solid_mech;

    //! The fluid mechanics associated with the IBM method from the input
    std::string _fluid_mechanics;

    //! The solid mechanics associated with the IBM method from the input
    std::string _solid_mechanics;

    //! The subdomain ids for the solid that is read from input
    std::set<libMesh::subdomain_id_type> _solid_subdomain_set;

    //! The subdomain ids for the fluid that are read from input
    std::set<libMesh::subdomain_id_type> _fluid_subdomain_set;

    std::unique_ptr<libMesh::PointLocatorBase> _point_locator;

    std::unique_ptr<OverlappingFluidSolidMap> _fluid_solid_overlap;

    std::unique_ptr<ImmersedBoundaryCouplingFunctor> _coupling_functor;

    std::unique_ptr<libMesh::CouplingMatrix> _coupling_matrix;

    std::unique_ptr<libMesh::FEMContext> _fluid_context;

    void setup_coupling_matrix( const VelocityVariable & flow_vars,
                                const DisplacementVariable & disp_vars,
                                const LagrangeMultVectorVariable & lambda_var,
                                const PressureFEVariable & press_var,
                                libMesh::CouplingMatrix & coupling_matrix );

    void diagonally_coupled_vars( const MultcomponentVectorVariable & var1,
                                  const MultcomponentVectorVariable & var2,
                                  libMesh::CouplingMatrix & coupling_matrix );

    void fully_coupled_vars( const MultcomponentVectorVariable & var1,
                             const MultcomponentVectorVariable & var2,
                             libMesh::CouplingMatrix & coupling_matrix );

    void compute_residuals( const AssemblyContext & solid_context,
                            const libMesh::FEMContext & fluid_context,
                            unsigned int sqp,
                            libMesh::DenseSubVector<libMesh::Number> & Fuf,
                            libMesh::DenseSubVector<libMesh::Number> & Fvf,
                            libMesh::DenseSubVector<libMesh::Number> & Fus,
                            libMesh::DenseSubVector<libMesh::Number> & Fvs,
                            libMesh::DenseSubVector<libMesh::Number> & Fulm,
                            libMesh::DenseSubVector<libMesh::Number> & Fvlm,
                            libMesh::DenseSubVector<libMesh::Number> & Fp);

    void compute_numerical_jacobians(const MultiphysicsSystem & system,
                                     AssemblyContext & solid_context,
                                     libMesh::FEMContext & fluid_context,
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
                                     libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs);

    void compute_analytic_jacobians(const MultiphysicsSystem & system,
                                    AssemblyContext & solid_context,
                                    libMesh::FEMContext & fluid_context,
                                    const std::vector<unsigned int> & quad_points,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
                                    libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf);

    void finite_difference_residuals(const MultiphysicsSystem & system,
                                     const std::vector<unsigned int> & quad_points,
                                     AssemblyContext & solid_context,
                                     libMesh::FEMContext & fluid_context,
                                     const libMesh::Real delta,
                                     libMesh::Number & coeff,
                                     libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                                     libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual);

    void compute_lambda_derivs(const MultiphysicsSystem & system,
                               const std::vector<unsigned int> & quad_points,
                               AssemblyContext & solid_context,
                               libMesh::FEMContext & fluid_context,
                               const libMesh::Real delta,
                               libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                               libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                               libMesh::DenseSubVector<libMesh::Number> & lambda_coeff,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kuf,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kvf,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kus,
                               libMesh::DenseSubMatrix<libMesh::Number> & Kvs);

    void compute_fluid_derivs(const MultiphysicsSystem & system,
                              const std::vector<unsigned int> & quad_points,
                              AssemblyContext & solid_context,
                              libMesh::FEMContext & fluid_context,
                              const libMesh::Real delta,
                              libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
                              libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
                              libMesh::DenseSubVector<libMesh::Number> & fluid_coeff,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kulm,
                              libMesh::DenseSubMatrix<libMesh::Number> & Kvlm);

    void compute_solid_derivs(const MultiphysicsSystem & system,
                              const std::vector<unsigned int> & quad_points,
                              AssemblyContext & solid_context,
                              libMesh::FEMContext & fluid_context,
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

    void compute_press_derivs(const MultiphysicsSystem & system,
			      const std::vector<unsigned int> & quad_points,
			      AssemblyContext & solid_context,
                              libMesh::FEMContext & fluid_context,
			      const libMesh::Real delta,
			      libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
			      libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
			      libMesh::DenseSubVector<libMesh::Number> & press_coeff,
			      libMesh::DenseSubMatrix<libMesh::Number> & Kus,
			      libMesh::DenseSubMatrix<libMesh::Number> & Kvs);

    void prepare_jacobians(unsigned int n_fluid_dofs,
                           unsigned int n_solid_dofs,
                           unsigned int n_lambda_dofs,
                           libMesh::DenseMatrix<libMesh::Number> & Kf_s,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_us,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vs,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_us,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs,
                           libMesh::DenseMatrix<libMesh::Number> & Klm_f,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_uf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf,
                           libMesh::DenseMatrix<libMesh::Number> & Kf_lm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_ulm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm) const;

    void assemble_fluid_jacobians(MultiphysicsSystem & system,
                                  const AssemblyContext & solid_context,
                                  const libMesh::FEMContext & fluid_context,
                                  unsigned int n_fluid_dofs,
                                  unsigned int n_solid_dofs,
                                  unsigned int n_lambda_dofs,
                                  libMesh::DenseMatrix<libMesh::Number> & Kf_s,
                                  libMesh::DenseMatrix<libMesh::Number> & Klm_f,
                                  libMesh::DenseMatrix<libMesh::Number> & Kf_lm) const;

    bool is_solid_elem( libMesh::subdomain_id_type elem_id );

    bool is_fluid_elem( libMesh::subdomain_id_type elem_id );

    void eval_first_Piola( const libMesh::Gradient & grad_u,
                      const libMesh::Gradient & grad_v,
                      libMesh::TensorValue<libMesh::Real> & F );

    void eval_deform_gradient( const libMesh::Gradient & grad_u,
			       const libMesh::Gradient & grad_v,
			       libMesh::TensorValue<libMesh::Real> & F );

    void compute_prev_timestep_deform_gradient( const AssemblyContext & solid_context_const,
                                                const unsigned int qp,
                                                libMesh::TensorValue<libMesh::Real> & Fold );

    void eval_deform_grad_rate( const libMesh::Gradient & grad_udot,
				const libMesh::Gradient & grad_vdot,
				libMesh::TensorValue<libMesh::Real> & Fdot );

    libMesh::Point compute_displaced_point( const MultiphysicsSystem & system,
                                            AssemblyContext & solid_context,
                                            const libMesh::Point & x_qp,
                                            const unsigned int qp ) const;

    const libMesh::Elem * get_fluid_elem( const MultiphysicsSystem & system,
                                          AssemblyContext & solid_context,
                                          const libMesh::Point & x_qp,
                                          const unsigned int qp ) const;

    void compute_ibm_residuals(const MultiphysicsSystem & system,
                               AssemblyContext & solid_context,
                               libMesh::FEMContext & fluid_context,
                               const std::vector<unsigned int> & quad_points);

    void reinit_ghosted_vectors( MultiphysicsSystem & system );

    void print_coupling_matrix( const libMesh::CouplingMatrix & coupling_matrix, const unsigned int n );
  };

  template<typename SolidMech>
  inline
  bool ImmersedBoundary<SolidMech>::is_solid_elem( libMesh::subdomain_id_type elem_id )
  {
    return _solid_subdomain_set.find(elem_id) != _solid_subdomain_set.end();
  }

  template<typename SolidMech>
  inline
  bool ImmersedBoundary<SolidMech>::is_fluid_elem( libMesh::subdomain_id_type elem_id )
  {
    return _fluid_subdomain_set.find(elem_id) != _fluid_subdomain_set.end();
  }

} //End namespace block
#endif //GRINS_IMMERSED_BOUNDARY_H
