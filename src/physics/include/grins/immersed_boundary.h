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
#include "grins/common.h"
#include "grins/solid_mechanics_abstract.h"
#include "grins/overlapping_fluid_solid_map.h"
#include "grins/immersed_boundary_augmented_sparsity.h"

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
    virtual void preassembly( MultiphysicsSystem & system );

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

    std::unique_ptr<ImmersedBoundaryAugmentedSparsity> _ibm_sparsity;

    std::unique_ptr<libMesh::FEMContext> _fluid_context;

    std::unique_ptr<libMesh::FEMContext> point_fluid_context;


    void prepare_fluid_context_batch( const MultiphysicsSystem & system,
                                      libMesh::dof_id_type fluid_elem_id,
                                      const AssemblyContext & solid_context,
                                      const std::vector<unsigned int> & solid_qpoint_indices,
                                      const std::vector<libMesh::Point> & solid_qpoints,
                                      std::vector<libMesh::Point> & solid_qpoints_subset,
                                      libMesh::FEMContext & fluid_context );

    void prepare_fluid_context( const MultiphysicsSystem & system,
                                const AssemblyContext & solid_context,
				const std::vector<libMesh::Point> & solid_qpoints,
                                unsigned int sqp, /* solid quadrature point */
                                libMesh::dof_id_type fluid_elem_id,
                                libMesh::FEMContext & fluid_context );

    void compute_residuals( AssemblyContext & solid_context,
                            libMesh::FEMContext & fluid_context,
                            unsigned int sqp,
                            libMesh::DenseSubVector<libMesh::Number> & Fuf,
                            libMesh::DenseSubVector<libMesh::Number> & Fvf,
                            libMesh::DenseSubVector<libMesh::Number> & Fus,
                            libMesh::DenseSubVector<libMesh::Number> & Fvs,
                            libMesh::DenseSubVector<libMesh::Number> & Fulm,
                            libMesh::DenseSubVector<libMesh::Number> & Fvlm);

    void fluid_residual_contribution( bool compute_jacobian, MultiphysicsSystem & system,
				      libMesh::FEMContext & fluid_context,
				      libMesh::dof_id_type fluid_elem_id,
				      AssemblyContext & solid_context,
				      const std::vector<libMesh::Point> & solid_qpoints,
				      unsigned int sqp,
				      libMesh::Real & jac,libMesh::Real delta,
				      libMesh::DenseSubVector<libMesh::Number> & Fuf,
				      libMesh::DenseSubVector<libMesh::Number> & Fvf,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kuf_us,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vs,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvf_us,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vlm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvf_ulm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm);

    void solid_residual_contribution( bool compute_jacobian,
				      AssemblyContext & solid_context,unsigned int sqp,
				      libMesh::Real & jac,libMesh::Real delta,
				      libMesh::DenseSubVector<libMesh::Number> & Fus,
				      libMesh::DenseSubVector<libMesh::Number> & Fvs,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kus_us,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_us,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kus_vs,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vs,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kus_ulm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_ulm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kus_vlm,
				      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vlm);

    void lambda_residual_contribution( bool compute_jacobian, MultiphysicsSystem & system,
				       libMesh::FEMContext & fluid_context,libMesh::dof_id_type fluid_elem_id,
				       AssemblyContext & solid_context,
				       const std::vector<libMesh::Point> & solid_qpoints,unsigned int sqp,
				       libMesh::Real & jac,libMesh::Real delta,
				       libMesh::DenseSubVector<libMesh::Number> & Fulm,
				       libMesh::DenseSubVector<libMesh::Number> & Fvlm,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_uf,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vf,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kulm_us,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_us,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vs,
				       libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vs);

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
