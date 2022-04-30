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

#ifndef GRINS_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_ABSTRACT_H
#define GRINS_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_ABSTRACT_H

//GRINS
#include "grins/physics.h"
#include "grins/single_variable.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/multiphysics_sys.h"
#include "grins/overlapping_fluid_solid_map.h"
#include "grins/overlapping_fluid_solid_coupling_functor.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/numeric_vector.h"

namespace GRINS
{
  class FictitiousDomainFluidStructureInteractionAbstract : public Physics
  {
  public:

    FictitiousDomainFluidStructureInteractionAbstract( const PhysicsName & physics_name,
                                                       const GetPot & input );

    FictitiousDomainFluidStructureInteractionAbstract() = delete;

    virtual ~FictitiousDomainFluidStructureInteractionAbstract() = default;

    //! Set time-evolving variables as time-evovling
    virtual void set_time_evolving_vars( libMesh::FEMSystem * system ) override;

    //! Init point locator, U_{n-1}
    virtual void auxiliary_init( MultiphysicsSystem & system ) override;

    virtual void reinit( MultiphysicsSystem & system ) override
    { this->reinit_point_locator(system); }

    //! Need to advance U_{n-1} before the other NumericVectors are updated
    virtual void preadvance_timestep( MultiphysicsSystem & system ) override;

  protected:

    //! FE variables for the flow
    VelocityVariable & _flow_vars;

    //! FE variables for the solid
    DisplacementVariable & _disp_vars;

    //! FE variables for the lagrange multiplier
    MultcomponentVectorVariable & _lambda_var;

    //! FE variable for the solid pressure
    PressureFEVariable & _solid_press_var;

    //! FE variable for the solid pressure
    PressureFEVariable & _fluid_press_var;

    //! The subdomain ids for the solid that is read from input
    std::set<libMesh::subdomain_id_type> _solid_subdomain_ids;

    //! The subdomain ids for the fluid that are read from input
    std::set<libMesh::subdomain_id_type> _fluid_subdomain_ids;

    //! Fluid density
    libMesh::Real _rho_fluid;

    //! Solid density
    libMesh::Real _rho_solid;

    //! Point locator on mesh for determining overlap
    std::unique_ptr<libMesh::PointLocatorBase> _point_locator;

    //! Data structure that holds current fluid/solid overlap
    std::unique_ptr<OverlappingFluidSolidMap> _fluid_solid_overlap;

    //! Coupling functor object that tells libMesh about additional coupling between fluid/solid
    std::unique_ptr<OverlappingFluidSolidCouplingFunctor> _coupling_functor;

    //! To contain the coupling between the variables
    std::unique_ptr<libMesh::CouplingMatrix> _coupling_matrix;

    std::unique_ptr<AssemblyContext> _fluid_context;

    //! Ghosted vector to store the previous time step solution (U_{n-1})
    std::unique_ptr<libMesh::NumericVector<libMesh::Number>> _prev_time_step_local_nonlinear_solution;

    void parse_subdomain_ids( const PhysicsName & physics_name,
                              const GetPot & input,
                              const std::string & subsection,
                              std::set<libMesh::subdomain_id_type> & subdomain_ids );

    bool is_solid_elem( libMesh::subdomain_id_type elem_id ) const
    { return _solid_subdomain_ids.find(elem_id) != _solid_subdomain_ids.end(); }

    bool is_fluid_elem( libMesh::subdomain_id_type elem_id ) const
    { return _fluid_subdomain_ids.find(elem_id) != _fluid_subdomain_ids.end(); }

    libMesh::Real parse_density( const PhysicsName & physics_name,
                                 const GetPot & input,
                                 const std::string & subsection );

    void reinit_point_locator( MultiphysicsSystem & system )
    {
      _point_locator.reset();
      _point_locator = system.get_mesh().sub_point_locator();
    }

    //! We need a parallel vector to store the U_{n-1} time step.
    void add_previous_time_step_parallel_vector_to_system( MultiphysicsSystem & system ) const;

    //! Rebuild and localize the ghosted vector from the parallel vector
    template<typename T>
    void reinit_single_ghosted_vector
    ( MultiphysicsSystem & system,
      const libMesh::NumericVector<libMesh::Number> & parallel_vector,
      T & ghosted_vector ) const
    {
      const libMesh::DofMap & dof_map = system.get_dof_map();

      ghosted_vector = libMesh::NumericVector<libMesh::Number>::build(system.comm());
      ghosted_vector->init(system.n_dofs(), system.n_local_dofs(),
                           dof_map.get_send_list(), false,
                           libMesh::GHOSTED);

      parallel_vector.localize( *ghosted_vector,dof_map.get_send_list());
    }

    //! Reinit all relevant ghosted vectors
    void reinit_all_ghosted_vectors( MultiphysicsSystem & system );

    //! Reinitializes all overlapping data structures and dependencies
    /*!
     * Rebuilds the OverlappingFluidSolidMap, OverlappingFluidSolidCouplingFunctor
     * and then reinitializes the appropriates parts of the DofMap and System objects.
     * If use_old_solution is true, this the OverlappingFluidSolidMap will use the solution
     * at the current time step (U_n) instead of the "next" (the one currently being sovle for)
     * time step (U_{n+1}). This will enable, essentially a semi-implicit method instead
     * of fully implicit method.
     */
    void reinit_overlapping_data( MultiphysicsSystem & system,
                                  bool use_old_solution );


    //! Build our local fluid context once and the reinit it when needed
    void build_fluid_context( MultiphysicsSystem & system );

    template<unsigned int Dim>
    void reinit_fluid_context( const libMesh::Point & x_qp,
                               const libMesh::Gradient & U,
                               const libMesh::Elem * fluid_elem,
                               AssemblyContext & fluid_context );

    //! Appropriately size, reposition coupled Jacobians for 2D case
    /*!
     *  Since the solid context wont have these Jacobians coupling fluid variables
     *  with any of the solid variables, we need to manually manage them. Two D
     *  version.
     */
    void prepare_jacobians(unsigned int n_fluid_dofs,
                           unsigned int n_solid_dofs,
                           unsigned int n_lambda_dofs,
                           unsigned int n_fluid_press_dofs,
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
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm,
                           libMesh::DenseMatrix<libMesh::Number> & Ks_pf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kus_pf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvs_pf ) const;

    //! Appropriately size, reposition coupled Jacobians for 3D case
    /*!
     *  Since the solid context wont have these Jacobians coupling fluid variables
     *  with any of the solid variables, we need to manually manage them. We reuse
     *  the 2D function for the 2D terms and then just handle the remaining
     *  3D terms.
     */
    void prepare_jacobians(unsigned int n_fluid_dofs,
                           unsigned int n_solid_dofs,
                           unsigned int n_lambda_dofs,
                           unsigned int n_fluid_press_dofs,
                           libMesh::DenseMatrix<libMesh::Number> & Kf_s,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_us,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vs,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ws,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_us,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_ws,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_us,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_vs,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_ws,
                           libMesh::DenseMatrix<libMesh::Number> & Klm_f,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kulm_wf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_uf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_wf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwlm_uf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwlm_vf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwlm_wf,
                           libMesh::DenseMatrix<libMesh::Number> & Kf_lm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_vlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kuf_wlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_ulm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvf_wlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_ulm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_vlm,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kwf_wlm,
                           libMesh::DenseMatrix<libMesh::Number> & Ks_pf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kus_pf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kvs_pf,
                           libMesh::DenseSubMatrix<libMesh::Number> & Kws_pf) const;

    void assemble_coupled_residuals( MultiphysicsSystem & system, AssemblyContext & fluid_context );

    //! Assemble local coupling terms into the global residual and Jacobian
    /*!
     *  Only handle the solid-fluid coupling terms. All the terms native/within
     *  the solid context will still be handled by FEMSystem.
     */
    template<unsigned int Dim>
    void assemble_coupled_jacobians( MultiphysicsSystem & system,
                                     const AssemblyContext & solid_context,
                                     AssemblyContext & fluid_context,
                                     unsigned int n_fluid_dofs,
                                     unsigned int n_solid_dofs,
                                     unsigned int n_lambda_dofs,
                                     unsigned int n_fluid_press_dofs,
                                     libMesh::DenseMatrix<libMesh::Number> & Kf_s,
                                     libMesh::DenseMatrix<libMesh::Number> & Klm_f,
                                     libMesh::DenseMatrix<libMesh::Number> & Kf_lm,
                                     libMesh::DenseMatrix<libMesh::Number> & Ks_pf );

    void get_prev_time_elem_solution( AssemblyContext & solid_context,
                                      libMesh::DenseVector<libMesh::Number> & prev_time_solution ) const;

    template<unsigned int Dim>
    void compute_displacement_accel( AssemblyContext & solid_context,
                                     const unsigned int qp,
                                     libMesh::Gradient & Uddot /*\ddot{U}*/ );

    template<unsigned int Dim, bool UseOldDisplacement>
    void compute_displaced_point( const MultiphysicsSystem & system,
                                  AssemblyContext & solid_context,
                                  const unsigned int qp,
                                  libMesh::Gradient & U) const;

  };

} // end namespace GRINS

#endif // GRINS_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_ABSTRACT_H
