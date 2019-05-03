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

    //! Need to advance U_{n-1} before the other NumericVectors are updated
    virtual void preadvance_timestep( MultiphysicsSystem & system ) override;

  protected:

    //! FE variables for the flow
    VelocityVariable & _flow_vars;

    //! FE variables for the solid
    DisplacementVariable & _disp_vars;

    //! FE variables for the lagrange multiplier
    MultcomponentVectorVariable & _lambda_var;

    //! FE variable for the fluid pressure
    PressureFEVariable & _fluid_press_var;

    //! FE variable for the solid pressure
    PressureFEVariable & _solid_press_var;

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
    void reinit_single_ghosted_vector
    ( MultiphysicsSystem & system,
      const libMesh::NumericVector<libMesh::Number> & parallel_vector,
      std::unique_ptr<libMesh::NumericVector<libMesh::Number>> & ghosted_vector ) const;

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

  };

} // end namespace GRINS

#endif // GRINS_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_ABSTRACT_H
