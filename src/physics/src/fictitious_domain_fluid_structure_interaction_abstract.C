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

// This class
#include "grins/fictitious_domain_fluid_structure_interaction_abstract.h"

// GRINS
#include "grins/variable_warehouse.h"
#include "grins/materials_parsing.h"
#include "grins/generic_ic_handler.h"


// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/unsteady_solver.h"
#include "libmesh/sparse_matrix.h"

namespace GRINS
{
  FictitiousDomainFluidStructureInteractionAbstract::
  FictitiousDomainFluidStructureInteractionAbstract( const PhysicsName & physics_name,
                                                     const GetPot & input )
    :Physics(physics_name, input),
     _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
     _disp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<DisplacementVariable>(VariablesParsing::disp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
     _lambda_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<MultcomponentVectorVariable>(VariablesParsing::vector_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
     _fluid_press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
     _solid_press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    _lambda_var.set_is_constraint_var(true);
    _solid_press_var.set_is_constraint_var(true);

    // Parse fluid data
    {
      std::string subsection("Fluid");
      this->parse_subdomain_ids(physics_name, input, subsection, _fluid_subdomain_ids );
      _rho_fluid = this->parse_density(physics_name, input, subsection);
    }

    // Parse solid data
    {
      std::string subsection("Solid");
      this->parse_subdomain_ids(physics_name, input, subsection, _solid_subdomain_ids );
      _rho_solid = this->parse_density(physics_name, input, subsection);
    }

    this->_ic_handler = new GenericICHandler(physics_name, input);
  }

  void FictitiousDomainFluidStructureInteractionAbstract::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Velocity are first order in time
    system->time_evolving(_flow_vars.u(),1);
    system->time_evolving(_flow_vars.v(),1);

    if ( _flow_vars.dim() == 3 )
      system->time_evolving(_flow_vars.w(),1);


    // In this formulation, we treat solid displacements
    // as first order in time.
    system->time_evolving(_disp_vars.u(),1);

    if( _disp_vars.dim() >= 2 )
      system->time_evolving(_disp_vars.v(),1);

    if ( _disp_vars.dim() == 3 )
      system->time_evolving(_disp_vars.w(),1);
  }

  void FictitiousDomainFluidStructureInteractionAbstract::auxiliary_init( MultiphysicsSystem & system )
  {
    this->reinit_point_locator(system);
    this->add_previous_time_step_parallel_vector_to_system(system);
    this->build_fluid_context(system);
  }

  void FictitiousDomainFluidStructureInteractionAbstract::preadvance_timestep( MultiphysicsSystem & system )
  {
    libMesh::NumericVector<libMesh::Number> & prev_time_step_nonlinear_soln =
      system.get_vector("_prev_time_step_nonlinear_solution");

    libMesh::NumericVector<libMesh::Number> & old_nonlinear_soln =
      system.get_vector("_old_nonlinear_solution");

    // timestep n copied to timestep n-1
    prev_time_step_nonlinear_soln = old_nonlinear_soln;

    prev_time_step_nonlinear_soln.localize( *_prev_time_step_local_nonlinear_solution,
                                            system.get_dof_map().get_send_list() );
  }


  void FictitiousDomainFluidStructureInteractionAbstract::parse_subdomain_ids
  ( const PhysicsName & physics_name,
    const GetPot & input,
    const std::string & subsection,
    std::set<libMesh::subdomain_id_type> & subdomain_ids)
  {
    std::string var_string = "Physics/"+physics_name+"/"+subsection+"/subdomains";

    if( !input.have_variable(var_string) )
      libmesh_error_msg("ERROR: Could not find required input variable: "+var_string+"!\n");

    unsigned int n_subdomains = input.vector_variable_size(var_string);
    if( n_subdomains == 0 )
      libmesh_error_msg("Error: Must specify at least one enabled subdomain for "+var_string+"!\n");

    for( unsigned int s = 0; s < n_subdomains; s++ )
      subdomain_ids.insert( input(var_string,std::numeric_limits<libMesh::subdomain_id_type>::max(), s) );
  }

  libMesh::Real FictitiousDomainFluidStructureInteractionAbstract::parse_density
  ( const PhysicsName & physics_name,
    const GetPot & input,
    const std::string & subsection )
  {
    std::string mat_string("Physics/"+physics_name+"/"+subsection+"/material");
    if( !input.have_variable(mat_string) )
      libmesh_error_msg("ERROR: Could not find required input variable: "+mat_string+"!\n");

    std::string material = input(mat_string,"DIE!");

    std::string option("Materials/"+material+"/Density/value");
    MaterialsParsing::check_for_input_option(input,option);
    return input(option,0.0);
  }

  void FictitiousDomainFluidStructureInteractionAbstract::add_previous_time_step_parallel_vector_to_system
  ( MultiphysicsSystem & system ) const
  {
    libMesh::NumericVector<libMesh::Number> & prev_time_step_nonlinear_soln =
      system.add_vector("_prev_time_step_nonlinear_solution");

    // This should be a parallel vector. The ghosted version we will store locally in this class.
    // We don't deal with localizing the ghosted version now because that will be taken care of
    // with all other ghosted vectors.
    if(!prev_time_step_nonlinear_soln.initialized())
      prev_time_step_nonlinear_soln.init(system.n_dofs(), system.n_local_dofs(), false, libMesh::PARALLEL);
  }

  void FictitiousDomainFluidStructureInteractionAbstract::reinit_single_ghosted_vector
    ( MultiphysicsSystem & system,
      const libMesh::NumericVector<libMesh::Number> & parallel_vector,
      std::unique_ptr<libMesh::NumericVector<libMesh::Number>> & ghosted_vector ) const
  {
    const libMesh::DofMap & dof_map = system.get_dof_map();

    ghosted_vector = libMesh::NumericVector<libMesh::Number>::build(system.comm());
    ghosted_vector->init(system.n_dofs(), system.n_local_dofs(),
                         dof_map.get_send_list(), false,
                         libMesh::GHOSTED);

    parallel_vector.localize( *ghosted_vector,dof_map.get_send_list());
  }

  void FictitiousDomainFluidStructureInteractionAbstract::reinit_all_ghosted_vectors( MultiphysicsSystem & system )
  {

    // Reinit current local solution, i.e. U_{n+1}
    this->reinit_single_ghosted_vector(system,
                                       (*system.solution),
                                       system.current_local_solution);

    // Also need to update unsteady system vectors
    libMesh::TimeSolver & time_solver_raw = system.get_time_solver();
    libMesh::UnsteadySolver * unsteady_solver = dynamic_cast<libMesh::UnsteadySolver *>(&time_solver_raw);
    if( unsteady_solver)
      {
        // Old nonlinear solution, i.e. U_n
        {
          const libMesh::NumericVector<libMesh::Number> & old_nonlinear_soln =
            system.get_vector("_old_nonlinear_solution");

          this->reinit_single_ghosted_vector(system,
                                             old_nonlinear_soln,
                                             unsteady_solver->old_local_nonlinear_solution);
        }

        // prev_time_step nonlinear solution, i.e. U_{n-1}
        {
          const libMesh::NumericVector<libMesh::Number> & prev_time_step_nonlinear_soln =
            system.get_vector("_prev_time_step_nonlinear_solution");

          this->reinit_single_ghosted_vector(system,
                                             prev_time_step_nonlinear_soln,
                                             _prev_time_step_local_nonlinear_solution);
        }

      } // if(unsteady_solver)
  }

  void FictitiousDomainFluidStructureInteractionAbstract::reinit_overlapping_data( MultiphysicsSystem & system,
                                                                                   bool use_old_solution )
  {
    // We need to rebuild the overlap map each time because the position
    // of the solid can change in between each Newton step.
    _fluid_solid_overlap.reset( new OverlappingFluidSolidMap(system,
                                                             (*_point_locator),
                                                             _solid_subdomain_ids,
                                                             _fluid_subdomain_ids,
                                                             _disp_vars,
                                                             use_old_solution) );

    libMesh::DofMap & dof_map = system.get_dof_map();

    if(!_coupling_functor)
      {
        _coupling_functor.reset( new OverlappingFluidSolidCouplingFunctor(system.get_mesh(),
                                                                          *_coupling_matrix,
                                                                          *_fluid_solid_overlap) );
        dof_map.add_coupling_functor(*_coupling_functor);
      }
    else
      {
        dof_map.remove_coupling_functor(*_coupling_functor);
        _coupling_functor.reset( new OverlappingFluidSolidCouplingFunctor(system.get_mesh(),
                                                                          *_coupling_matrix,
                                                                          *_fluid_solid_overlap) );
        dof_map.add_coupling_functor(*_coupling_functor);
      }

    // Handle algebraic coupling
    // We have to manually make these calls since we couldn't attach this
    // functor before the mesh and dof_map was prepared.
    dof_map.reinit_send_list(system.get_mesh());

    this->reinit_all_ghosted_vectors( system );

    // Handle full coupling (sparsity pattern)
    dof_map.clear_sparsity();
    dof_map.compute_sparsity(system.get_mesh());

    // Now to reinit the matrix since we changed the sparsity pattern
    libMesh::SparseMatrix<libMesh::Number> & matrix = system.get_matrix("System Matrix");
    libmesh_assert(dof_map.is_attached(matrix));
    matrix.init();
  }

  void FictitiousDomainFluidStructureInteractionAbstract::build_fluid_context( MultiphysicsSystem & system )
  {
    std::unique_ptr<libMesh::DiffContext> raw_context = system.build_context();
    AssemblyContext * context = libMesh::cast_ptr<AssemblyContext*>(raw_context.release());
    _fluid_context.reset(context);
  }

} // end namespace GRINS
