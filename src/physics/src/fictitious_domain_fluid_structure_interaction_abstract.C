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

  void FictitiousDomainFluidStructureInteractionAbstract::prepare_jacobians
  (unsigned int n_fluid_dofs, unsigned int n_solid_dofs, unsigned int n_lambda_dofs,
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
   libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm) const
  {
    libmesh_assert_equal_to( this->_flow_vars.dim(),this->_disp_vars.dim() );
    libmesh_assert_equal_to( this->_lambda_var.dim(),this->_flow_vars.dim() );

    Kf_s.resize( this->_flow_vars.dim()*n_fluid_dofs, this->_disp_vars.dim()*n_solid_dofs );

    // We need to manually manage the indexing since we're working only on this particular subblock
    Kuf_us.reposition( 0, 0, n_fluid_dofs, n_solid_dofs );
    Kuf_vs.reposition( 0, n_solid_dofs, n_fluid_dofs, n_solid_dofs );
    Kvf_us.reposition( n_fluid_dofs, 0, n_fluid_dofs, n_solid_dofs );
    Kvf_vs.reposition( n_fluid_dofs, n_solid_dofs, n_fluid_dofs, n_solid_dofs );

    Klm_f.resize( this->_lambda_var.dim()*n_lambda_dofs, this->_flow_vars.dim()*n_fluid_dofs );

    Kulm_uf.reposition( 0, 0, n_lambda_dofs, n_fluid_dofs );
    Kulm_vf.reposition( 0, n_fluid_dofs, n_lambda_dofs, n_fluid_dofs );
    Kvlm_uf.reposition( n_lambda_dofs, 0, n_lambda_dofs, n_fluid_dofs );
    Kvlm_vf.reposition( n_lambda_dofs, n_fluid_dofs, n_lambda_dofs, n_fluid_dofs );

    Kf_lm.resize( this->_flow_vars.dim()*n_fluid_dofs, this->_lambda_var.dim()*n_lambda_dofs );

    Kuf_ulm.reposition( 0, 0, n_fluid_dofs, n_lambda_dofs );
    Kuf_vlm.reposition( 0, n_lambda_dofs, n_fluid_dofs, n_lambda_dofs );
    Kvf_ulm.reposition( n_fluid_dofs, 0, n_fluid_dofs, n_lambda_dofs );
    Kvf_vlm.reposition( n_fluid_dofs, n_lambda_dofs, n_fluid_dofs, n_lambda_dofs );
  }

  void FictitiousDomainFluidStructureInteractionAbstract::prepare_jacobians
  (unsigned int n_fluid_dofs, unsigned int n_solid_dofs, unsigned int n_lambda_dofs,
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
   libMesh::DenseSubMatrix<libMesh::Number> & Kwf_wlm) const
  {
    // We can reuse all the 2-D bits and just handle the 3-D terms
    this->prepare_jacobians(n_fluid_dofs, n_solid_dofs, n_lambda_dofs,
                            Kf_s,
                            Kuf_us,Kuf_vs,Kvf_us,Kvf_vs,
                            Klm_f,
                            Kulm_uf,Kulm_vf, Kvlm_uf,Kvlm_vf,
                            Kf_lm,
                            Kuf_ulm,Kuf_vlm,Kvf_ulm,Kvf_vlm);

    // Fluid-Solid
    Kuf_ws.reposition( 0,              2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
    Kvf_ws.reposition( n_fluid_dofs,   2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
    Kwf_us.reposition( 2*n_fluid_dofs, 0,              n_fluid_dofs, n_solid_dofs );
    Kwf_vs.reposition( 2*n_fluid_dofs, n_solid_dofs,   n_fluid_dofs, n_solid_dofs );
    Kwf_ws.reposition( 2*n_fluid_dofs, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );

    // Lambda-fluid
    Kulm_wf.reposition( 0,               2*n_fluid_dofs, n_lambda_dofs, n_fluid_dofs );
    Kvlm_wf.reposition( n_lambda_dofs,   2*n_fluid_dofs,  n_lambda_dofs, n_fluid_dofs );
    Kwlm_uf.reposition( 2*n_lambda_dofs, 0,              n_lambda_dofs, n_fluid_dofs );
    Kwlm_vf.reposition( 2*n_lambda_dofs, n_fluid_dofs,   n_lambda_dofs, n_fluid_dofs );
    Kwlm_wf.reposition( 2*n_lambda_dofs, 2*n_fluid_dofs, n_lambda_dofs, n_fluid_dofs );

    // Fluid-lambda
    Kuf_wlm.reposition( 0,              2*n_lambda_dofs, n_fluid_dofs, n_lambda_dofs );
    Kvf_wlm.reposition( n_fluid_dofs,   2*n_lambda_dofs, n_fluid_dofs, n_lambda_dofs );
    Kwf_ulm.reposition( 2*n_fluid_dofs, 0,               n_fluid_dofs, n_lambda_dofs );
    Kwf_vlm.reposition( 2*n_fluid_dofs, n_lambda_dofs,   n_fluid_dofs, n_lambda_dofs );
    Kwf_wlm.reposition( 2*n_fluid_dofs, 2*n_lambda_dofs, n_fluid_dofs, n_lambda_dofs );
  }

  template<unsigned int Dim>
  void FictitiousDomainFluidStructureInteractionAbstract::assemble_coupled_terms
  ( bool compute_jacobian,
    MultiphysicsSystem & system,
    const AssemblyContext & solid_context,
    AssemblyContext & fluid_context,
    unsigned int n_fluid_dofs,
    unsigned int n_solid_dofs,
    unsigned int n_lambda_dofs,
    libMesh::DenseMatrix<libMesh::Number> & Kf_s,
    libMesh::DenseMatrix<libMesh::Number> & Klm_f,
    libMesh::DenseMatrix<libMesh::Number> & Kf_lm )
  {
    // Now that we've looped over all the quadrature points on this fluid element
    // we can now assemble the fluid residual (the residuals in the solid context
    // automatically get assembled by the FEMSystem)
    system.get_dof_map().constrain_element_vector
      ( fluid_context.get_elem_residual(),
        fluid_context.get_dof_indices(), false );

    system.rhs->add_vector( fluid_context.get_elem_residual(),
                            fluid_context.get_dof_indices() );

    if( compute_jacobian )
      {
        std::vector<libMesh::dof_id_type> velocity_dof_indices;

        const std::vector<libMesh::dof_id_type> & uf_dof_indices =
          fluid_context.get_dof_indices(this->_flow_vars.u());

        const std::vector<libMesh::dof_id_type> & vf_dof_indices =
          fluid_context.get_dof_indices(this->_flow_vars.v());

        const std::vector<libMesh::dof_id_type> * wf_dof_indices = nullptr;
        if(Dim==3)
          wf_dof_indices = &fluid_context.get_dof_indices(this->_flow_vars.w());

        libmesh_assert_equal_to( uf_dof_indices.size(), n_fluid_dofs );
        libmesh_assert_equal_to( vf_dof_indices.size(), n_fluid_dofs );
        if(Dim==3)
          libmesh_assert_equal_to( wf_dof_indices->size(), n_fluid_dofs );

        velocity_dof_indices.insert( velocity_dof_indices.end(),
                                     uf_dof_indices.begin(),
                                     uf_dof_indices.end() );

        velocity_dof_indices.insert( velocity_dof_indices.end(),
                                     vf_dof_indices.begin(),
                                     vf_dof_indices.end() );

        if(Dim==3)
          velocity_dof_indices.insert( velocity_dof_indices.end(),
                                       wf_dof_indices->begin(),
                                       wf_dof_indices->end() );

        //Build up solid dof indices
        std::vector<libMesh::dof_id_type> solid_dof_indices;

        const std::vector<libMesh::dof_id_type>& us_dof_indices =
          solid_context.get_dof_indices(this->_disp_vars.u());

        const std::vector<libMesh::dof_id_type>& vs_dof_indices =
          solid_context.get_dof_indices(this->_disp_vars.v());

        const std::vector<libMesh::dof_id_type> * ws_dof_indices = nullptr;
        if(Dim==3)
          ws_dof_indices = &solid_context.get_dof_indices(this->_disp_vars.w());

        libmesh_assert_equal_to( us_dof_indices.size(), n_solid_dofs );
        libmesh_assert_equal_to( vs_dof_indices.size(), n_solid_dofs );
        if(Dim==3)
          libmesh_assert_equal_to( ws_dof_indices->size(), n_solid_dofs );

        solid_dof_indices.insert( solid_dof_indices.end(),
                                  us_dof_indices.begin(),
                                  us_dof_indices.end() );

        solid_dof_indices.insert( solid_dof_indices.end(),
                                  vs_dof_indices.begin(),
                                  vs_dof_indices.end() );

        if(Dim==3)
          solid_dof_indices.insert( solid_dof_indices.end(),
                                    ws_dof_indices->begin(),
                                    ws_dof_indices->end() );

        //Build up lambda dof indices
        std::vector<libMesh::dof_id_type> lambda_dof_indices;

        const std::vector<libMesh::dof_id_type>& ulm_dof_indices =
          solid_context.get_dof_indices(this->_lambda_var.u());

        const std::vector<libMesh::dof_id_type>& vlm_dof_indices =
          solid_context.get_dof_indices(this->_lambda_var.v());

        const std::vector<libMesh::dof_id_type> * wlm_dof_indices = nullptr;
        if(Dim==3)
          wlm_dof_indices = &solid_context.get_dof_indices(this->_lambda_var.w());

        libmesh_assert_equal_to( ulm_dof_indices.size(), n_lambda_dofs );
        libmesh_assert_equal_to( vlm_dof_indices.size(), n_lambda_dofs );
        if(Dim==3)
          libmesh_assert_equal_to( wlm_dof_indices->size(), n_lambda_dofs );

        lambda_dof_indices.insert( lambda_dof_indices.end(),
                                   ulm_dof_indices.begin(),
                                   ulm_dof_indices.end() );

        lambda_dof_indices.insert( lambda_dof_indices.end(),
                                   vlm_dof_indices.begin(),
                                   vlm_dof_indices.end() );

        if(Dim==3)
          lambda_dof_indices.insert( lambda_dof_indices.end(),
                                     wlm_dof_indices->begin(),
                                     wlm_dof_indices->end() );

        // Since we manually built the fluid context, we have to manually
        // constrain and add the residuals and Jacobians.
        //! \todo  We're hardcoding to the case that the residual is always
        //  assembled and homogeneous constraints.

        system.get_dof_map().constrain_element_matrix
          ( Kf_s,
            velocity_dof_indices,
            solid_dof_indices, false );

        system.matrix->add_matrix( Kf_s,
                                   velocity_dof_indices,
                                   solid_dof_indices );

        system.get_dof_map().constrain_element_matrix
          ( Klm_f,
            lambda_dof_indices,
            velocity_dof_indices,
            false );

        system.matrix->add_matrix( Klm_f,
                                   lambda_dof_indices,
                                   velocity_dof_indices );

        system.get_dof_map().constrain_element_matrix
          ( Kf_lm,
            velocity_dof_indices,
            lambda_dof_indices, false );

        system.matrix->add_matrix( Kf_lm,
                                   velocity_dof_indices,
                                   lambda_dof_indices );
      } // compute_jacobian
  }

  // Instantiate
  template void FictitiousDomainFluidStructureInteractionAbstract::assemble_coupled_terms<2>
  ( bool, MultiphysicsSystem &, const AssemblyContext &, AssemblyContext &,
    unsigned int, unsigned int, unsigned int,
    libMesh::DenseMatrix<libMesh::Number> &,
    libMesh::DenseMatrix<libMesh::Number> &,
    libMesh::DenseMatrix<libMesh::Number> &);

  template void FictitiousDomainFluidStructureInteractionAbstract::assemble_coupled_terms<3>
  ( bool, MultiphysicsSystem &, const AssemblyContext &, AssemblyContext &,
    unsigned int, unsigned int, unsigned int,
    libMesh::DenseMatrix<libMesh::Number> &,
    libMesh::DenseMatrix<libMesh::Number> &,
    libMesh::DenseMatrix<libMesh::Number> &);

} // end namespace GRINS
