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

// this class
#include "grins/immersed_boundary.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/elasticity_tensor.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"
#include "grins/generic_ic_handler.h"

// includes for IBM instantiation
#include "grins/elastic_membrane.h"
#include "grins/elastic_cable.h"
#include "grins/hookes_law_1d.h"
#include "grins/hookes_law.h"
#include "grins/mooney_rivlin.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/unsteady_solver.h"

// PETSc includes
# include <petscsnes.h>

namespace GRINS
{

  template<typename SolidMech>
  ImmersedBoundary<SolidMech>::ImmersedBoundary(const std::string& physics_name,
                                                std::unique_ptr<SolidMech> & solid_mech_ptr,
                                                const GetPot& input)
    : Physics(physics_name, input),
      _flow_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<VelocityVariable>(VariablesParsing::velocity_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _disp_vars(GRINSPrivate::VariableWarehouse::get_variable_subclass<DisplacementVariable>(VariablesParsing::disp_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _lambda_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<LagrangeMultVectorVariable>(VariablesParsing::lagrange_mult_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _solid_press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _solid_mech(std::move(solid_mech_ptr)),
      _fluid_mechanics(input("Physics/ImmersedBoundary/fluid_mechanics","DIE!")),
      _solid_mechanics(input("Physics/ImmersedBoundary/solid_mechanics","DIE!")),
      _coupling_matrix(nullptr),
      _old_old_local_nonlinear_solution(nullptr)
  {
    _lambda_var.set_is_constraint_var(true);
    _solid_press_var.set_is_constraint_var(true);

    // Need to sanity check the dim() of the Variables

    // Get the fluid mechanics from the input, this is needed to query the fluid subdomain ids
    if( !input.have_variable( "Physics/ImmersedBoundary/fluid_mechanics" ) )
      libmesh_error_msg("Error: Must specify fluid_mechanics to identify the fluid in the IBM physics.\n");

    // Get the solid mechanics from the input, this is needed to prerequest data in init_context
    if( !input.have_variable( "Physics/ImmersedBoundary/solid_mechanics" ) )
      libmesh_error_msg("Error: Must specify solid_mechanics to identify the solid in the IBM physics.\n");

    // Get the solid subdomain ids from the input
    std::string solid_id_str = "Physics/"+_solid_mechanics+"/enabled_subdomains";
    unsigned int n_solid_subdomains = input.vector_variable_size(solid_id_str);
    if( n_solid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify at least one enabled solid subdomain for Immersed Boundary Physics!");

    for( unsigned int i = 0; i < n_solid_subdomains; i++)
      _solid_subdomain_set.insert( input(solid_id_str, -1, i) );

    // Get the subdomain ids for the fluid
    std::string fluid_id_str = "Physics/"+_fluid_mechanics+"/enabled_subdomains";
    unsigned int n_fluid_subdomains = input.vector_variable_size(fluid_id_str);
    if( n_fluid_subdomains == 0 )
      libmesh_error_msg("Error: Must specify at least one enabled fluid subdomain for Immersed Boundary Physics!");

    for( unsigned int i = 0; i < n_fluid_subdomains; i++)
      _fluid_subdomain_set.insert( input(fluid_id_str, -1, i) );

    // TODO: Need to check that Mesh has all the fluid and solid subdomain ids

    this->_ic_handler = new GenericICHandler(physics_name, input);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::auxiliary_init( MultiphysicsSystem & system )
  {
    // Setup the CouplingMatrix that will be used repeatedly in the
    // CouplingFunctor
    _coupling_matrix.reset( new libMesh::CouplingMatrix(system.n_vars()) );
    this->setup_coupling_matrix( _flow_vars, _disp_vars, _lambda_var, _press_var, *_coupling_matrix );

    // Get the point locator object that will find the right fluid element
    _point_locator = system.get_mesh().sub_point_locator();

    // Build helper FEMContexts. We'll use this to handle
    // supplemental finite element data for variables that
    // we need, but are not defined on the "current" subdomain.
    {
      std::unique_ptr<libMesh::DiffContext> raw_context = system.build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      _fluid_context.reset(context);
    }

    libMesh::NumericVector<libMesh::Number> & old_old_nonlinear_soln =
      system.add_vector("_old_old_nonlinear_solution");

    if(!old_old_nonlinear_soln.initialized())
      old_old_nonlinear_soln.init(system.n_dofs(), system.n_local_dofs(), false, libMesh::PARALLEL);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_ics( libMesh::FEMSystem* system,
                                              libMesh::CompositeFunction<libMesh::Number>& all_ics )
  {
    libmesh_assert(this->_solid_mech);
    _solid_mech->init_ics(system,all_ics);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Velocity are first order in time
    system->time_evolving(_flow_vars.u(),1);
    system->time_evolving(_flow_vars.v(),1);

    if ( _flow_vars.dim() == 3 )
      system->time_evolving(_flow_vars.w(),1);


    // Displacements are second order in time
    system->time_evolving(_disp_vars.u(),1);

    if( _disp_vars.dim() >= 2 )
      system->time_evolving(_disp_vars.v(),1);

    if ( _disp_vars.dim() == 3 )
      system->time_evolving(_disp_vars.w(),1);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.

    if (_solid_mechanics == "ElasticMembrane")
      {
        //membrane context inits
        context.get_element_fe(_disp_vars.u(),2)->get_JxW();
        context.get_element_fe(_disp_vars.u(),2)->get_phi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphidxi();
        context.get_element_fe(_disp_vars.u(),2)->get_dphideta();
        //for metric tensors of the membrane
        context.get_element_fe(_disp_vars.u(),2)->get_dxyzdxi();
        context.get_element_fe(_disp_vars.u(),2)->get_dxyzdeta();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidx();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidy();
        context.get_element_fe(_disp_vars.u(),2)->get_dxidz();
        context.get_element_fe(_disp_vars.u(),2)->get_detadx();
        context.get_element_fe(_disp_vars.u(),2)->get_detady();
        context.get_element_fe(_disp_vars.u(),2)->get_detadz();
      }
    else if (_solid_mechanics == "ElasticCable")
      {
        //cable context inits
        context.get_element_fe(_disp_vars.u(),1)->get_JxW();
        context.get_element_fe(_disp_vars.u(),1)->get_phi();
        context.get_element_fe(_disp_vars.u(),1)->get_dphidxi();
        //for metric tensors of cable
        context.get_element_fe(_disp_vars.u(),1)->get_dxyzdxi();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidx();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidy();
        context.get_element_fe(_disp_vars.u(),1)->get_dxidz();
      }
    else
      {
        std::string err = "ERROR: solid_mechanics not properly specified";
        libmesh_error_msg(err);
      }

    context.get_element_fe(_solid_press_var.p(),2)->get_JxW();
    context.get_element_fe(_solid_press_var.p(),2)->get_phi();

    context.get_element_fe( _lambda_var.u(),2)->get_dphi();
    context.get_element_fe( _lambda_var.u(),2)->get_phi();
    context.get_element_fe( _lambda_var.u(),2)->get_JxW();
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::presolve( MultiphysicsSystem & system )
  {
    // We need to rebuild the overlap map each time because the position
    // of the solid can change in between each Newton step.
    _fluid_solid_overlap.reset( new OverlappingFluidSolidMap(system,
                                                             (*_point_locator),
                                                             _solid_subdomain_set,
                                                             _fluid_subdomain_set,
                                                             _disp_vars) );

    libMesh::DofMap & dof_map = system.get_dof_map();

    if(!_coupling_functor)
      {
        _coupling_functor.reset( new ImmersedBoundaryCouplingFunctor(system.get_mesh(),
                                                                     *_coupling_matrix,
                                                                     *_fluid_solid_overlap) );
        dof_map.add_coupling_functor(*_coupling_functor);
      }
    else
      {
        dof_map.remove_coupling_functor(*_coupling_functor);
        _coupling_functor.reset( new ImmersedBoundaryCouplingFunctor(system.get_mesh(),
                                                                     *_coupling_matrix,
                                                                     *_fluid_solid_overlap) );
        dof_map.add_coupling_functor(*_coupling_functor);
      }

    // Handle algebraic coupling
    // We have to manually make these calls since we couldn't attach this
    // functor before the mesh and dof_map was prepared.
    dof_map.reinit_send_list(system.get_mesh());

    this->reinit_ghosted_vectors( system );

    // Handle full coupling (sparsity pattern)
    dof_map.clear_sparsity();
    dof_map.compute_sparsity(system.get_mesh());

    // Now to reinit the matrix since we changed the sparsity pattern
    libMesh::SparseMatrix<libMesh::Number> & matrix = system.get_matrix("System Matrix");
    libmesh_assert(dof_map.is_attached(matrix));
    matrix.init();
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::preadvance_timestep( MultiphysicsSystem & system )
  {
    libMesh::NumericVector<libMesh::Number> & old_old_nonlinear_soln =
      system.get_vector("_old_old_nonlinear_solution");

    libMesh::NumericVector<libMesh::Number> & old_nonlinear_soln =
      system.get_vector("_old_nonlinear_solution");

    // timestep n copied to timestep n-1 now
    old_old_nonlinear_soln = old_nonlinear_soln;

    old_old_nonlinear_soln.localize( *_old_old_local_nonlinear_solution,
                                     system.get_dof_map().get_send_list() );
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::reinit_ghosted_vectors( MultiphysicsSystem & system )
  {
    const libMesh::DofMap & dof_map = system.get_dof_map();

    // Update current local solution
    system.current_local_solution = libMesh::NumericVector<libMesh::Number>::build(system.comm());

    system.current_local_solution->init(system.n_dofs(), system.n_local_dofs(),
                                        dof_map.get_send_list(), false,
                                        libMesh::GHOSTED);

    system.solution->localize(*(system.current_local_solution),dof_map.get_send_list());

    // Also need to update unsteady system vectors
    libMesh::TimeSolver & time_solver_raw = system.get_time_solver();

    libMesh::UnsteadySolver * unsteady_solver = dynamic_cast<libMesh::UnsteadySolver *>(&time_solver_raw);
    if( unsteady_solver)
      {
        // Old nonlinear solution
        {
          unsteady_solver->old_local_nonlinear_solution =
            libMesh::NumericVector<libMesh::Number>::build(system.comm());

          unsteady_solver->old_local_nonlinear_solution->init(system.n_dofs(), system.n_local_dofs(),
                                                              dof_map.get_send_list(), false,
                                                              libMesh::GHOSTED);

          libMesh::NumericVector<libMesh::Number> & old_nonlinear_soln =
            system.get_vector("_old_nonlinear_solution");

          old_nonlinear_soln.localize(*(unsteady_solver->old_local_nonlinear_solution), dof_map.get_send_list());
        }

        // Old Old nonlinear solution
        {
          _old_old_local_nonlinear_solution =
            libMesh::NumericVector<libMesh::Number>::build(system.comm());

          _old_old_local_nonlinear_solution->init(system.n_dofs(), system.n_local_dofs(),
                                                  dof_map.get_send_list(), false,
                                                  libMesh::GHOSTED);

          libMesh::NumericVector<libMesh::Number> & old_old_nonlinear_soln =
            system.get_vector("_old_old_nonlinear_solution");

          old_old_nonlinear_soln.localize(*_old_old_local_nonlinear_solution, dof_map.get_send_list());
        }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::reinit( MultiphysicsSystem & system )
  {
    _point_locator.reset();
    _point_locator = system.get_mesh().sub_point_locator();
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::element_time_derivative( bool compute_jacobian,
                                                             AssemblyContext & context )
  {
    const libMesh::Elem & solid_elem = context.get_elem();

    // Only compute this if we are on a solid element
    if( this->is_solid_elem( solid_elem.subdomain_id() ) &&
        this->_fluid_solid_overlap->has_overlapping_fluid_elem(solid_elem.id()) )
      {
        // For clarity
        AssemblyContext & solid_context = context;

        unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
        unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

        MultiphysicsSystem & system = context.get_multiphysics_system();

        libMesh::DenseMatrix<libMesh::Number> Kf_s;
        libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kf_s), Kuf_vs(Kf_s);
        libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kf_s), Kvf_vs(Kf_s);

        libMesh::DenseMatrix<libMesh::Number> Klm_f;
        libMesh::DenseSubMatrix<libMesh::Number> Kulm_uf(Klm_f), Kulm_vf(Klm_f);
        libMesh::DenseSubMatrix<libMesh::Number> Kvlm_uf(Klm_f), Kvlm_vf(Klm_f);

        libMesh::DenseMatrix<libMesh::Number> Kf_lm;
        libMesh::DenseSubMatrix<libMesh::Number> Kuf_ulm(Kf_lm), Kuf_vlm(Kf_lm);
        libMesh::DenseSubMatrix<libMesh::Number> Kvf_ulm(Kf_lm), Kvf_vlm(Kf_lm);

        // We need to grab the fluid elements that are overlapping with this solid elem.
        // Then, for that fluid element, extract the indices of the *solid* quadrature points
        // that are in that fluid element
        const std::set<libMesh::dof_id_type> & fluid_elem_ids =
          _fluid_solid_overlap->get_overlapping_fluid_elems(solid_elem.id());

        // Loop over those fluid elements
        for( const auto & fluid_id : fluid_elem_ids )
          {
            const libMesh::Elem * fluid_elem = system.get_mesh().elem_ptr(fluid_id);

            const std::vector<unsigned int> & quad_points =
              _fluid_solid_overlap->get_solid_qps(solid_elem.id(), fluid_elem->id() );

            (this->_fluid_context)->pre_fe_reinit(system,fluid_elem);

            unsigned int n_fluid_dofs = (this->_fluid_context)->get_dof_indices(this->_flow_vars.u()).size();

            // Prepare the space for Jacobians first in case we're going to compute them
            // analytically
            if ( compute_jacobian )
              {
                this->prepare_jacobians(n_fluid_dofs, n_solid_dofs, n_lambda_dofs,
                                        Kf_s,
                                        Kuf_us,Kuf_vs,Kvf_us,Kvf_vs,
                                        Klm_f,
                                        Kulm_uf,Kulm_vf, Kvlm_uf,Kvlm_vf,
                                        Kf_lm,
                                        Kuf_ulm,Kuf_vlm,Kvf_ulm,Kvf_vlm);
              }

            // Compute and assemble residuals for all the quadrature points in this fluid element
            // and also assemble analytical Jacobians if we're using analytical Jacobians
            this->compute_ibm_residuals(system,solid_context,*(this->_fluid_context),quad_points);


            // Now that we've looped over all the quadrature points on this fluid element
            // we can now assemble the fluid residual (the residuals in the solid context
            // automatically get assembled by the FEMSystem)
            system.get_dof_map().constrain_element_vector
              ( (this->_fluid_context)->get_elem_residual(),
                (this->_fluid_context)->get_dof_indices(), false );

            system.rhs->add_vector( (this->_fluid_context)->get_elem_residual(),
                                    (this->_fluid_context)->get_dof_indices() );


            // If need be, we can now also compute numerical jacobians and assemble them
            if ( compute_jacobian )
              {
                 this->compute_numerical_jacobians( system,solid_context,*(this->_fluid_context),
                                                    quad_points,
                                                    Kuf_ulm,Kuf_vlm,Kvf_ulm,Kvf_vlm,
                                                    Kulm_uf,Kulm_vf,Kvlm_uf,Kvlm_vf,
                                                    Kuf_us,Kuf_vs,Kvf_us,Kvf_vs);

                 /*
                this->compute_analytic_jacobians(system,solid_context,*(this->_fluid_context),
                                                 quad_points,
                                                 Kuf_ulm,Kvf_vlm,
                                                 Kulm_uf,Kvlm_vf);
                 */

                 this->assemble_fluid_jacobians(system,solid_context,*(this->_fluid_context),
                                                n_fluid_dofs, n_solid_dofs, n_lambda_dofs,
                                                Kf_s,Klm_f,Kf_lm);

              } // if compute_jacobian


          } // loop over fluid elems


      } // if(solid_elem)
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_ibm_residuals(const MultiphysicsSystem & system,
                                                          AssemblyContext & solid_context,
                                                          libMesh::FEMContext & fluid_context,
                                                          const std::vector<unsigned int> & quad_points )
  {
    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();

    unsigned int lambda_x = this->_lambda_var.u();
    unsigned int lambda_y = this->_lambda_var.v();

    unsigned int p_var = this->_solid_press_var.p();

    // Fluid residual data structures
    libMesh::DenseSubVector<libMesh::Number> & Fuf =
      fluid_context.get_elem_residual(this->_flow_vars.u());

    libMesh::DenseSubVector<libMesh::Number> & Fvf =
      fluid_context.get_elem_residual(this->_flow_vars.v());

    libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(u_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(v_var);

    libMesh::DenseSubVector<libMesh::Number> & Fulm = solid_context.get_elem_residual(lambda_x);
    libMesh::DenseSubVector<libMesh::Number> & Fvlm = solid_context.get_elem_residual(lambda_y);

    libMesh::DenseSubVector<libMesh::Number> & Fp = solid_context.get_elem_residual(p_var);

    const std::vector<libMesh::Point> & solid_qpoints = solid_context.get_element_fe(u_var,2)->get_xyz();

    const libMesh::Elem & fluid_elem = fluid_context.get_elem();

    for( const auto & qp : quad_points )
      {

        // Reinit fluid context for current solid quadrature point
        {
          const libMesh::Point & x_qp =  solid_qpoints[qp];
          libMesh::Point x = this->compute_displaced_point(system,solid_context,x_qp,qp);

          libMesh::FEBase * fe = fluid_context.get_element_fe(this->_flow_vars.u());
          libMesh::FEType fetype = fe->get_fe_type();

          //We need to hand *reference* element points to the FEMContext to reinit
          unsigned int dim = 2;
          libMesh::Point x_ref = libMesh::FEInterface::inverse_map(dim,fetype,&fluid_elem,x);

          std::vector<libMesh::Point> x_ref_vec(1,x_ref);

          fluid_context.elem_fe_reinit(&x_ref_vec);
        }

        this->compute_residuals(solid_context,fluid_context,qp,
                                Fuf,Fvf,Fus,Fvs,Fulm,Fvlm,Fp);

      } // loop over quadrature points on the current fluid element
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_residuals( AssemblyContext & solid_context,
                                                       const libMesh::FEMContext & fluid_context,
                                                       unsigned int sqp,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fuf,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fvf,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fus,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fvs,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fulm,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fvlm,
                                                       libMesh::DenseSubVector<libMesh::Number> & Fp)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();
    unsigned int n_press_dofs = solid_context.get_dof_indices(this->_solid_press_var.p()).size();

    libMesh::Real udot, vdot;
    solid_context.interior_rate(this->_disp_vars.u(), sqp, udot);
    solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot);

    libMesh::Gradient grad_u, grad_v;
    solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_u);
    solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_v);

    libMesh::Gradient grad_udot, grad_vdot;
    solid_context.interior_rate_gradient(this->_disp_vars.u(), sqp, grad_udot);
    solid_context.interior_rate_gradient(this->_disp_vars.v(), sqp, grad_vdot);

    libMesh::Real lambda_x, lambda_y;
    solid_context.interior_value(this->_lambda_var.u(), sqp, lambda_x);
    solid_context.interior_value(this->_lambda_var.v(), sqp, lambda_y);

    libMesh::Gradient grad_lambda_x, grad_lambda_y;
    solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lambda_x);
    solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lambda_y);

    libMesh::Real Vx, Vy;
    fluid_context.interior_value(this->_flow_vars.u(), 0, Vx);
    fluid_context.interior_value(this->_flow_vars.v(), 0, Vy);

    libMesh::Gradient grad_Vx, grad_Vy;
    fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx);
    fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy);

    //Solid Pressure value and shape function
    libMesh::Real p;
    solid_context.interior_value(this->_solid_press_var.p(), sqp, p);

    const std::vector<std::vector<libMesh::Real> > p_phi =
      solid_context.get_element_fe(this->_solid_press_var.p())->get_phi();

    libMesh::TensorValue<libMesh::Real> F;
    this->eval_deform_gradient(grad_u,grad_v,F);

    libMesh::TensorValue<libMesh::Real> Fold;
    this->compute_prev_timestep_deform_gradient(solid_context,sqp,Fold);

    libMesh::TensorValue<libMesh::Real> Fdot;
    this->eval_deform_grad_rate(grad_udot,grad_vdot,Fdot);

    libMesh::TensorValue<libMesh::Real> P;
    this->eval_first_Piola(grad_u,grad_v,P);

    const std::vector<std::vector<libMesh::Real> > fluid_phi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

    const std::vector<std::vector<libMesh::Real> > solid_phi =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > solid_dphi =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

    const std::vector<std::vector<libMesh::Real> > lambda_phi =
      solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & lambda_dphi =
      solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

    const std::vector<libMesh::Real> & solid_JxW =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_JxW();

    libMesh::Real jac = solid_JxW[sqp];

    libMesh::TensorValue<libMesh::Real> fdphi_times_F;

    libMesh::Tensor grad_lam( grad_lambda_x(0), grad_lambda_x(1), 0,
                              grad_lambda_y(0), grad_lambda_y(1), 0,
                              0, 0, 0);

    libMesh::Tensor grad_lam_timesFT( grad_lam*(Fold.transpose()) );

    libMesh::Tensor gradV( grad_Vx(0),  grad_Vx(1), 0,
                           grad_Vy(0), grad_Vy(1), 0,
                           0, 0, 0);

    libMesh::TensorValue<libMesh::Real> gradV_times_F( gradV*Fold );

    libMesh::Tensor FT_times_gradV( (Fold.transpose())*gradV );

    libMesh::TensorValue<libMesh::Real> Ftrans = F.transpose();
    libMesh::TensorValue<libMesh::Real> C(Ftrans*F);
    libMesh::TensorValue<libMesh::Real> Cinv = C.inverse();
    libMesh::Number I1 = C.tr();
    libMesh::Number J = F.det();
    libMesh::Number Up = std::log(J)/J;

    libMesh::Tensor F_times_Cinv( F*Cinv );

    libMesh::Real J23 = std::pow(J,-2/3);
    libMesh::Real mus = 5;

    // Fluid residual
    for (unsigned int i=0; i != n_fluid_dofs; i++)
      {
	libmesh_assert_equal_to( fluid_phi[i].size(), 1 );

        libMesh::Gradient fluid_term = grad_lam_timesFT*fluid_dphi[i][0];

        // L2 Norm
	Fuf(i) -= lambda_x*fluid_phi[i][0]*jac;
	Fvf(i) -= lambda_y*fluid_phi[i][0]*jac;

        // H1 Term
        Fuf(i) -= fluid_term(0)*jac;
        Fvf(i) -= fluid_term(1)*jac;

	//Computing fdphi_times_F
        /*
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    for( unsigned int beta = 0; beta < 2; beta++ )
	      {
		fdphi_times_F(0,alpha) += fluid_dphi[i][0](beta)*Fold(beta,alpha);
		fdphi_times_F(1,alpha) += fluid_dphi[i][0](beta)*Fold(beta,alpha);
	      }
   	  }
        */

	// H1 Norm
        /*
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fuf(i) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha))*jac;
	    Fvf(i) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha))*jac;
	  }
        */
      }

    // Solid residual
    for (unsigned int i=0; i != n_solid_dofs; i++)
      {
        //L2 Norm
	Fus(i) += lambda_x*solid_phi[i][sqp]*jac;
	Fvs(i) += lambda_y*solid_phi[i][sqp]*jac;
	/*
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) -= P(0,alpha)*solid_dphi[i][sqp](alpha)*jac;
	    Fvs(i) -= P(1,alpha)*solid_dphi[i][sqp](alpha)*jac;
	  }
	*/

	//pressure term
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) -= p*J*F_times_Cinv(0,alpha)*solid_dphi[i][sqp](alpha)*jac;
	    Fvs(i) -= p*J*F_times_Cinv(1,alpha)*solid_dphi[i][sqp](alpha)*jac;
	  }

	//dWdC term
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) -= mus*J23*(F(0,alpha)-(1/3)*I1*F_times_Cinv(0,alpha))*solid_dphi[i][sqp](alpha)*jac;
	    Fvs(i) -= mus*J23*(F(1,alpha)-(1/3)*I1*F_times_Cinv(1,alpha))*solid_dphi[i][sqp](alpha)*jac;
	  }

	//H1 Norm
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) += grad_lambda_x(alpha)*solid_dphi[i][sqp](alpha)*jac;
	    Fvs(i) += grad_lambda_y(alpha)*solid_dphi[i][sqp](alpha)*jac;
	  }
      }

    // Lambda residual
    for( unsigned int i = 0; i < n_lambda_dofs; i++ )
      {
        //L2 Norm
        Fulm(i) += lambda_phi[i][sqp]*(Vx - udot)*jac;
	Fvlm(i) += lambda_phi[i][sqp]*(Vy - vdot)*jac;

        //libMesh::Gradient fluid_term = (FT_times_gradV - Fdot)*lambda_dphi[i][sqp];
        libMesh::Gradient fluid_term = (gradV_times_F-Fdot)*lambda_dphi[i][sqp];

        Fulm(i) += fluid_term(0)*jac;
        Fvlm(i) += fluid_term(1)*jac;

	//Computing gradV_times_F
        /*
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    for( unsigned int beta = 0; beta < 2; beta++ )
	      {
		gradV_times_F(0,alpha) += grad_Vx(beta)*Fold(beta,alpha);
		gradV_times_F(1,alpha) += grad_Vy(beta)*Fold(beta,alpha);
	      }
	  }
        */

	//H1 Norm
        /*
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	{
	  Fulm(i) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha)))*jac;
	  Fvlm(i) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha)))*jac;
	}
        */

      }

    // Solid pressure residual
    for( unsigned int i = 0; i < n_press_dofs; i++ )
      {
        Fp(i) += Up*p_phi[i][sqp]*jac;
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_numerical_jacobians(const MultiphysicsSystem & system,
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
                                                                libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vs)
  {
    // Cache original residuals
    libMesh::DenseVector<libMesh::Number> original_solid_residual(solid_context.get_elem_residual());
    libMesh::DenseVector<libMesh::Number> original_fluid_residual(fluid_context.get_elem_residual());

    // Create space for storing residuals with backward-perturbed solutions
    libMesh::DenseVector<libMesh::Number> backwards_solid_residual(solid_context.get_elem_residual());
    libMesh::DenseVector<libMesh::Number> backwards_fluid_residual(fluid_context.get_elem_residual());

    libMesh::Real delta = 1.0e-8;

    // Compute lambdax derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & lambda_xcoeff =
        solid_context.get_elem_solution(this->_lambda_var.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kus_ulm =
        solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_lambda_var.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_ulm =
        solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_lambda_var.u());

      this->compute_lambda_derivs(system,quad_points,solid_context,fluid_context,delta,
                                  backwards_solid_residual,backwards_fluid_residual,
                                  lambda_xcoeff,
                                  Kuf_ulm,Kvf_ulm,Kus_ulm,Kvs_ulm);
    }

    // Compute lambday derivs

    {
      libMesh::DenseSubVector<libMesh::Number> & lambda_ycoeff =
        solid_context.get_elem_solution(this->_lambda_var.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kus_vlm =
        solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_lambda_var.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vlm =
        solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_lambda_var.v());

      this->compute_lambda_derivs(system,quad_points,solid_context,fluid_context,delta,
                                  backwards_solid_residual,backwards_fluid_residual,
                                  lambda_ycoeff,
                                  Kuf_vlm,Kvf_vlm,Kus_vlm,Kvs_vlm);
    }


    // Compute fluidx derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & fluid_ucoeff =
        fluid_context.get_elem_solution(this->_flow_vars.u());

      this->compute_fluid_derivs(system,quad_points,solid_context,fluid_context,delta,
                                 backwards_solid_residual,backwards_fluid_residual,
                                 fluid_ucoeff,
                                 Kulm_uf,Kvlm_uf);
    }

    // Compute fluidy derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & fluid_vcoeff =
        fluid_context.get_elem_solution(this->_flow_vars.v());

      this->compute_fluid_derivs(system,quad_points,solid_context,fluid_context,delta,
                                 backwards_solid_residual,backwards_fluid_residual,
                                 fluid_vcoeff,
                                 Kulm_vf,Kvlm_vf);
    }

    // Compute solidx derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & u_coeffs =
        solid_context.get_elem_solution(this->_disp_vars.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kus_us =
        solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_us =
        solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kulm_us =
        solid_context.get_elem_jacobian(this->_lambda_var.u(),this->_disp_vars.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_us =
        solid_context.get_elem_jacobian(this->_lambda_var.v(),this->_disp_vars.u());

      libMesh::DenseSubMatrix<libMesh::Number> & Kp_us =
        solid_context.get_elem_jacobian(this->_solid_press_var.p(),this->_disp_vars.u());

      this->compute_solid_derivs(system,quad_points,solid_context,fluid_context,delta,
                                 backwards_solid_residual,backwards_fluid_residual,
                                 u_coeffs,
                                 Kuf_us,Kvf_us,Kus_us,Kvs_us,Kulm_us,Kvlm_us,Kp_us);
    }


    // Compute solidy derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & v_coeffs =
        solid_context.get_elem_solution(this->_disp_vars.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kus_vs =
        solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vs =
        solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vs =
        solid_context.get_elem_jacobian(this->_lambda_var.u(),this->_disp_vars.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vs =
        solid_context.get_elem_jacobian(this->_lambda_var.v(),this->_disp_vars.v());

      libMesh::DenseSubMatrix<libMesh::Number> & Kp_vs =
        solid_context.get_elem_jacobian(this->_solid_press_var.p(),this->_disp_vars.v());

      this->compute_solid_derivs(system,quad_points,solid_context,fluid_context,delta,
                                 backwards_solid_residual,backwards_fluid_residual,
                                 v_coeffs,
                                 Kuf_vs,Kvf_vs,Kus_vs,Kvs_vs,Kulm_vs,Kvlm_vs,Kp_vs);

    }

    // Compute p derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & p_coeffs =
        solid_context.get_elem_solution(this->_solid_press_var.p());

      libMesh::DenseSubMatrix<libMesh::Number> & Kus_p =
        solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_solid_press_var.p());

      libMesh::DenseSubMatrix<libMesh::Number> & Kvs_p =
        solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_solid_press_var.p());

      this->compute_press_derivs(system,quad_points,solid_context,fluid_context,delta,
                                 backwards_solid_residual,backwards_fluid_residual,
                                 p_coeffs,Kus_p,Kvs_p);

    }

    // Restore original residuals
    solid_context.get_elem_residual() = original_solid_residual;
    fluid_context.get_elem_residual() = original_fluid_residual;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_analytic_jacobians(const MultiphysicsSystem & system,
                                                               AssemblyContext & solid_context,
                                                               libMesh::FEMContext & fluid_context,
                                                               const std::vector<unsigned int> & quad_points,
                                                               libMesh::DenseSubMatrix<libMesh::Number> & Kuf_ulm,
                                                               libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm,
                                                               libMesh::DenseSubMatrix<libMesh::Number> & Kulm_uf,
                                                               libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vf)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    libMesh::DenseSubMatrix<libMesh::Number> & Kus_ulm =
      solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_lambda_var.u());

    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vlm =
      solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_lambda_var.v());

    libMesh::DenseSubMatrix<libMesh::Number> & Kus_us =
      solid_context.get_elem_jacobian(this->_disp_vars.u(),this->_disp_vars.u());

    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vs =
      solid_context.get_elem_jacobian(this->_disp_vars.v(),this->_disp_vars.v());

    libMesh::DenseSubMatrix<libMesh::Number> & Kulm_us =
      solid_context.get_elem_jacobian(this->_lambda_var.u(),this->_disp_vars.u());

    libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vs =
      solid_context.get_elem_jacobian(this->_lambda_var.v(),this->_disp_vars.v());

    const std::vector<libMesh::Point> & solid_qpoints =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_xyz();

    const libMesh::Elem & fluid_elem = fluid_context.get_elem();

    for( const auto & qp : quad_points )
      {
        // Reinit fluid context for current solid quadrature point
        {
          const libMesh::Point & x_qp =  solid_qpoints[qp];
          libMesh::Point x = this->compute_displaced_point(system,solid_context,x_qp,qp);

          libMesh::FEBase * fe = fluid_context.get_element_fe(this->_flow_vars.u());
          libMesh::FEType fetype = fe->get_fe_type();

          //We need to hand *reference* element points to the FEMContext to reinit
          unsigned int dim = 2;
          libMesh::Point x_ref = libMesh::FEInterface::inverse_map(dim,fetype,&fluid_elem,x);

          std::vector<libMesh::Point> x_ref_vec(1,x_ref);

          fluid_context.elem_fe_reinit(&x_ref_vec);
        }

        const std::vector<libMesh::Real> & solid_JxW =
          solid_context.get_element_fe(this->_disp_vars.u(),2)->get_JxW();

        libMesh::Real jac = solid_JxW[qp];

        const std::vector<std::vector<libMesh::Real> > fluid_phi =
          fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();

        const std::vector<std::vector<libMesh::Real> > solid_phi =
          solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();

        const std::vector<std::vector<libMesh::RealGradient> > solid_dphi =
          solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

        const std::vector<std::vector<libMesh::Real> > lambda_phi =
          solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();

        // Solid derivatives
        for( unsigned int j = 0; j < n_solid_dofs; j++ )
          {
            // Lambda residual
            for( unsigned int i = 0; i < n_lambda_dofs; i++ )
              {
                const libMesh::Real value =
                  lambda_phi[i][qp]*solid_phi[j][qp]*jac*solid_context.get_elem_solution_rate_derivative();

                Kulm_us(i,j) -= value;
                Kvlm_vs(i,j) -= value;
              }

            // Solid residual
            for (unsigned int i = 0; i != n_solid_dofs; i++)
              {
                const libMesh::Real value =
                  5*solid_dphi[j][qp]*solid_dphi[i][qp]*jac*solid_context.get_elem_solution_derivative();

                Kus_us(i,j) -= value;
                Kvs_vs(i,j) -= value;
              }

          }

        // Fluid derivatives
        for( unsigned int j = 0; j < n_fluid_dofs; j++ )
          {
            libmesh_assert_equal_to( fluid_phi[j].size(), 1 );

            // Lambda residual
            for (unsigned int i=0; i != n_lambda_dofs; i++)
              {
                Kulm_uf(i,j) += lambda_phi[i][qp]*fluid_phi[j][0]*jac*solid_context.get_elem_solution_derivative();
                Kvlm_vf(i,j) += lambda_phi[i][qp]*fluid_phi[j][0]*jac*solid_context.get_elem_solution_derivative();
              }
          }

        // Lambda deriviatives
        for( unsigned int j = 0; j < n_lambda_dofs; j++ )
          {
            const libMesh::Real value = lambda_phi[j][qp]*jac*solid_context.get_elem_solution_derivative();

            // Fluid residual
            for (unsigned int i=0; i != n_fluid_dofs; i++)
              {
                libmesh_assert_equal_to( fluid_phi[i].size(), 1 );

                Kuf_ulm(i,j) += -value*fluid_phi[i][0];
                Kvf_vlm(i,j) += -value*fluid_phi[i][0];
              }

            // Solid residual
            for (unsigned int i=0; i != n_solid_dofs; i++)
              {
                Kus_ulm(i,j) += value*solid_phi[i][qp];
                Kvs_vlm(i,j) += value*solid_phi[i][qp];
              }
          }

      } // end qp loop
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::finite_difference_residuals
  (const MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   libMesh::FEMContext & fluid_context,
   const libMesh::Real delta,
   libMesh::Number & coeff,
   libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
   libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual)
  {
    // Zero out then compute backward residuals
    solid_context.get_elem_residual().zero();
    fluid_context.get_elem_residual().zero();

    coeff -= delta;
    solid_context.recompute_elem_solution_rate(system);

    this->compute_ibm_residuals(system,solid_context,fluid_context,quad_points);

    backwards_solid_residual = solid_context.get_elem_residual();
    backwards_fluid_residual = fluid_context.get_elem_residual();

    // Zero out and compute forward residuals in place
    solid_context.get_elem_residual().zero();
    fluid_context.get_elem_residual().zero();

    coeff += 2.0*delta;
    solid_context.recompute_elem_solution_rate(system);

    this->compute_ibm_residuals(system,solid_context,fluid_context,quad_points);

    // Now store the finite difference residuals in place
    solid_context.get_elem_residual().add( -1.0, backwards_solid_residual);
    solid_context.get_elem_residual().scale( 1.0/(2.0*delta) );

    fluid_context.get_elem_residual().add( -1.0, backwards_fluid_residual);
    fluid_context.get_elem_residual().scale( 1.0/(2.0*delta) );

    // Restore coeff
    coeff -= delta;
    solid_context.recompute_elem_solution_rate(system);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_lambda_derivs
  (const MultiphysicsSystem & system,
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
   libMesh::DenseSubMatrix<libMesh::Number> & Kvs)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    libmesh_assert_equal_to(Kuf.m(),n_fluid_dofs);
    libmesh_assert_equal_to(Kvf.m(),n_fluid_dofs);
    libmesh_assert_equal_to(Kus.m(),n_solid_dofs);
    libmesh_assert_equal_to(Kvs.m(),n_solid_dofs);

    libmesh_assert_equal_to(Kuf.n(),n_lambda_dofs);
    libmesh_assert_equal_to(Kvf.n(),n_lambda_dofs);
    libmesh_assert_equal_to(Kus.n(),n_lambda_dofs);
    libmesh_assert_equal_to(Kvs.n(),n_lambda_dofs);

    for( unsigned int j = 0; j < n_lambda_dofs; j++ )
      {
        // Finite differenced residual is sitting in the solid and fluid context elem_residuals
        // after this call
        this->finite_difference_residuals(system,quad_points,solid_context,fluid_context,delta,lambda_coeff(j),
                                          backwards_solid_residual,backwards_fluid_residual);

        for (unsigned int i=0; i != n_solid_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(this->_disp_vars.u());
            libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(this->_disp_vars.v());

            Kus(i,j) += Fus(i);
            Kvs(i,j) += Fvs(i);
          }

        for (unsigned int i=0; i != n_fluid_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fuf = fluid_context.get_elem_residual(this->_flow_vars.u());
            libMesh::DenseSubVector<libMesh::Number> & Fvf = fluid_context.get_elem_residual(this->_flow_vars.v());

            Kuf(i,j) += Fuf(i);
            Kvf(i,j) += Fvf(i);
          }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_fluid_derivs
  (const MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   libMesh::FEMContext & fluid_context,
   const libMesh::Real delta,
   libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
   libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
   libMesh::DenseSubVector<libMesh::Number> & fluid_coeff,
   libMesh::DenseSubMatrix<libMesh::Number> & Kulm,
   libMesh::DenseSubMatrix<libMesh::Number> & Kvlm)
  {
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    libmesh_assert_equal_to(Kulm.m(),n_lambda_dofs);
    libmesh_assert_equal_to(Kvlm.m(),n_lambda_dofs);

    libmesh_assert_equal_to(Kulm.n(),n_fluid_dofs);
    libmesh_assert_equal_to(Kvlm.n(),n_fluid_dofs);

    for( unsigned int j = 0; j < n_fluid_dofs; j++ )
      {
        // Finite differenced residual is sitting in the solid and fluid context elem_residuals
        // after this call
        this->finite_difference_residuals(system,quad_points,solid_context,fluid_context,delta,fluid_coeff(j),
                                          backwards_solid_residual,backwards_fluid_residual);

        for (unsigned int i=0; i != n_lambda_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fulm =
              solid_context.get_elem_residual(this->_lambda_var.u());

            libMesh::DenseSubVector<libMesh::Number> & Fvlm =
              solid_context.get_elem_residual(this->_lambda_var.v());

            Kulm(i,j) += Fulm(i);
            Kvlm(i,j) += Fvlm(i);
          }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_solid_derivs
  (const MultiphysicsSystem & system,
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
   libMesh::DenseSubMatrix<libMesh::Number> & Kp)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();
    unsigned int n_press_dofs = solid_context.get_dof_indices(this->_solid_press_var.p()).size();

    libmesh_assert_equal_to(Kuf.m(),n_fluid_dofs);
    libmesh_assert_equal_to(Kvf.m(),n_fluid_dofs);
    libmesh_assert_equal_to(Kus.m(),n_solid_dofs);
    libmesh_assert_equal_to(Kvs.m(),n_solid_dofs);
    libmesh_assert_equal_to(Kulm.m(),n_lambda_dofs);
    libmesh_assert_equal_to(Kvlm.m(),n_lambda_dofs);
    libmesh_assert_equal_to(Kp.m(),n_press_dofs);

    libmesh_assert_equal_to(Kuf.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kvf.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kus.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kvs.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kulm.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kvlm.n(),n_solid_dofs);
    libmesh_assert_equal_to(Kp.n(),n_solid_dofs);

    for( unsigned int j = 0; j < n_solid_dofs; j++ )
      {
        // Finite differenced residual is sitting in the solid and fluid context elem_residuals
        // after this call
        this->finite_difference_residuals(system,quad_points,solid_context,fluid_context,delta,solid_coeff(j),
                                          backwards_solid_residual,backwards_fluid_residual);

        for (unsigned int i=0; i != n_fluid_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fuf = fluid_context.get_elem_residual(this->_flow_vars.u());
            libMesh::DenseSubVector<libMesh::Number> & Fvf = fluid_context.get_elem_residual(this->_flow_vars.v());

            Kuf(i,j) += Fuf(i);
            Kvf(i,j) += Fvf(i);
          }

        for (unsigned int i=0; i != n_solid_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(this->_disp_vars.u());
            libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(this->_disp_vars.v());

            Kus(i,j) += Fus(i);
            Kvs(i,j) += Fvs(i);
          }

        for (unsigned int i=0; i != n_lambda_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fulm =
              solid_context.get_elem_residual(this->_lambda_var.u());

            libMesh::DenseSubVector<libMesh::Number> & Fvlm =
              solid_context.get_elem_residual(this->_lambda_var.v());

            Kulm(i,j) += Fulm(i);
            Kvlm(i,j) += Fvlm(i);
          }

	for (unsigned int i=0; i != n_press_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fp = solid_context.get_elem_residual(this->_solid_press_var.p());

            Kp(i,j) += Fp(i);
          }
      }
  }

template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_press_derivs
  (const MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   libMesh::FEMContext & fluid_context,
   const libMesh::Real delta,
   libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
   libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
   libMesh::DenseSubVector<libMesh::Number> & press_coeff,
   libMesh::DenseSubMatrix<libMesh::Number> & Kus,
   libMesh::DenseSubMatrix<libMesh::Number> & Kvs)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_press_dofs = solid_context.get_dof_indices(this->_solid_press_var.p()).size();

    libmesh_assert_equal_to(Kus.m(),n_solid_dofs);
    libmesh_assert_equal_to(Kvs.m(),n_solid_dofs);

    libmesh_assert_equal_to(Kus.n(),n_press_dofs);
    libmesh_assert_equal_to(Kvs.n(),n_press_dofs);

    for( unsigned int j = 0; j < n_press_dofs; j++ )
      {
        // Finite differenced residual is sitting in the solid and fluid context elem_residuals
        // after this call
        this->finite_difference_residuals(system,quad_points,solid_context,fluid_context,delta,press_coeff(j),
                                          backwards_solid_residual,backwards_fluid_residual);

        for (unsigned int i=0; i != n_solid_dofs; i++)
          {
            libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(this->_disp_vars.u());
            libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(this->_disp_vars.v());

            Kus(i,j) += Fus(i);
            Kvs(i,j) += Fvs(i);
          }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::prepare_jacobians(unsigned int n_fluid_dofs,
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
                                                      libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm) const
  {
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

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::assemble_fluid_jacobians(MultiphysicsSystem & system,
                                                             const AssemblyContext & solid_context,
                                                             const libMesh::FEMContext & fluid_context,
                                                             unsigned int n_fluid_dofs,
                                                             unsigned int n_solid_dofs,
                                                             unsigned int n_lambda_dofs,
                                                             libMesh::DenseMatrix<libMesh::Number> & Kf_s,
                                                             libMesh::DenseMatrix<libMesh::Number> & Klm_f,
                                                             libMesh::DenseMatrix<libMesh::Number> & Kf_lm) const
  {
    std::vector<libMesh::dof_id_type> velocity_dof_indices;

    const std::vector<libMesh::dof_id_type>& uf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.u());

    const std::vector<libMesh::dof_id_type>& vf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.v());

    libmesh_assert_equal_to( uf_dof_indices.size(), n_fluid_dofs );
    libmesh_assert_equal_to( vf_dof_indices.size(), n_fluid_dofs );

    velocity_dof_indices.insert( velocity_dof_indices.end(),
                                 uf_dof_indices.begin(),
                                 uf_dof_indices.end() );

    velocity_dof_indices.insert( velocity_dof_indices.end(),
                                 vf_dof_indices.begin(),
                                 vf_dof_indices.end() );

    //Build up solid dof indices
    std::vector<libMesh::dof_id_type> solid_dof_indices;

    const std::vector<libMesh::dof_id_type>& us_dof_indices =
      solid_context.get_dof_indices(this->_disp_vars.u());

    const std::vector<libMesh::dof_id_type>& vs_dof_indices =
      solid_context.get_dof_indices(this->_disp_vars.v());

    libmesh_assert_equal_to( us_dof_indices.size(), n_solid_dofs );
    libmesh_assert_equal_to( vs_dof_indices.size(), n_solid_dofs );

    solid_dof_indices.insert( solid_dof_indices.end(),
                                 us_dof_indices.begin(),
                                 us_dof_indices.end() );

    solid_dof_indices.insert( solid_dof_indices.end(),
                                 vs_dof_indices.begin(),
                                 vs_dof_indices.end() );

    //Build up lambda dof indices
    std::vector<libMesh::dof_id_type> lambda_dof_indices;

    const std::vector<libMesh::dof_id_type>& ulm_dof_indices =
      solid_context.get_dof_indices(this->_lambda_var.u());

    const std::vector<libMesh::dof_id_type>& vlm_dof_indices =
      solid_context.get_dof_indices(this->_lambda_var.v());

    libmesh_assert_equal_to( ulm_dof_indices.size(), n_lambda_dofs );
    libmesh_assert_equal_to( vlm_dof_indices.size(), n_lambda_dofs );

    lambda_dof_indices.insert( lambda_dof_indices.end(),
                                 ulm_dof_indices.begin(),
                                 ulm_dof_indices.end() );

    lambda_dof_indices.insert( lambda_dof_indices.end(),
                                 vlm_dof_indices.begin(),
                                 vlm_dof_indices.end() );

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

  }

  template<typename SolidMech>
  libMesh::Point ImmersedBoundary<SolidMech>::compute_displaced_point( const MultiphysicsSystem & system,
                                                                       AssemblyContext & solid_context,
                                                                       const libMesh::Point & x_qp,
                                                                       const unsigned int qp ) const
  {
    // Compute position of displaced quadrature point at the previous time step
    libMesh::DenseVector<libMesh::Number> old_elem_solution(solid_context.get_elem_solution().size());
    libMesh::DenseVector<libMesh::Number> elem_solution_copy(solid_context.get_elem_solution().size());

    // This populates the old_elem_solution vector for the solid element
    // using the solution values from the prevous timestep
    solid_context.get_old_elem_solution(system,old_elem_solution);

    // Put in the old_elem_solution so we use the previous timestep values
    elem_solution_copy = solid_context.get_elem_solution();
    solid_context.get_elem_solution() = old_elem_solution;

    libMesh::Number u, v;
    solid_context.interior_value(this->_disp_vars.u(), qp, u );
    solid_context.interior_value(this->_disp_vars.v(), qp, v );

    // Copy back
    solid_context.get_elem_solution() = elem_solution_copy;

    return libMesh::Point(x_qp(0) + u, x_qp(1) + v);
  }

  template<typename SolidMech>
  const libMesh::Elem * ImmersedBoundary<SolidMech>::get_fluid_elem( const MultiphysicsSystem & system,
                                                                     AssemblyContext & solid_context,
                                                                     const libMesh::Point & x_qp,
                                                                     const unsigned int qp ) const
  {
    libMesh::Point x = this->compute_displaced_point(system,solid_context,x_qp,qp);

    // Now find the fluid element that contains the point x
    return (*_point_locator)( x, &_fluid_subdomain_set );
  }

 template<typename SolidMech>
 void ImmersedBoundary<SolidMech>::eval_deform_gradient( const libMesh::Gradient & grad_u,
							  const libMesh::Gradient & grad_v,
							  libMesh::TensorValue<libMesh::Real> & F )
  {

    F(0,0) = grad_u(0);
    F(0,1) = grad_u(1);
    //   grad_U(0,2) = grad_u(2);
    F(1,0) = grad_v(0);
    F(1,1) = grad_v(1);
    //    grad_U(1,2) = grad_v(2);

    F(0,0) += 1;
    F(1,1) += 1;
    F(2,2) = 1;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::compute_prev_timestep_deform_gradient
  (AssemblyContext & solid_context,
   const unsigned int qp,
   libMesh::TensorValue<libMesh::Real> & Fold)
  {
    MultiphysicsSystem & system = solid_context.get_multiphysics_system();

    libMesh::DenseVector<libMesh::Number> old_elem_solution(solid_context.get_elem_solution().size());
    libMesh::DenseVector<libMesh::Number> elem_solution_copy(solid_context.get_elem_solution().size());

    // This populates the old_elem_solution vector for the solid element
    // using the solution values from the prevous timestep
    solid_context.get_old_elem_solution(system,old_elem_solution);

    // Put in the old_elem_solution so we use the previous timestep values
    elem_solution_copy = solid_context.get_elem_solution();
    solid_context.get_elem_solution() = old_elem_solution;

    libMesh::Gradient grad_u, grad_v;
    solid_context.interior_gradient(this->_disp_vars.u(), qp, grad_u);
    solid_context.interior_gradient(this->_disp_vars.v(), qp, grad_v);

    this->eval_deform_gradient(grad_u,grad_v,Fold);

    // Copy back
    solid_context.get_elem_solution() = elem_solution_copy;
  }

 template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::eval_deform_grad_rate( const libMesh::Gradient & grad_udot,
							   const libMesh::Gradient & grad_vdot,
							   libMesh::TensorValue<libMesh::Real> & Fdot )
  {

    Fdot(0,0) = grad_udot(0);
    Fdot(0,1) = grad_udot(1);
    Fdot(1,0) = grad_vdot(0);
    Fdot(1,1) = grad_vdot(1);

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::eval_first_Piola( const libMesh::Gradient & grad_u,
                                                      const libMesh::Gradient & grad_v,
                                                      libMesh::TensorValue<libMesh::Real> & P )
  {
    /*
    libMesh::TensorValue<libMesh::Real> grad_U;
    grad_U(0,0) = grad_u(0);
    grad_U(0,1) = grad_u(1);
    grad_U(0,2) = grad_u(2);
    grad_U(1,0) = grad_v(0);
    grad_U(1,1) = grad_v(1);
    grad_U(1,2) = grad_v(2);

    libMesh::TensorValue<libMesh::Real> F(grad_U);
    F(0,0) += 1;
    F(1,1) += 1;
    F(2,2) = 1;
    */

    libMesh::TensorValue<libMesh::Real> F;
    this->eval_deform_gradient(grad_u,grad_v,F);

    libMesh::Tensor Finv( F.inverse() );

    libMesh::Number lndetF = std::log(F.det());

    libMesh::Real nus = 0.4;
    libMesh::Real mus = 5;
    libMesh::Real kappas = 2*mus*(1+nus)/(3*(1-2*nus));

    P = mus*F + kappas*lndetF*(Finv.transpose());

    /*
    // We need to use F^T a few times so just cache it.
    libMesh::TensorValue<libMesh::Real> Ftrans = F.transpose();

    libMesh::TensorValue<libMesh::Real> C(Ftrans*F);

    //libMesh::TensorValue<libMesh::Real> Cinv = C.inverse();

    //libMesh::Real I1 = C.tr();

    libMesh::TensorValue<libMesh::Real> I(1,0,0,0,1,0,0,0,1);
    //libMesh::TensorValue<libMesh::Real> S( I + (I1*I-C) );
    //S *= 2.0;

    //libMesh::TensorValue<libMesh::Real> E(grad_U + grad_U.transpose());
    libMesh::TensorValue<libMesh::Real> E(Ftrans*F);
    E(0,0) -= 1;
    E(1,1) -= 1;
    E *= 0.5;

    libMesh::Real Em = 100000;
    libMesh::Real nu = 0.3;
    libMesh::Real lambda = Em*nu/(1+nu)*(1-2*nu);
    libMesh::Real mu = Em/(2*(1+nu));

    libMesh::Real trE = E.tr();
    libMesh::TensorValue<libMesh::Real> S(2.0*E*mu);
    S(0,0) += lambda*trE;
    S(1,1) += lambda*trE;

    P = F*S;
    */
    // The F^T comes from needing the derivative of the fluid
    // shape function w.r.t. solid coordinates
    //tau = P*Ftrans;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::setup_coupling_matrix( const VelocityVariable & flow_vars,
                                                           const DisplacementVariable & disp_vars,
                                                           const LagrangeMultVectorVariable & lambda_vars,
                                                           const PressureFEVariable & press_var,
                                                           libMesh::CouplingMatrix & coupling_matrix )
  {
    libmesh_assert_equal_to(flow_vars.dim(),disp_vars.dim());
    libmesh_assert_equal_to(flow_vars.dim(),lambda_vars.dim());
    libmesh_assert_greater_equal( flow_vars.dim(), 2 );
    libmesh_assert_greater_equal( disp_vars.dim(), 2 );
    libmesh_assert_greater_equal( lambda_vars.dim(), 2 );

    // All the default coupling will be present, so we just need to indicate
    // the extra coupling from the ImmersedBoundary terms.

    // Fluid-lambda
    this->fully_coupled_vars(flow_vars,lambda_vars,coupling_matrix);

    // Fluid-solid
    this->fully_coupled_vars(flow_vars,disp_vars,coupling_matrix);

    // lambda-fluid
    this->fully_coupled_vars(lambda_vars,flow_vars,coupling_matrix);


    // We need this for algebraic ghosting. Currently, when we build the FEMContext
    // for the fluid element, it will build all variables so we need to ghost the
    // pressure even though we don't ever use it. If we someday update the FEMContext
    // to only build for a subset of variables, we should be able to remove this.
    coupling_matrix(disp_vars.u(), press_var.var()) = true;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::diagonally_coupled_vars( const MultcomponentVectorVariable & var1,
                                                             const MultcomponentVectorVariable & var2,
                                                             libMesh::CouplingMatrix & coupling_matrix )
  {
    coupling_matrix( var1.u(), var2.u() ) = true;
    coupling_matrix( var1.v(), var2.v() ) = true;
    if( var1.dim() == 3 )
      coupling_matrix( var1.w(), var2.w() ) = true;
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::fully_coupled_vars( const MultcomponentVectorVariable & var1,
                                                        const MultcomponentVectorVariable & var2,
                                                        libMesh::CouplingMatrix & coupling_matrix )
  {
    coupling_matrix( var1.u(), var2.u() ) = true;
    coupling_matrix( var1.u(), var2.v() ) = true;
    coupling_matrix( var1.v(), var2.u() ) = true;
    coupling_matrix( var1.v(), var2.v() ) = true;
    if( var1.dim() == 3 )
      {
        coupling_matrix( var1.u(), var2.w() ) = true;
        coupling_matrix( var1.v(), var2.w() ) = true;
        coupling_matrix( var1.w(), var2.u() ) = true;
        coupling_matrix( var1.w(), var2.v() ) = true;
        coupling_matrix( var1.w(), var2.w() ) = true;
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::print_coupling_matrix( const libMesh::CouplingMatrix & coupling_matrix,
                                                           const unsigned int n )
  {
    libMesh::out << std::endl << "Coupling Matrix:" << std::endl;

    for( unsigned int i = 0; i < n; i++ )
      {
        for( unsigned int j = 0; j < n; j++ )
          libMesh::out << coupling_matrix(i,j) << " ";

        libMesh::out << std::endl;
      }

    libMesh::out << std::endl << std::endl;
  }


  //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

} // namespace GRINS
