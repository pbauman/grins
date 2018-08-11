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
      _solid_mech(std::move(solid_mech_ptr)),
      _fluid_mechanics(input("Physics/ImmersedBoundary/fluid_mechanics","DIE!")),
      _solid_mechanics(input("Physics/ImmersedBoundary/solid_mechanics","DIE!"))
  {
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
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::preassembly( MultiphysicsSystem & system )
  {
    // We need to rebuild the overlap map each time because the position
    // of the solid can change in between each Newton step.
    _fluid_solid_overlap.reset( new OverlappingFluidSolidMap(system,
                                                             (*_point_locator),
                                                             _solid_subdomain_set,
                                                             _fluid_subdomain_set,
                                                             _disp_vars) );

    // Since we need to rebuild the overlap, we need to rebuild the sparsity too
    _ibm_sparsity.reset( new ImmersedBoundaryAugmentedSparsity(system,
                                                               _disp_vars,
                                                               _flow_vars,
                                                               *_fluid_solid_overlap) );
    libMesh::DofMap & dof_map = system.get_dof_map();
    dof_map.clear_sparsity();
    dof_map.attach_extra_sparsity_object(*_ibm_sparsity);
    dof_map.compute_sparsity(system.get_mesh());

    // New to reinit the matrix since we changed the sparsity pattern
    libMesh::SparseMatrix<libMesh::Number> & matrix = system.get_matrix("System Matrix");
    libmesh_assert(dof_map.is_attached(matrix));
    matrix.init();
    matrix.zero();

    libMesh::PetscMatrix<libMesh::Number> & petsc_matrix =
      libMesh::cast_ref<libMesh::PetscMatrix<libMesh::Number> &>(matrix);

    Mat A = petsc_matrix.mat();

    PetscErrorCode ierr =0;
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRABORT(system.comm().get(), ierr);
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
    if( this->is_solid_elem( context.get_elem().subdomain_id() ) )
      {
        //this->assemble_accel_term( compute_jacobian, context );
        this->assemble_solid_var_residual_contributions(compute_jacobian,context);
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::assemble_accel_term( bool compute_jacobian, AssemblyContext & context )
  {
    // For clarity
    AssemblyContext & solid_context = context;

    MultiphysicsSystem & system = context.get_multiphysics_system();

    const unsigned int n_u_dofs = solid_context.get_dof_indices(_disp_vars.u()).size();

    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();
    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);
    unsigned int v_dot_var = system.get_second_order_dot_var(v_var);

    const std::vector<libMesh::Real> &JxW =
      solid_context.get_element_fe(u_var,2)->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& u_phi =
      solid_context.get_element_fe(u_var,2)->get_phi();

    // Now assemble fluid part of velocity matching term
    // into solid residual.
    libMesh::DenseSubVector<libMesh::Number> & Fud = solid_context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvd = solid_context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> * Fwd = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> & Kudud = solid_context.get_elem_jacobian(u_dot_var,u_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvdvd = solid_context.get_elem_jacobian(v_dot_var,v_dot_var);
    libMesh::DenseSubMatrix<libMesh::Number> * Kwdwd = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> & Kudu = solid_context.get_elem_jacobian(u_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvdv = solid_context.get_elem_jacobian(v_dot_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number> * Kwdw = NULL;

    unsigned int w_var = libMesh::invalid_uint;
    unsigned int w_dot_var = libMesh::invalid_uint;
    if ( this->_disp_vars.dim() == 3 )
      {
        w_var = this->_disp_vars.w();
        w_dot_var = system.get_second_order_dot_var(w_var);

        Fwd = &solid_context.get_elem_residual(w_dot_var);
        Kwdwd = &solid_context.get_elem_jacobian(w_dot_var,w_dot_var);
        Kwdw  = &solid_context.get_elem_jacobian(w_dot_var,w_var);
      }

    unsigned int n_qpoints = solid_context.get_element_qrule().n_points();

    for( unsigned int qp = 0; qp < n_qpoints; qp++ )
      {
        libMesh::Real jac = JxW[qp];

        // Time derivative of auxillary solid velocity variable
        libMesh::Real u_dot_rate, v_dot_rate;
        solid_context.interior_rate( u_dot_var, qp, u_dot_rate );
        solid_context.interior_rate( v_dot_var, qp, v_dot_rate );

        libMesh::Real w_dot_rate = 0.0;
        if( this->_disp_vars.dim() == 3 )
          solid_context.interior_rate( w_dot_var, qp, w_dot_rate );

        // Displacement acceleration
        libMesh::Real u_ddot, v_ddot;
        solid_context.interior_accel( u_var, qp, u_ddot );
        solid_context.interior_accel( v_var, qp, v_ddot );

        libMesh::Real w_ddot = 0.0;
        if( this->_disp_vars.dim() == 3 )
          solid_context.interior_accel( w_var, qp, w_ddot );

        for (unsigned int i=0; i != n_u_dofs; i++)
          {
            libMesh::Real phi_times_jac = u_phi[i][qp]*jac;

            Fud(i) += (u_dot_rate - u_ddot)*phi_times_jac;
            Fvd(i) += (v_dot_rate - v_ddot)*phi_times_jac;

            if( this->_disp_vars.dim() == 3 )
              (*Fwd)(i) += (w_dot_rate - w_ddot)*phi_times_jac;

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    libMesh::Real udterm =
                      (u_phi[j][qp]*solid_context.get_elem_solution_rate_derivative())*phi_times_jac;

                    Kudud(i,j) += udterm;
                    Kvdvd(i,j) += udterm;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kwdwd)(i,j) += udterm;

                    libMesh::Real uterm =
                      (u_phi[j][qp]*solid_context.get_elem_solution_accel_derivative())*phi_times_jac;

                    Kudu(i,j) -= uterm;
                    Kvdv(i,j) -= uterm;

                    if( this->_disp_vars.dim() == 3 )
                      (*Kwdw)(i,j) -= uterm;
                  }
              }
          }
      }
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::assemble_solid_var_residual_contributions( bool compute_jacobian,
                                                                               AssemblyContext & context )
  {
    // For clarity
    AssemblyContext & solid_context = context;

    MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_var = this->_disp_vars.u();

    // Prepare solid info needed
    const std::vector<libMesh::Point> & solid_qpoints =
      solid_context.get_element_fe(u_var,2)->get_xyz();

    const std::vector<libMesh::Real> & solid_JxW =
      solid_context.get_element_fe(u_var,2)->get_JxW();

    const std::vector<std::vector<libMesh::Real> > & solid_phi =
      solid_context.get_element_fe(u_var,2)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & solid_dphi =
      solid_context.get_element_fe(u_var,2)->get_dphi();

    // Prepare fluid info needed
    const std::vector<std::vector<libMesh::Real> > & fluid_phi =
      _fluid_context->get_element_fe(this->_flow_vars.u())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & fluid_dphi =
      _fluid_context->get_element_fe(this->_flow_vars.u())->get_dphi();


    // We need to grab the fluid elements that are overlapping with this solid elem.
    // Then, for that fluid element, extract the indices of the *solid* quadrature points
    // that are in that fluid element
    std::map<libMesh::dof_id_type,std::map<libMesh::dof_id_type,std::vector<unsigned int> > >::const_iterator
      solid_elem_map_it = _fluid_solid_overlap->solid_map().find(solid_context.get_elem().id());

    // We should've had at least one overlapping fluid element
    libmesh_assert( solid_elem_map_it != _fluid_solid_overlap->solid_map().end() );

    const std::map<libMesh::dof_id_type,std::vector<unsigned int> > &
      fluid_elem_map = solid_elem_map_it->second;

    std::vector<libMesh::Point> solid_qpoints_subset;

    std::vector<unsigned int> qps_visted;

    for( std::map<libMesh::dof_id_type,std::vector<unsigned int> >::const_iterator
           fluid_elem_map_it = fluid_elem_map.begin();
         fluid_elem_map_it != fluid_elem_map.end();
         ++fluid_elem_map_it )
      {


        // Grab the current fluid element id
        libMesh::dof_id_type fluid_elem_id = fluid_elem_map_it->first;

        // Extract out the subset of solid quadrature points we're dealing with on
        // this fluid element
        const std::vector<unsigned int> & solid_qpoint_indices = fluid_elem_map_it->second;

        for( auto qp : solid_qpoint_indices )
          qps_visted.push_back(qp);

        // Prepare the fluid context for all the things that we'll need for the current
        // fluid element
        this->prepare_fluid_context( system,
                                     fluid_elem_id,
                                     solid_context,
                                     solid_qpoint_indices,
                                     solid_qpoints,
                                     solid_qpoints_subset,
                                     *(this->_fluid_context) );

        libmesh_assert_equal_to( solid_qpoint_indices.size(), solid_qpoints_subset.size() );

        this->add_source_term_to_fluid_residual(compute_jacobian,
                                                system,
                                                *(this->_fluid_context),
                                                solid_context,
                                                solid_qpoint_indices,
                                                solid_qpoints_subset,
                                                solid_JxW,solid_dphi,fluid_dphi);

        this->add_velocity_coupling_term_to_solid_residual(compute_jacobian,
                                                           system,
                                                           *(this->_fluid_context),
                                                           solid_context,
                                                           solid_qpoint_indices,
                                                           solid_qpoints_subset,
                                                           solid_JxW,solid_phi,fluid_phi);

      } // end loop over overlapping fluid elements

    libmesh_assert_equal_to( solid_qpoints.size(), qps_visted.size() );
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::prepare_fluid_context
  ( const MultiphysicsSystem & system,
    libMesh::dof_id_type fluid_elem_id,
    const AssemblyContext & solid_context,
    const std::vector<unsigned int> & solid_qpoint_indices,
    const std::vector<libMesh::Point> & solid_qpoints,
    std::vector<libMesh::Point> & solid_qpoints_subset,
    libMesh::FEMContext & fluid_context )
  {
    solid_qpoints_subset.clear();
    solid_qpoints_subset.reserve(solid_qpoint_indices.size());

    libMesh::Real u,v,w;

    for( unsigned int i = 0; i < solid_qpoint_indices.size(); i++ )
      {
        unsigned int qp = solid_qpoint_indices[i];

        solid_context.interior_value(this->_disp_vars.u(), qp, u );
        solid_context.interior_value(this->_disp_vars.v(), qp, v );
        if(this->_disp_vars.dim() == 3)
          solid_context.interior_value(this->_disp_vars.w(), qp, w );

        const libMesh::Point & xqp = solid_qpoints[qp];

        libMesh::Point xpu_qp(xqp(0)+u, xqp(1)+v, xqp(2)+w);

        solid_qpoints_subset.push_back( xpu_qp );
      }

    // Prepare the fluid context for things we're evaluating on the fluid
    // element at the *deformed* solid quadrature point locations within
    // the fluid element.
    const libMesh::Elem * fluid_elem = system.get_mesh().elem(fluid_elem_id);

    fluid_context.pre_fe_reinit(system,fluid_elem);
    fluid_context.elem_fe_reinit(&solid_qpoints_subset);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::add_source_term_to_fluid_residual
  ( bool compute_jacobian,
    MultiphysicsSystem & system,
    libMesh::FEMContext & fluid_context,
    AssemblyContext & solid_context,
    const std::vector<unsigned int> & solid_qpoint_indices,
    const std::vector<libMesh::Point> & solid_qpoints_subset,
    const std::vector<libMesh::Real> & solid_JxW,
    const std::vector<std::vector<libMesh::RealGradient> > & solid_dphi,
    const std::vector<std::vector<libMesh::RealGradient> > & fluid_dphi )
  {
    libMesh::DenseMatrix<libMesh::Number> Kmat;

    libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kmat), Kuf_vs(Kmat), Kuf_ws(Kmat);
    libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kmat), Kvf_vs(Kmat), Kvf_ws(Kmat);
    libMesh::DenseSubMatrix<libMesh::Number> Kwf_us(Kmat), Kwf_vs(Kmat), Kwf_ws(Kmat);



    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();

    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

    libMesh::DenseSubVector<libMesh::Number> & Fuf = fluid_context.get_elem_residual(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & Fvf = fluid_context.get_elem_residual(this->_flow_vars.v());

    if( compute_jacobian)
      {
        Kmat.resize( this->_flow_vars.dim()*n_fluid_dofs, this->_disp_vars.dim()*n_solid_dofs );

        // We need to manually manage the indexing since we're working only on this particular subblock
        Kuf_us.reposition( 0, 0, n_fluid_dofs, n_solid_dofs );
        Kuf_vs.reposition( 0, n_solid_dofs, n_fluid_dofs, n_solid_dofs );
        Kvf_us.reposition( n_fluid_dofs, 0, n_fluid_dofs, n_solid_dofs );
        Kvf_vs.reposition( n_fluid_dofs, n_solid_dofs, n_fluid_dofs, n_solid_dofs );

        if( _flow_vars.dim() == 3 )
          {
            Kuf_ws.reposition( 0, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
            Kvf_ws.reposition( n_fluid_dofs, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
            Kwf_us.reposition( 2*n_fluid_dofs, 0, n_fluid_dofs, n_solid_dofs );
            Kwf_vs.reposition( 2*n_fluid_dofs, n_solid_dofs, n_fluid_dofs, n_solid_dofs );
            Kwf_ws.reposition( 2*n_fluid_dofs, 2*n_solid_dofs, n_fluid_dofs, n_solid_dofs );
          }
      }

    libMesh::DenseSubVector<libMesh::Number> & u_coeffs = solid_context.get_elem_solution(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & v_coeffs = solid_context.get_elem_solution(this->_disp_vars.v());

    libMesh::Real delta = 1.0e-8;

    for( unsigned int qp = 0; qp < solid_qpoints_subset.size(); qp++ )
      {
        // Gradients w.r.t. the master element coordinates
        /*
        _solid_mech->get_grad_disp(solid_context, sqp,
                                   grad_u, grad_v, grad_w);

        libMesh::TensorValue<libMesh::Real> blah;
        ElasticityTensor C;
        _solid_mech->get_stress_and_elasticity(solid_context,sqp,
                                               grad_u,grad_v,grad_w,blah,C);
        */

        unsigned int sqp = solid_qpoint_indices[qp];

        libMesh::Real jac = solid_JxW[sqp];

        libMesh::Gradient grad_u, grad_v;
        solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_u);
        solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_v);


        libMesh::TensorValue<libMesh::Real> tau0;
        this->eval_stress(grad_u,grad_v,tau0);

        for (unsigned int i=0; i != n_fluid_dofs; i++)
          {

            // Zero index for fluid dphi/JxW since we only requested one quad. point.
            for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
              {
                Fuf(i) -= tau0(0,alpha)*fluid_dphi[i][qp](alpha)*jac;
                Fvf(i) -= tau0(1,alpha)*fluid_dphi[i][qp](alpha)*jac;
              }

            if( compute_jacobian )
              {
                for (unsigned int j=0; j != n_solid_dofs; j++)
                  {
                    u_coeffs(j) += delta;
                    v_coeffs(j) += delta;

                    libMesh::Gradient grad_upe, grad_vpe;
                    for (unsigned int k = 0; k != n_solid_dofs; k++ )
                      {
                        grad_upe += u_coeffs(k)*solid_dphi[k][sqp];
                        grad_vpe += v_coeffs(k)*solid_dphi[k][sqp];
                      }

                    u_coeffs(j) -= 2*delta;
                    v_coeffs(j) -= 2*delta;

                    libMesh::Gradient grad_ume, grad_vme;
                    for (unsigned int k = 0; k != n_solid_dofs; k++ )
                      {
                        grad_ume += u_coeffs(k)*solid_dphi[k][sqp];
                        grad_vme += v_coeffs(k)*solid_dphi[k][sqp];
                      }

                    u_coeffs(j) += delta;
                    v_coeffs(j) += delta;


                    libMesh::TensorValue<libMesh::Real> tau_upe;
                    this->eval_stress(grad_upe,grad_v,tau_upe);

                    libMesh::TensorValue<libMesh::Real> tau_ume;
                    this->eval_stress(grad_ume,grad_v,tau_ume);

                    libMesh::TensorValue<libMesh::Real> tau_vpe;
                    this->eval_stress(grad_u,grad_vpe,tau_vpe);

                    libMesh::TensorValue<libMesh::Real> tau_vme;
                    this->eval_stress(grad_u,grad_vme,tau_vme);

                    for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
                      {
                        Kuf_us(i,j) -= (tau_upe(0,alpha)-tau_ume(0,alpha))/(2*delta)*fluid_dphi[i][qp](alpha)*jac;
                        Kvf_us(i,j) -= (tau_upe(1,alpha)-tau_ume(1,alpha))/(2*delta)*fluid_dphi[i][qp](alpha)*jac;
                        Kuf_vs(i,j) -= (tau_vpe(0,alpha)-tau_vme(0,alpha))/(2*delta)*fluid_dphi[i][qp](alpha)*jac;
                        Kvf_vs(i,j) -= (tau_vpe(1,alpha)-tau_vme(1,alpha))/(2*delta)*fluid_dphi[i][qp](alpha)*jac;
                      }

                  }

              }// compute_jacobian

          } //fluid dof loop

      } // end solid_qpoints_subset loop

    std::vector<libMesh::dof_id_type> velocity_dof_indices;
    velocity_dof_indices.resize(_flow_vars.dim()*n_fluid_dofs);

    const std::vector<libMesh::dof_id_type>& uf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.u());

    for( unsigned int i = 0; i < uf_dof_indices.size(); i++ )
      velocity_dof_indices[i] = uf_dof_indices[i];

    const std::vector<libMesh::dof_id_type>& vf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.v());

    for( unsigned int i = 0; i < vf_dof_indices.size(); i++ )
      velocity_dof_indices[i+n_fluid_dofs] = vf_dof_indices[i];

    //Build up solid dof indices
    std::vector<libMesh::dof_id_type> solid_dof_indices;
    solid_dof_indices.resize(_disp_vars.dim()*n_solid_dofs);

    const std::vector<libMesh::dof_id_type>& us_dof_indices =
      solid_context.get_dof_indices(this->_disp_vars.u());

    for( unsigned int i = 0; i < us_dof_indices.size(); i++ )
      solid_dof_indices[i] = us_dof_indices[i];

    const std::vector<libMesh::dof_id_type>& vs_dof_indices =
      solid_context.get_dof_indices(this->_disp_vars.v());

    for( unsigned int i = 0; i < vs_dof_indices.size(); i++ )
      solid_dof_indices[i+n_solid_dofs] = vs_dof_indices[i];

    // Since we manually built the fluid context, we have to manually
    // constrain and add the residuals and Jacobians.
    //! \todo  We're hardcoding to the case that the residual is always
    //  assembled and homogeneous constraints.
    if( compute_jacobian )
      {
        system.get_dof_map().constrain_element_matrix
          ( Kmat,
            velocity_dof_indices,
            solid_dof_indices, false );

        system.matrix->add_matrix( Kmat,
                                   velocity_dof_indices,
                                   solid_dof_indices );
      }

    system.get_dof_map().constrain_element_vector
      ( fluid_context.get_elem_residual(),
        fluid_context.get_dof_indices(), false );

    system.rhs->add_vector( fluid_context.get_elem_residual(),
                            fluid_context.get_dof_indices() );

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::add_velocity_coupling_term_to_solid_residual
  (bool compute_jacobian,
   MultiphysicsSystem & system,
   libMesh::FEMContext & fluid_context,
   AssemblyContext & solid_context,
   const std::vector<unsigned int> & solid_qpoint_indices,
   const std::vector<libMesh::Point> & solid_qpoints_subset,
   const std::vector<libMesh::Real> & solid_JxW,
   const std::vector<std::vector<libMesh::Real> > & solid_phi,
   const std::vector<std::vector<libMesh::Real> > & fluid_phi)
  {
    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();

    /*
    unsigned int w_var = libMesh::invalid_uint;

    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);
    unsigned int v_dot_var = system.get_second_order_dot_var(v_var);
    unsigned int w_dot_var = libMesh::invalid_uint;

    // Shift if we're not using first order time solver
    unsigned int sshift = 1;
    if( u_var != u_dot_var )
      sshift = 2;
    */

    const unsigned int n_solid_dofs = solid_context.get_dof_indices(u_var).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

    libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(u_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(v_var);
    //libMesh::DenseSubVector<libMesh::Number> * Fws = NULL;

    libMesh::DenseSubMatrix<libMesh::Number> & Kus_us = solid_context.get_elem_jacobian(u_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vs = solid_context.get_elem_jacobian(v_var,v_var);
    //libMesh::DenseSubMatrix<libMesh::Number> * Kws_ws = NULL;

    /*
    if ( this->_disp_vars.dim() == 3 )
      {
        w_var = this->_disp_vars.w();
        w_dot_var = system.get_second_order_dot_var(w_var);

        Fws = &solid_context.get_elem_residual(w_var);

        Kws_ws = &solid_context.get_elem_jacobian(w_var,w_dot_var);
      }
    */
    // Setup solid-fluid coupling Jacobian
    // The residual equation is (\dot{U} - V)
    // Thus, this depends on wheter or not we're using a first order time solver
    // If we are *NOT* using a first order time solver, then the coupling
    // is U-V.
    // If we *ARE* using a first order time solver, then the coupling
    // is the to the \dot{U} variable (added automatically by FEMSystem) and V.
    libMesh::DenseMatrix<libMesh::Number> K;
    libMesh::DenseSubMatrix<libMesh::Number> Kus_uf(K), Kvs_vf(K), Kws_wf(K);
    if( compute_jacobian)
      {
        K.resize( this->_disp_vars.dim()*n_solid_dofs, this->_flow_vars.dim()*n_fluid_dofs );

         Kus_uf.reposition( 0, 0, n_solid_dofs, n_fluid_dofs );
         Kvs_vf.reposition( n_solid_dofs, n_fluid_dofs, n_solid_dofs, n_fluid_dofs );
      }

    for( unsigned int qp = 0; qp < solid_qpoints_subset.size(); qp++ )
      {
        unsigned int solid_qp_idx = solid_qpoint_indices[qp];

        // Compute the fluid velocity at the solid element quadrature points.
        libMesh::Real Vx, Vy;

        fluid_context.interior_value(this->_flow_vars.u(), qp, Vx);
        fluid_context.interior_value(this->_flow_vars.v(), qp, Vy);

        // Compute the dipslacement velocity
        libMesh::Real udot, vdot;
        solid_context.interior_rate(u_var, solid_qp_idx, udot);
        solid_context.interior_rate(v_var, solid_qp_idx, vdot);

        libMesh::Real jac = solid_JxW[solid_qp_idx];

        for( unsigned int i = 0; i < n_solid_dofs; i++ )
          {
            libMesh::Real sphi_times_jac = solid_phi[i][solid_qp_idx]*jac;

            Fus(i) += (-udot + Vx)*sphi_times_jac;
            Fvs(i) += (-vdot + Vy)*sphi_times_jac;

            if( compute_jacobian )
              {
                // Solid-solid block
                for( unsigned int j = 0; j < n_solid_dofs; j++ )
                  {
                    libMesh::Real sphij = solid_phi[j][solid_qp_idx];

                    libMesh::Real diag_value = sphij*sphi_times_jac;
                    diag_value *= solid_context.get_elem_solution_rate_derivative();

                    Kus_us(i,j) -= diag_value;
                    Kvs_vs(i,j) -= diag_value;
                  }

                // Solid-fluid block
                for( unsigned int j = 0; j < n_fluid_dofs; j++ )
                  {
                    libMesh::Real diag_value = fluid_phi[j][qp]*sphi_times_jac*
                      solid_context.get_elem_solution_derivative();

                    Kus_uf(i,j) += diag_value;
                    Kvs_vf(i,j) += diag_value;
                  }

              } // if compute jacobian

          } // solid dof loop

      } // end solid_qpoints_subset loop

    // Prepare dof indice vectors for constraining the coupled solid-fluid Jacobian
    std::vector<libMesh::dof_id_type> solid_dof_indices;
    solid_dof_indices.resize(_disp_vars.dim()*n_solid_dofs);

    const std::vector<libMesh::dof_id_type>& us_dof_indices =
      solid_context.get_dof_indices(u_var);

    for( unsigned int i = 0; i < us_dof_indices.size(); i++ )
      solid_dof_indices[i] = us_dof_indices[i];

    const std::vector<libMesh::dof_id_type>& vs_dof_indices =
      solid_context.get_dof_indices(v_var);

    for( unsigned int i = 0; i < vs_dof_indices.size(); i++ )
      solid_dof_indices[i+n_solid_dofs] = vs_dof_indices[i];

    std::vector<libMesh::dof_id_type> velocity_dof_indices;
    velocity_dof_indices.resize(_flow_vars.dim()*n_fluid_dofs);

    const std::vector<libMesh::dof_id_type>& uf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.u());

    for( unsigned int i = 0; i < uf_dof_indices.size(); i++ )
      velocity_dof_indices[i] = uf_dof_indices[i];

    const std::vector<libMesh::dof_id_type>& vf_dof_indices =
      fluid_context.get_dof_indices(this->_flow_vars.v());

    for( unsigned int i = 0; i < vf_dof_indices.size(); i++ )
      velocity_dof_indices[i+n_fluid_dofs] = vf_dof_indices[i];

    if( compute_jacobian )
      {
        system.get_dof_map().constrain_element_matrix
          ( K,
            solid_dof_indices,
            velocity_dof_indices,
            false );

        system.matrix->add_matrix( K,
                                   solid_dof_indices,
                                   velocity_dof_indices );
      }

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::eval_stress( const libMesh::Gradient & grad_u,
                                                 const libMesh::Gradient & grad_v,
                                                 libMesh::TensorValue<libMesh::Real> & tau )
  {
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

    // We need to use F^T a few times so just cache it.
    libMesh::TensorValue<libMesh::Real> Ftrans = F.transpose();

    libMesh::TensorValue<libMesh::Real> C(Ftrans*F);

    //libMesh::TensorValue<libMesh::Real> Cinv = C.inverse();

    libMesh::Real I1 = C.tr();

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

    libMesh::TensorValue<libMesh::Real> P(F*S);

    // The F^T comes from needing the derivative of the fluid
    // shape function w.r.t. solid coordinates
    tau = P*Ftrans;
  }

    //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

} // namespace GRINS
