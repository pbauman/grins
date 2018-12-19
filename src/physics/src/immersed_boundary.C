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
      _lambda_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<LagrangeMultVectorVariable>(VariablesParsing::lagrange_mult_variable_name(input,physics_name,VariablesParsing::PHYSICS))),
      _solid_mech(std::move(solid_mech_ptr)),
      _fluid_mechanics(input("Physics/ImmersedBoundary/fluid_mechanics","DIE!")),
      _solid_mechanics(input("Physics/ImmersedBoundary/solid_mechanics","DIE!"))
  {
    _lambda_var.set_is_constraint_var(true);

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
    {
      std::unique_ptr<libMesh::DiffContext> raw_context = system.build_context();
      libMesh::FEMContext * context = libMesh::cast_ptr<libMesh::FEMContext *>(raw_context.release());
      point_fluid_context.reset(context);
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
    // For clarity
    AssemblyContext & solid_context = context;

    MultiphysicsSystem & system = context.get_multiphysics_system();

    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();

    unsigned int lambda_x = this->_lambda_var.u();
    unsigned int lambda_y = this->_lambda_var.v();

    // Prepare solid info needed
    const std::vector<libMesh::Point> & solid_qpoints =
      solid_context.get_element_fe(u_var,2)->get_xyz();

    const std::vector<libMesh::Real> & solid_JxW =
      solid_context.get_element_fe(u_var,2)->get_JxW();

    libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(u_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(v_var);

    // For computing numerical Jacobians
    libMesh::DenseSubVector<libMesh::Number> Fusp(Fus), Fusm(Fus);
    libMesh::DenseSubVector<libMesh::Number> Fvsp(Fvs), Fvsm(Fvs);

    libMesh::DenseSubMatrix<libMesh::Number> & Kus_us = solid_context.get_elem_jacobian(u_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_us = solid_context.get_elem_jacobian(v_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kus_vs = solid_context.get_elem_jacobian(u_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vs = solid_context.get_elem_jacobian(v_var,v_var);

    libMesh::DenseSubVector<libMesh::Number> & Fulm = solid_context.get_elem_residual(lambda_x);
    libMesh::DenseSubVector<libMesh::Number> & Fvlm = solid_context.get_elem_residual(lambda_y);

    // For computing numerical Jacobians
    libMesh::DenseSubVector<libMesh::Number> Fulmp(Fulm), Fulmm(Fulm);
    libMesh::DenseSubVector<libMesh::Number> Fvlmp(Fvlm), Fvlmm(Fvlm);

    libMesh::DenseSubMatrix<libMesh::Number> & Kulm_us = solid_context.get_elem_jacobian(lambda_x,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_us = solid_context.get_elem_jacobian(lambda_y,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kulm_vs = solid_context.get_elem_jacobian(lambda_x,v_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vs = solid_context.get_elem_jacobian(lambda_y,v_var);

    libMesh::DenseSubMatrix<libMesh::Number> & Kus_ulm = solid_context.get_elem_jacobian(u_var,lambda_x);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_ulm = solid_context.get_elem_jacobian(v_var,lambda_x);
    libMesh::DenseSubMatrix<libMesh::Number> & Kus_vlm = solid_context.get_elem_jacobian(u_var,lambda_y);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vlm = solid_context.get_elem_jacobian(v_var,lambda_y);

    libMesh::DenseMatrix<libMesh::Number> Kf_s;
    libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kf_s), Kuf_vs(Kf_s);
    libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kf_s), Kvf_vs(Kf_s);

    libMesh::DenseMatrix<libMesh::Number> Klm_f;
    libMesh::DenseSubMatrix<libMesh::Number> Kulm_uf(Klm_f), Kulm_vf(Klm_f);
    libMesh::DenseSubMatrix<libMesh::Number> Kvlm_uf(Klm_f), Kvlm_vf(Klm_f);

    libMesh::DenseMatrix<libMesh::Number> Kf_lm;
    libMesh::DenseSubMatrix<libMesh::Number> Kuf_ulm(Kf_lm), Kuf_vlm(Kf_lm);
    libMesh::DenseSubMatrix<libMesh::Number> Kvf_ulm(Kf_lm), Kvf_vlm(Kf_lm);

    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

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

	libMesh::Real delta = 1.0e-8;

        // Grab the current fluid element id
        libMesh::dof_id_type fluid_elem_id = fluid_elem_map_it->first;

        // Extract out the subset of solid quadrature points we're dealing with on
        // this fluid element
        const std::vector<unsigned int> & solid_qpoint_indices = fluid_elem_map_it->second;

        for( auto qp : solid_qpoint_indices )
          qps_visted.push_back(qp);

        // Prepare the fluid context for all the things that we'll need for the current
        // fluid element
        this->prepare_fluid_context_batch( system,
                                           fluid_elem_id,
                                           solid_context,
                                           solid_qpoint_indices,
                                           solid_qpoints,
                                           solid_qpoints_subset,
                                           *(this->_fluid_context) );

        libmesh_assert_equal_to( solid_qpoint_indices.size(), solid_qpoints_subset.size() );

	libMesh::DenseSubVector<libMesh::Number> & Fuf =
          (this->_fluid_context)->get_elem_residual(this->_flow_vars.u());

	libMesh::DenseSubVector<libMesh::Number> & Fvf =
          (this->_fluid_context)->get_elem_residual(this->_flow_vars.v());

        // For computing numerical Jacobians
        libMesh::DenseSubVector<libMesh::Number> Fufp(Fuf), Fufm(Fuf);
        libMesh::DenseSubVector<libMesh::Number> Fvfp(Fvf), Fvfm(Fvf);

	unsigned int n_fluid_dofs = (this->_fluid_context)->get_dof_indices(this->_flow_vars.u()).size();

	if ( compute_jacobian )
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

	for( unsigned int qp = 0; qp < solid_qpoints_subset.size(); qp++ )
	  {
	    unsigned int sqp = solid_qpoint_indices[qp];

	    libMesh::Real jac = solid_JxW[sqp];

	    this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,*(this->point_fluid_context));

	    this->fluid_residual_contribution(compute_jacobian,system,
					      *(this->point_fluid_context),fluid_elem_id,
					      solid_context,solid_qpoints,sqp,
					      jac,delta,Fuf,Fvf,
					      Kuf_us,Kuf_vs,Kvf_us,Kvf_vs,
					      Kuf_ulm,Kuf_vlm,Kvf_ulm,Kvf_vlm);

	    this->solid_residual_contribution(compute_jacobian,
					      solid_context,sqp,
					      jac,delta,Fus,Fvs,
					      Kus_us,Kvs_us,Kus_vs,Kvs_vs,
					      Kus_ulm,Kvs_ulm,Kus_vlm,Kvs_vlm);

	    this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,*(this->point_fluid_context));

	    this->lambda_residual_contribution(compute_jacobian,system,
					       *(this->point_fluid_context),fluid_elem_id,
					       solid_context,solid_qpoints,sqp,
					       jac,delta,Fulm,Fvlm,
					       Kulm_uf,Kvlm_uf,Kulm_vf,Kvlm_vf,
					       Kulm_us,Kvlm_us,Kulm_vs,Kvlm_vs);

	  } // end solid_qpoints_subset loop

	std::vector<libMesh::dof_id_type> velocity_dof_indices;
	velocity_dof_indices.resize(_flow_vars.dim()*n_fluid_dofs);

	const std::vector<libMesh::dof_id_type>& uf_dof_indices =
	  (this->_fluid_context)->get_dof_indices(this->_flow_vars.u());

	for( unsigned int i = 0; i < uf_dof_indices.size(); i++ )
	  velocity_dof_indices[i] = uf_dof_indices[i];

	const std::vector<libMesh::dof_id_type>& vf_dof_indices =
	  (this->_fluid_context)->get_dof_indices(this->_flow_vars.v());

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

	//Build up lambda dof indices
	std::vector<libMesh::dof_id_type> lambda_dof_indices;
	lambda_dof_indices.resize(_lambda_var.dim()*n_lambda_dofs);

	const std::vector<libMesh::dof_id_type>& ulm_dof_indices =
	  solid_context.get_dof_indices(this->_lambda_var.u());

	for( unsigned int i = 0; i < ulm_dof_indices.size(); i++ )
	  lambda_dof_indices[i] = ulm_dof_indices[i];

	const std::vector<libMesh::dof_id_type>& vlm_dof_indices =
	  solid_context.get_dof_indices(this->_lambda_var.v());

	for( unsigned int i = 0; i < vlm_dof_indices.size(); i++ )
	  lambda_dof_indices[i+n_lambda_dofs] = vlm_dof_indices[i];

	// Since we manually built the fluid context, we have to manually
	// constrain and add the residuals and Jacobians.
	//! \todo  We're hardcoding to the case that the residual is always
	//  assembled and homogeneous constraints.
	if( compute_jacobian )
	  {
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

	system.get_dof_map().constrain_element_vector
	  ( (this->_fluid_context)->get_elem_residual(),
	    (this->_fluid_context)->get_dof_indices(), false );

	system.rhs->add_vector( (this->_fluid_context)->get_elem_residual(),
				(this->_fluid_context)->get_dof_indices() );

      } // end loop over overlapping fluid elements

    libmesh_assert_equal_to( solid_qpoints.size(), qps_visted.size() );
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::prepare_fluid_context_batch( const MultiphysicsSystem & system,
								 libMesh::dof_id_type fluid_elem_id,
								 const AssemblyContext & solid_context,
								 const std::vector<unsigned int> & solid_qpoint_indices,
								 const std::vector<libMesh::Point> & solid_qpoints,
								 std::vector<libMesh::Point> & solid_qpoints_subset,
								 libMesh::FEMContext & fluid_context )
  {
    solid_qpoints_subset.clear();
    solid_qpoints_subset.reserve(solid_qpoint_indices.size());

    // We compute the *physical* location of the points we need
    std::vector<libMesh::Point> solid_qpoints_subset_xyz(solid_qpoint_indices.size());

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

        solid_qpoints_subset_xyz[i] = xpu_qp;
      }

    // Prepare the fluid context for things we're evaluating on the fluid
    // element at the *deformed* solid quadrature point locations within
    // the fluid element.
    const libMesh::Elem * fluid_elem = system.get_mesh().elem(fluid_elem_id);

    fluid_context.pre_fe_reinit(system,fluid_elem);

    libMesh::FEBase * fe = fluid_context.get_element_fe(this->_flow_vars.u());
    libMesh::FEType fetype = fe->get_fe_type();

    // But we need to hand *reference* element points to the FEMContext to reinit
    unsigned int dim = 2;
    libMesh::FEInterface::inverse_map(dim,fetype,fluid_elem,solid_qpoints_subset_xyz,solid_qpoints_subset);

    fluid_context.elem_fe_reinit(&solid_qpoints_subset);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::prepare_fluid_context( const MultiphysicsSystem & system,
							   const AssemblyContext & solid_context,
							   const std::vector<libMesh::Point> & solid_qpoints,
							   unsigned int sqp, /* solid quadrature point */
							   libMesh::dof_id_type fluid_elem_id,
							   libMesh::FEMContext & fluid_context )
  {
    libMesh::Real Ux, Uy, Uz;
    solid_context.interior_value(this->_disp_vars.u(), sqp, Ux );
    solid_context.interior_value(this->_disp_vars.v(), sqp, Uy );
    if(this->_disp_vars.dim() == 3)
      solid_context.interior_value(this->_disp_vars.w(), sqp, Uz );

    const libMesh::Point & xqp = solid_qpoints[sqp];

    libMesh::Point xpu_qp( xqp(0)+Ux, xqp(1)+Uy, xqp(2)+Uz );

    const libMesh::Elem * fluid_elem = system.get_mesh().elem(fluid_elem_id);

    fluid_context.pre_fe_reinit(system,fluid_elem);

    libMesh::FEBase * fe = fluid_context.get_element_fe(this->_flow_vars.u());
    libMesh::FEType fetype = fe->get_fe_type();

    // But we need to hand *reference* element points to the FEMContext to reinit
    unsigned int dim = 2;
    libMesh::Point x_ref =
      libMesh::FEInterface::inverse_map(dim,fetype,fluid_elem,xpu_qp);

    std::vector<libMesh::Point> x_ref_vec(1,x_ref);

    fluid_context.elem_fe_reinit(&x_ref_vec);
  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::fluid_residual_contribution( bool compute_jacobian, MultiphysicsSystem & system,
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
								 libMesh::DenseSubMatrix<libMesh::Number> & Kvf_vlm)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    libMesh::DenseSubVector<libMesh::Number> & u_coeffs = solid_context.get_elem_solution(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & v_coeffs = solid_context.get_elem_solution(this->_disp_vars.v());

    libMesh::Real lambda_x, lambda_y;
    solid_context.interior_value(this->_lambda_var.u(), sqp, lambda_x);
    solid_context.interior_value(this->_lambda_var.v(), sqp, lambda_y);

    libMesh::DenseSubVector<libMesh::Number> & lambda_xcoeff = solid_context.get_elem_solution(this->_lambda_var.u());
    libMesh::DenseSubVector<libMesh::Number> & lambda_ycoeff = solid_context.get_elem_solution(this->_lambda_var.v());
    /*
    libMesh::Gradient grad_u, grad_v;
    solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_u);
    solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_v);

    libMesh::Gradient grad_lambda_x, grad_lambda_y;
    solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lambda_x);
    solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lambda_y);

    libMesh::TensorValue<libMesh::Real> F;
    this->eval_deform_gradient(grad_u,grad_v,F);
    */

    const std::vector<std::vector<libMesh::Real> > fluid_phi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();
    //const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi =
    //  fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

    //libMesh::TensorValue<libMesh::Real> fdphi_times_F;

    for (unsigned int i=0; i != n_fluid_dofs; i++)
      {
	libmesh_assert_equal_to( fluid_phi[i].size(), 1 );

	/*
	//Computing fdphi_times_F
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    for( unsigned int beta = 0; beta < 2; beta++ )
	      {
		fdphi_times_F(0,alpha) += fluid_dphi[i][0](beta)*F(beta,alpha);
		fdphi_times_F(1,alpha) += fluid_dphi[i][0](beta)*F(beta,alpha);
	      }
   	  }
	*/
	// Zero index for fluid dphi/JxW since we only requested one quad. point.
	/*
	// H1 Norm
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fuf(i) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha) + lambda_x*fluid_phi[i][0])*jac;
	    Fvf(i) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha) + lambda_y*fluid_phi[i][0])*jac;
	  }
	*/
	// L2 Norm
	Fuf(i) -= lambda_x*fluid_phi[i][0]*jac;
	Fvf(i) -= lambda_y*fluid_phi[i][0]*jac;

	if( compute_jacobian )
	  {
	    // fluid-solid block
	    for (unsigned int j=0; j != n_solid_dofs; j++)
	      {

		//Lambda terms don't vary with solid so we don't need this

		/*
		// Computing the grad_lambda and lambda derivative terms w.r.t solid

		u_coeffs(j) += delta;

		libMesh::Gradient grad_lmx_upd, grad_lmy_upd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_upd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_upd);

		libMesh::Real lmx_upd, lmy_upd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_upd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_upd);

		u_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_lmx_umd, grad_lmy_umd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_umd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_umd);

		libMesh::Real lmx_umd, lmy_umd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_umd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_umd);

		u_coeffs(j) += delta;

		v_coeffs(j) += delta;

		libMesh::Gradient grad_lmx_vpd, grad_lmy_vpd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_vpd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_vpd);

		libMesh::Real lmx_vpd, lmy_vpd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_vpd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_vpd);

		v_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_lmx_vmd, grad_lmy_vmd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_vmd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_vmd);

		libMesh::Real lmx_vmd, lmy_vmd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_vmd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_vmd);

		v_coeffs(j) += delta;

		//H1 Norm

		// Finite differencing the grad_lambda terms w.r.t solid
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= (((grad_lmx_upd(alpha)-grad_lmx_umd(alpha))/(2*delta))*fdphi_times_F(0,alpha)
				    + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_us(i,j) -= (((grad_lmy_upd(alpha)-grad_lmy_umd(alpha))/(2*delta))*fdphi_times_F(1,alpha)
				    + lambda_y*fluid_phi[i][0])*jac;

		    Kuf_vs(i,j) -= (((grad_lmx_vpd(alpha)-grad_lmx_vmd(alpha))/(2*delta))*fdphi_times_F(0,alpha)
				    + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_vs(i,j) -= (((grad_lmy_vpd(alpha)-grad_lmy_vmd(alpha))/(2*delta))*fdphi_times_F(1,alpha)
				    + lambda_y*fluid_phi[i][0])*jac;
		  }

		// Finite differencing the lambda terms w.r.t solid
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + ((lmx_upd-lmx_umd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kvf_us(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + ((lmy_upd-lmy_umd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kuf_vs(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + ((lmx_vpd-lmx_vmd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kvf_vs(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + ((lmy_vpd-lmy_vmd)/(2*delta))*fluid_phi[i][0])*jac;
		  }


		//L2 Norm

		// Finite differencing the lambda terms w.r.t solid
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= ((lmx_upd-lmx_umd)/(2*delta))*fluid_phi[i][0]*jac;

		    Kvf_us(i,j) -= ((lmy_upd-lmy_umd)/(2*delta))*fluid_phi[i][0]*jac;

		    Kuf_vs(i,j) -= ((lmx_vpd-lmx_vmd)/(2*delta))*fluid_phi[i][0]*jac;

		    Kvf_vs(i,j) -= ((lmy_vpd-lmy_vmd)/(2*delta))*fluid_phi[i][0]*jac;
		  }
		*/

		/*

		// Finite differencing F terms

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::Gradient grad_upd, grad_vpd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_upd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vpd);

		u_coeffs(j) -= 2*delta;
		v_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_umd, grad_vmd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_umd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vmd);

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::TensorValue<libMesh::Real> F_upd;
		this->eval_deform_gradient(grad_upd,grad_v,F_upd);

		libMesh::TensorValue<libMesh::Real> F_umd;
		this->eval_deform_gradient(grad_umd,grad_v,F_umd);

		libMesh::TensorValue<libMesh::Real> F_vpd;
		this->eval_deform_gradient(grad_u,grad_vpd,F_vpd);

		libMesh::TensorValue<libMesh::Real> F_vmd;
		this->eval_deform_gradient(grad_u,grad_vmd,F_vmd);

		libMesh::TensorValue<libMesh::Real> ftF_us, ftF_vs;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    for( unsigned int beta = 0; beta < 2; beta++ )
		      {
			ftF_us(0,alpha) += fluid_dphi[i][0](beta)*((F_upd(beta,alpha)-F_umd(beta,alpha))/(2*delta));
			ftF_us(1,alpha) += fluid_dphi[i][0](beta)*((F_upd(beta,alpha)-F_umd(beta,alpha))/(2*delta));
			ftF_vs(0,alpha) += fluid_dphi[i][0](beta)*((F_vpd(beta,alpha)-F_vmd(beta,alpha))/(2*delta));
			ftF_vs(1,alpha) += fluid_dphi[i][0](beta)*((F_vpd(beta,alpha)-F_vmd(beta,alpha))/(2*delta));
		      }
		  }

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= (grad_lambda_x(alpha)*ftF_us(0,alpha) + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_us(i,j) -= (grad_lambda_y(alpha)*ftF_us(1,alpha) + lambda_y*fluid_phi[i][0])*jac;

		    Kuf_vs(i,j) -= (grad_lambda_x(alpha)*ftF_vs(0,alpha) + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_vs(i,j) -= (grad_lambda_y(alpha)*ftF_vs(1,alpha) + lambda_y*fluid_phi[i][0])*jac;
		  }
		*/

		// Compute fluid_phi and fluid_dphi derivative terms

		u_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		const std::vector<std::vector<libMesh::Real> > fluid_phi_upd =
		  fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi_upd =
		//  fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

		u_coeffs(j) -= 2*delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		const std::vector<std::vector<libMesh::Real> > fluid_phi_umd =
		  fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi_umd =
		//  fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

		u_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		v_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		const std::vector<std::vector<libMesh::Real> > fluid_phi_vpd =
		  fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi_vpd =
		//  fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

		v_coeffs(j) -= 2*delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		const std::vector<std::vector<libMesh::Real> > fluid_phi_vmd =
		  fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > fluid_dphi_vmd =
		//  fluid_context.get_element_fe(this->_flow_vars.u())->get_dphi();

		v_coeffs(j) += delta;

		//H1 Norm
		/*
		//Finite differencing the fluid_dphi terms
		libMesh::TensorValue<libMesh::Real> fdphi_F_us, fdphi_F_vs;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    for( unsigned int beta = 0; beta < 2; beta++ )
		      {
			fdphi_F_us(0,alpha) += ((fluid_dphi_upd[i][0](beta)-fluid_dphi_umd[i][0](beta))/(2*delta))*F(beta,alpha);
			fdphi_F_us(1,alpha) += ((fluid_dphi_upd[i][0](beta)-fluid_dphi_umd[i][0](beta))/(2*delta))*F(beta,alpha);
			fdphi_F_vs(0,alpha) += ((fluid_dphi_vpd[i][0](beta)-fluid_dphi_vmd[i][0](beta))/(2*delta))*F(beta,alpha);
			fdphi_F_vs(1,alpha) += ((fluid_dphi_vpd[i][0](beta)-fluid_dphi_vmd[i][0](beta))/(2*delta))*F(beta,alpha);
		      }
		  }


		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= (grad_lambda_x(alpha)*fdphi_F_us(0,alpha) + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_us(i,j) -= (grad_lambda_y(alpha)*fdphi_F_us(1,alpha) + lambda_y*fluid_phi[i][0])*jac;

		    Kuf_vs(i,j) -= (grad_lambda_x(alpha)*fdphi_F_vs(0,alpha) + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_vs(i,j) -= (grad_lambda_y(alpha)*fdphi_F_vs(1,alpha) + lambda_y*fluid_phi[i][0])*jac;
		  }

		//Finite differencing the fluid_phi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_us(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + lambda_x*((fluid_phi_upd[i][0]-fluid_phi_umd[i][0])/(2*delta)))*jac;

		    Kvf_us(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + lambda_y*((fluid_phi_upd[i][0]-fluid_phi_umd[i][0])/(2*delta)))*jac;

		    Kuf_vs(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + lambda_x*((fluid_phi_vpd[i][0]-fluid_phi_vmd[i][0])/(2*delta)))*jac;

		    Kvf_vs(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + lambda_y*((fluid_phi_vpd[i][0]-fluid_phi_vmd[i][0])/(2*delta)))*jac;
		  }
		*/
		//L2 Norm

		//Finite differencing the fluid_phi terms

		Kuf_us(i,j) -= lambda_x*((fluid_phi_upd[i][0]-fluid_phi_umd[i][0])/(2*delta))*jac;

		Kvf_us(i,j) -= lambda_y*((fluid_phi_upd[i][0]-fluid_phi_umd[i][0])/(2*delta))*jac;

		Kuf_vs(i,j) -= lambda_x*((fluid_phi_vpd[i][0]-fluid_phi_vmd[i][0])/(2*delta))*jac;

		Kvf_vs(i,j) -= lambda_y*((fluid_phi_vpd[i][0]-fluid_phi_vmd[i][0])/(2*delta))*jac;


		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

	      } //solid dof loop

	    // fluid-lambda block
	    for (unsigned int j=0; j != n_lambda_dofs; j++)
	      {

		// Computing the grad_lambda and lambda derivative terms w.r.t lambda

		lambda_xcoeff(j) += delta;
		/*
		libMesh::Gradient grad_lmx_xpd, grad_lmy_xpd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_xpd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_xpd);
		*/
		libMesh::Real lmx_xpd, lmy_xpd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_xpd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_xpd);

		lambda_xcoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_lmx_xmd, grad_lmy_xmd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_xmd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_xmd);
		*/
		libMesh::Real lmx_xmd, lmy_xmd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_xmd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_xmd);

		lambda_xcoeff(j) += delta;

		lambda_ycoeff(j) += delta;
		/*
		libMesh::Gradient grad_lmx_ypd, grad_lmy_ypd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_ypd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_ypd);
		*/
		libMesh::Real lmx_ypd, lmy_ypd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_ypd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_ypd);

		lambda_ycoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_lmx_ymd, grad_lmy_ymd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_ymd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_ymd);
		*/
		libMesh::Real lmx_ymd, lmy_ymd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_ymd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_ymd);

		lambda_ycoeff(j) += delta;

		//H1 Norm
		/*
		//Finite differencing the grad_lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_ulm(i,j) -= (((grad_lmx_xpd(alpha)-grad_lmx_xmd(alpha))/(2*delta))*fdphi_times_F(0,alpha)
				    + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_ulm(i,j) -= (((grad_lmy_xpd(alpha)-grad_lmy_xmd(alpha))/(2*delta))*fdphi_times_F(1,alpha)
				    + lambda_y*fluid_phi[i][0])*jac;

		    Kuf_vlm(i,j) -= (((grad_lmx_ypd(alpha)-grad_lmx_ymd(alpha))/(2*delta))*fdphi_times_F(0,alpha)
				    + lambda_x*fluid_phi[i][0])*jac;

		    Kvf_vlm(i,j) -= (((grad_lmy_ypd(alpha)-grad_lmy_ymd(alpha))/(2*delta))*fdphi_times_F(1,alpha)
				    + lambda_y*fluid_phi[i][0])*jac;
		  }

		//Finite differencing the lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kuf_ulm(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + ((lmx_xpd-lmx_xmd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kvf_ulm(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + ((lmy_xpd-lmy_xmd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kuf_vlm(i,j) -= (grad_lambda_x(alpha)*fdphi_times_F(0,alpha)
				    + ((lmx_ypd-lmx_ymd)/(2*delta))*fluid_phi[i][0])*jac;

		    Kvf_vlm(i,j) -= (grad_lambda_y(alpha)*fdphi_times_F(1,alpha)
				    + ((lmy_ypd-lmy_ymd)/(2*delta))*fluid_phi[i][0])*jac;
		  }
		*/

		//L2 Norm
		//Finite differencing the lambda terms

		Kuf_ulm(i,j) -= ((lmx_xpd-lmx_xmd)/(2*delta))*fluid_phi[i][0]*jac;

		Kvf_ulm(i,j) -= ((lmy_xpd-lmy_xmd)/(2*delta))*fluid_phi[i][0]*jac;

		Kuf_vlm(i,j) -= ((lmx_ypd-lmx_ymd)/(2*delta))*fluid_phi[i][0]*jac;

		Kvf_vlm(i,j) -= ((lmy_ypd-lmy_ymd)/(2*delta))*fluid_phi[i][0]*jac;


	      } //lambda dof loop

	  }// compute_jacobian

      } //fluid dof loop

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::solid_residual_contribution( bool compute_jacobian,
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
								 libMesh::DenseSubMatrix<libMesh::Number> & Kvs_vlm)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    const std::vector<std::vector<libMesh::Real> > solid_phi =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > solid_dphi =
      solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

    libMesh::DenseSubVector<libMesh::Number> & u_coeffs = solid_context.get_elem_solution(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & v_coeffs = solid_context.get_elem_solution(this->_disp_vars.v());

    libMesh::Real lambda_x, lambda_y;
    solid_context.interior_value(this->_lambda_var.u(), sqp, lambda_x);
    solid_context.interior_value(this->_lambda_var.v(), sqp, lambda_y);

    libMesh::DenseSubVector<libMesh::Number> & lambda_xcoeff = solid_context.get_elem_solution(this->_lambda_var.u());
    libMesh::DenseSubVector<libMesh::Number> & lambda_ycoeff = solid_context.get_elem_solution(this->_lambda_var.v());

    libMesh::Gradient grad_u, grad_v;
    solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_u);
    solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_v);
    /*
    libMesh::Gradient grad_lambda_x, grad_lambda_y;
    solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lambda_x);
    solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lambda_y);
    */
    libMesh::TensorValue<libMesh::Real> P;
    this->eval_first_Piola(grad_u,grad_v,P);

    for (unsigned int i=0; i != n_solid_dofs; i++)
      {
	/*
	//H1 Norm
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
		       + lambda_x*solid_phi[i][sqp])*jac;

	    Fvs(i) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
		       + lambda_y*solid_phi[i][sqp])*jac;
	  }
	*/
	//L2 Norm
	Fus(i) += lambda_x*solid_phi[i][sqp]*jac;

	Fvs(i) += lambda_y*solid_phi[i][sqp]*jac;

	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fus(i) -= P(0,alpha)*solid_dphi[i][sqp](alpha)*jac;

	    Fvs(i) -= P(1,alpha)*solid_dphi[i][sqp](alpha)*jac;
	  }

	if( compute_jacobian )
	  {
	    // solid-solid block

	    for (unsigned int j=0; j != n_solid_dofs; j++)
	      {

		//lambda terms don't vary with solid. So we don't need this.

		/*
		// Computing the grad_lambda and lambda derivative terms w.r.t solid

		u_coeffs(j) += delta;

		libMesh::Gradient grad_lmx_upd, grad_lmy_upd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_upd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_upd);

		libMesh::Real lmx_upd, lmy_upd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_upd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_upd);

		u_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_lmx_umd, grad_lmy_umd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_umd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_umd);

		libMesh::Real lmx_umd, lmy_umd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_umd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_umd);

		u_coeffs(j) += delta;

		v_coeffs(j) += delta;

		libMesh::Gradient grad_lmx_vpd, grad_lmy_vpd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_vpd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_vpd);

		libMesh::Real lmx_vpd, lmy_vpd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_vpd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_vpd);

		v_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_lmx_vmd, grad_lmy_vmd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_vmd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_vmd);

		libMesh::Real lmx_vmd, lmy_vmd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_vmd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_vmd);

		v_coeffs(j) += delta;

		//H1 Norm

		// Finite differencing the grad_lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += ((((grad_lmx_upd(alpha)-grad_lmx_umd(alpha))/(2*delta))
				     -P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += ((((grad_lmy_upd(alpha)-grad_lmy_umd(alpha))/(2*delta))
				     -P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_y*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += ((((grad_lmx_vpd(alpha)-grad_lmx_vmd(alpha))/(2*delta))
				     -P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += ((((grad_lmy_vpd(alpha)-grad_lmy_vmd(alpha))/(2*delta))
				     -P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_y*solid_phi[i][sqp])*jac;
		  }

		// Finite differencing the lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + ((lmx_upd-lmx_umd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + ((lmy_upd-lmy_umd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + ((lmx_vpd-lmx_vmd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + ((lmy_vpd-lmy_vmd)/(2*delta))*solid_phi[i][sqp])*jac;
		  }

		//L2 Norm
		// Finite differencing the lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += (-P(0,alpha)*solid_dphi[i][sqp](alpha)
				    + ((lmx_upd-lmx_umd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += (-P(1,alpha)*solid_dphi[i][sqp](alpha)
				    + ((lmy_upd-lmy_umd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += (-P(0,alpha)*solid_dphi[i][sqp](alpha)
				    + ((lmx_vpd-lmx_vmd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += (-P(1,alpha)*solid_dphi[i][sqp](alpha)
				    + ((lmy_vpd-lmy_vmd)/(2*delta))*solid_phi[i][sqp])*jac;
		  }
		*/

		// Finite differencing P terms

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::Gradient grad_upd, grad_vpd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_upd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vpd);

		u_coeffs(j) -= 2*delta;
		v_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_umd, grad_vmd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_umd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vmd);

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::TensorValue<libMesh::Real> P_upd;
		this->eval_first_Piola(grad_upd,grad_v,P_upd);

		libMesh::TensorValue<libMesh::Real> P_umd;
		this->eval_first_Piola(grad_umd,grad_v,P_umd);

		libMesh::TensorValue<libMesh::Real> P_vpd;
		this->eval_first_Piola(grad_u,grad_vpd,P_vpd);

		libMesh::TensorValue<libMesh::Real> P_vmd;
		this->eval_first_Piola(grad_u,grad_vmd,P_vmd);
		/*
		//H1 Norm
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += ((grad_lambda_x(alpha)
				     -((P_upd(0,alpha)-P_umd(0,alpha))/(2*delta)))*solid_dphi[i][sqp](alpha)
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += ((grad_lambda_y(alpha)
				     -((P_upd(1,alpha)-P_umd(1,alpha))/(2*delta)))*solid_dphi[i][sqp](alpha)
				    + lambda_y*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += ((grad_lambda_x(alpha)
				     -((P_vpd(0,alpha)-P_vmd(0,alpha))/(2*delta)))*solid_dphi[i][sqp](alpha)
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += ((grad_lambda_y(alpha)
				     -((P_vpd(1,alpha)-P_vmd(1,alpha))/(2*delta)))*solid_dphi[i][sqp](alpha)
				    + lambda_y*solid_phi[i][sqp])*jac;
		  }
		*/
		//L2 Norm

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) -= ((P_upd(0,alpha)-P_umd(0,alpha))/(2*delta))*solid_dphi[i][sqp](alpha)*jac;

		    Kvs_us(i,j) -= ((P_upd(1,alpha)-P_umd(1,alpha))/(2*delta))*solid_dphi[i][sqp](alpha)*jac;

		    Kus_vs(i,j) -= ((P_vpd(0,alpha)-P_vmd(0,alpha))/(2*delta))*solid_dphi[i][sqp](alpha)*jac;

		    Kvs_vs(i,j) -= ((P_vpd(1,alpha)-P_vmd(1,alpha))/(2*delta))*solid_dphi[i][sqp](alpha)*jac;
		  }

		//Solid phi and dphi wont change. No need to differentiate.

		/*
		// Compute solid_phi and solid_dphi derivative terms

		u_coeffs(j) += delta;

		const std::vector<std::vector<libMesh::Real> > solid_phi_upd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();
		const std::vector<std::vector<libMesh::RealGradient> > solid_dphi_upd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

		u_coeffs(j) -= 2*delta;

		const std::vector<std::vector<libMesh::Real> > solid_phi_umd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();
		const std::vector<std::vector<libMesh::RealGradient> > solid_dphi_umd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

		u_coeffs(j) += delta;

		v_coeffs(j) += delta;

		const std::vector<std::vector<libMesh::Real> > solid_phi_vpd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();
		const std::vector<std::vector<libMesh::RealGradient> > solid_dphi_vpd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

		v_coeffs(j) -= 2*delta;

		const std::vector<std::vector<libMesh::Real> > solid_phi_vmd =
		  solid_context.get_element_fe(this->_disp_vars.u(),2)->get_phi();
		const std::vector<std::vector<libMesh::RealGradient> > solid_dphi_vmd
		  = solid_context.get_element_fe(this->_disp_vars.u(),2)->get_dphi();

		v_coeffs(j) += delta;

		//H1 Norm
		//Finite differencing the solid_dphi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))
				    *((solid_dphi_upd[i][sqp](alpha)-solid_dphi_umd[i][sqp](alpha))/(2*delta))
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))
				    *((solid_dphi_upd[i][sqp](alpha)-solid_dphi_umd[i][sqp](alpha))/(2*delta))
				    + lambda_y*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))
				    *((solid_dphi_vpd[i][sqp](alpha)-solid_dphi_vmd[i][sqp](alpha))/(2*delta))
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))
				    *((solid_dphi_vpd[i][sqp](alpha)-solid_dphi_vmd[i][sqp](alpha))/(2*delta))
				    + lambda_y*solid_phi[i][sqp])*jac;
		  }

		//Finite differencing the solid_phi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_x*((solid_phi_upd[i][sqp]-solid_phi_umd[i][sqp])/(2*delta)))*jac;

		    Kvs_us(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_y*((solid_phi_upd[i][sqp]-solid_phi_umd[i][sqp])/(2*delta)))*jac;

		    Kus_vs(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_x*((solid_phi_vpd[i][sqp]-solid_phi_vmd[i][sqp])/(2*delta)))*jac;

		    Kvs_vs(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				    + lambda_y*((solid_phi_vpd[i][sqp]-solid_phi_vmd[i][sqp])/(2*delta)))*jac;
		  }

		//L2 Norm
		//Finite differencing the solid_dphi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += (-P(0,alpha)*((solid_dphi_upd[i][sqp](alpha)-solid_dphi_umd[i][sqp](alpha))/(2*delta))
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_us(i,j) += (-P(1,alpha)*((solid_dphi_upd[i][sqp](alpha)-solid_dphi_umd[i][sqp](alpha))/(2*delta))
				    + lambda_y*solid_phi[i][sqp])*jac;

		    Kus_vs(i,j) += (-P(0,alpha)*((solid_dphi_vpd[i][sqp](alpha)-solid_dphi_vmd[i][sqp](alpha))/(2*delta))
				    + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_vs(i,j) += (-P(1,alpha)*((solid_dphi_vpd[i][sqp](alpha)-solid_dphi_vmd[i][sqp](alpha))/(2*delta))
				    + lambda_y*solid_phi[i][sqp])*jac;
		  }

		//Finite differencing the solid_phi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_us(i,j) += (-P(0,alpha)*solid_dphi[i][sqp](alpha)
				    + lambda_x*((solid_phi_upd[i][sqp]-solid_phi_umd[i][sqp])/(2*delta)))*jac;

		    Kvs_us(i,j) += (-P(1,alpha)*solid_dphi[i][sqp](alpha)
				    + lambda_y*((solid_phi_upd[i][sqp]-solid_phi_umd[i][sqp])/(2*delta)))*jac;

		    Kus_vs(i,j) += (-P(0,alpha)*solid_dphi[i][sqp](alpha)
				    + lambda_x*((solid_phi_vpd[i][sqp]-solid_phi_vmd[i][sqp])/(2*delta)))*jac;

		    Kvs_vs(i,j) += (-P(1,alpha)*solid_dphi[i][sqp](alpha)
				    + lambda_y*((solid_phi_vpd[i][sqp]-solid_phi_vmd[i][sqp])/(2*delta)))*jac;
		  }
		*/

	      } //solid dof loop

	    // solid-lambda block

	    for (unsigned int j=0; j != n_lambda_dofs; j++)
	      {

		// Computing the grad_lambda and lambda derivative terms w.r.t lambda

		lambda_xcoeff(j) += delta;
		/*
		libMesh::Gradient grad_lmx_xpd, grad_lmy_xpd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_xpd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_xpd);
		*/
		libMesh::Real lmx_xpd, lmy_xpd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_xpd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_xpd);

		lambda_xcoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_lmx_xmd, grad_lmy_xmd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_xmd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_xmd);
		*/
		libMesh::Real lmx_xmd, lmy_xmd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_xmd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_xmd);

		lambda_xcoeff(j) += delta;

		lambda_ycoeff(j) += delta;
		/*
		libMesh::Gradient grad_lmx_ypd, grad_lmy_ypd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_ypd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_ypd);
		*/
		libMesh::Real lmx_ypd, lmy_ypd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_ypd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_ypd);

		lambda_ycoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_lmx_ymd, grad_lmy_ymd;
		solid_context.interior_gradient(this->_lambda_var.u(), sqp, grad_lmx_ymd);
		solid_context.interior_gradient(this->_lambda_var.v(), sqp, grad_lmy_ymd);
		*/
		libMesh::Real lmx_ymd, lmy_ymd;
		solid_context.interior_value(this->_lambda_var.u(), sqp, lmx_ymd);
		solid_context.interior_value(this->_lambda_var.v(), sqp, lmy_ymd);

		lambda_ycoeff(j) += delta;
		/*
		//H1 Norm
		//Finite differencing the grad_lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_ulm(i,j) += ((((grad_lmx_xpd(alpha)-grad_lmx_xmd(alpha))/(2*delta))
				      - P(0,alpha))*solid_dphi[i][sqp](alpha)
				     + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_ulm(i,j) += ((((grad_lmy_xpd(alpha)-grad_lmy_xmd(alpha))/(2*delta))
				      - P(1,alpha))*solid_dphi[i][sqp](alpha)
				     + lambda_y*solid_phi[i][sqp])*jac;

		    Kus_vlm(i,j) += ((((grad_lmx_ypd(alpha)-grad_lmx_ymd(alpha))/(2*delta))
				      - P(0,alpha))*solid_dphi[i][sqp](alpha)
				     + lambda_x*solid_phi[i][sqp])*jac;

		    Kvs_vlm(i,j) += ((((grad_lmy_ypd(alpha)-grad_lmy_ymd(alpha))/(2*delta))
				      - P(1,alpha))*solid_dphi[i][sqp](alpha)
				     + lambda_y*solid_phi[i][sqp])*jac;
		  }

		//Finite differencing the lambda terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kus_ulm(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				     + ((lmx_xpd-lmx_xmd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_ulm(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				     + ((lmy_xpd-lmy_xmd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kus_vlm(i,j) += ((grad_lambda_x(alpha)-P(0,alpha))*solid_dphi[i][sqp](alpha)
				     + ((lmx_ypd-lmx_ymd)/(2*delta))*solid_phi[i][sqp])*jac;

		    Kvs_vlm(i,j) += ((grad_lambda_y(alpha)-P(1,alpha))*solid_dphi[i][sqp](alpha)
				     + ((lmy_ypd-lmy_ymd)/(2*delta))*solid_phi[i][sqp])*jac;
		  }
		*/
		//L2 Norm
		//Finite differencing the lambda terms

		Kus_ulm(i,j) += ((lmx_xpd-lmx_xmd)/(2*delta))*solid_phi[i][sqp]*jac;

		Kvs_ulm(i,j) += ((lmy_xpd-lmy_xmd)/(2*delta))*solid_phi[i][sqp]*jac;

		Kus_vlm(i,j) += ((lmx_ypd-lmx_ymd)/(2*delta))*solid_phi[i][sqp]*jac;

		Kvs_vlm(i,j) += ((lmy_ypd-lmy_ymd)/(2*delta))*solid_phi[i][sqp]*jac;

	      } //lambda dof loop

	  }// compute_jacobian

      } //solid dof loop

  }

  template<typename SolidMech>
  void ImmersedBoundary<SolidMech>::lambda_residual_contribution( bool compute_jacobian, MultiphysicsSystem & system,
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
								  libMesh::DenseSubMatrix<libMesh::Number> & Kvlm_vs)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
    unsigned int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

    libMesh::DenseSubVector<libMesh::Number> & u_coeffs = solid_context.get_elem_solution(this->_disp_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & v_coeffs = solid_context.get_elem_solution(this->_disp_vars.v());

    libMesh::DenseSubVector<libMesh::Number> & fluid_ucoeff = fluid_context.get_elem_solution(this->_flow_vars.u());
    libMesh::DenseSubVector<libMesh::Number> & fluid_vcoeff = fluid_context.get_elem_solution(this->_flow_vars.v());

    // Prepare lagrange multiplier info needed
    const std::vector<std::vector<libMesh::Real> > lambda_phi =
      solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();

    const std::vector<std::vector<libMesh::Real> > fluid_phi =
      fluid_context.get_element_fe(this->_flow_vars.u())->get_phi();

    //const std::vector<std::vector<libMesh::RealGradient> > & lambda_dphi =
    //  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

    // Compute the fluid velocity at the solid element quadrature points.
    libMesh::Real Vx, Vy;
    fluid_context.interior_value(this->_flow_vars.u(), 0, Vx);
    fluid_context.interior_value(this->_flow_vars.v(), 0, Vy);
    /*
    libMesh::Gradient grad_Vx, grad_Vy;
    fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx);
    fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy);

    libMesh::Gradient grad_u, grad_v;
    solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_u);
    solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_v);
    */
    libMesh::Real udot, vdot;
    solid_context.interior_rate(this->_disp_vars.u(), sqp, udot);
    solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot);
    /*
    libMesh::TensorValue<libMesh::Real> F;
    this->eval_deform_gradient(grad_u,grad_v,F);

    //Computing Fdot and gradV_times_F
    libMesh::TensorValue<libMesh::Real> gradV_times_F;
    libMesh::TensorValue<libMesh::Real> Fdot;

    for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    solid_context.interior_rate(F(0,alpha), sqp, Fdot(0,alpha));
	    solid_context.interior_rate(F(1,alpha), sqp, Fdot(1,alpha));

	    for( unsigned int beta = 0; beta < 2; beta++ )
	      {
		gradV_times_F(0,alpha) += grad_Vx(beta)*F(beta,alpha);
		gradV_times_F(1,alpha) += grad_Vy(beta)*F(beta,alpha);
	      }
   	  }
    */

    for( unsigned int i = 0; i < n_lambda_dofs; i++ )
      {
	/*
	//H1 Norm
	for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
	  {
	    Fulm(i) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
			+ lambda_phi[i][sqp]*(Vx - udot))*jac;

	    Fvlm(i) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
			+ lambda_phi[i][sqp]*(Vy - vdot))*jac;
	  }
	*/
	//L2 Norm

	Fulm(i) += lambda_phi[i][sqp]*(Vx - udot)*jac;

	Fvlm(i) += lambda_phi[i][sqp]*(Vy - vdot)*jac;


	if( compute_jacobian )
	  {
	    // lambda-fluid block
	    for( unsigned int j = 0; j < n_fluid_dofs; j++ )
	      {
		libmesh_assert_equal_to( fluid_phi[j].size(), 1 );

		// Computing V and grad_V terms fluid derivative terms

		fluid_ucoeff(j) += delta;
		/*
		libMesh::Gradient grad_Vx_upd, grad_Vy_upd;
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_upd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_upd);
		*/
		libMesh::Real Vx_upd, Vy_upd;
		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_upd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_upd);

		fluid_ucoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_Vx_umd, grad_Vy_umd;
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_umd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_umd);
		*/
		libMesh::Real Vx_umd, Vy_umd;
		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_umd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_umd);

		fluid_ucoeff(j) += delta;

		fluid_vcoeff(j) += delta;
		/*
		libMesh::Gradient grad_Vx_vpd, grad_Vy_vpd;
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_vpd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_vpd);
		*/
		libMesh::Real Vx_vpd, Vy_vpd;
		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_vpd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_vpd);

		fluid_vcoeff(j) -= 2*delta;
		/*
		libMesh::Gradient grad_Vx_vmd, grad_Vy_vmd;
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_vmd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_vmd);
		*/
		libMesh::Real Vx_vmd, Vy_vmd;
		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_vmd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_vmd);

		fluid_vcoeff(j) += delta;
		/*
		//H1 Norm
		// Finite differencing the grad_V terms w.r.t fluid
		libMesh::TensorValue<libMesh::Real> gradVtF_uf, gradVtF_vf;

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    for( unsigned int beta = 0; beta < 2; beta++ )
		      {
			gradVtF_uf(0,alpha) += ((grad_Vx_upd(beta)-grad_Vx_umd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_uf(1,alpha) += ((grad_Vy_upd(beta)-grad_Vy_umd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_vf(0,alpha) += ((grad_Vx_vpd(beta)-grad_Vx_vmd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_vf(1,alpha) += ((grad_Vy_vpd(beta)-grad_Vy_vmd(beta))/(2*delta))*F(beta,alpha);
		      }
		  }


		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_uf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_uf(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_uf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_uf(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;

		    Kulm_vf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_vf(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_vf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_vf(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;
		  }

		//Finite differencing the V terms w.r.t fluid
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_uf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(((Vx_upd-Vx_umd)/(2*delta)) - udot))*jac;

		    Kvlm_uf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(((Vy_upd-Vy_umd)/(2*delta)) - vdot))*jac;

		    Kulm_vf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(((Vx_vpd-Vx_vmd)/(2*delta)) - udot))*jac;

		    Kvlm_vf(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(((Vy_vpd-Vy_vmd)/(2*delta)) - vdot))*jac;
		  }
		*/
		//L2 Norm
		//Finite differencing the V terms w.r.t fluid

		Kulm_uf(i,j) += lambda_phi[i][sqp]*((Vx_upd-Vx_umd)/(2*delta))*jac;

		Kvlm_uf(i,j) += lambda_phi[i][sqp]*((Vy_upd-Vy_umd)/(2*delta))*jac;

		Kulm_vf(i,j) += lambda_phi[i][sqp]*((Vx_vpd-Vx_vmd)/(2*delta))*jac;

		Kvlm_vf(i,j) += lambda_phi[i][sqp]*((Vy_vpd-Vy_vmd)/(2*delta))*jac;


	      } //fluid dof loop

	    // lambda-solid block
	    for( unsigned int j = 0; j < n_solid_dofs; j++ )
	      {
		/*
		// Computing Fdot and F derivative terms

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::Gradient grad_upd, grad_vpd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_upd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vpd);

		u_coeffs(j) -= 2*delta;
		v_coeffs(j) -= 2*delta;

		libMesh::Gradient grad_umd, grad_vmd;
		solid_context.interior_gradient(this->_disp_vars.u(), sqp, grad_umd);
		solid_context.interior_gradient(this->_disp_vars.v(), sqp, grad_vmd);

		u_coeffs(j) += delta;
		v_coeffs(j) += delta;

		libMesh::TensorValue<libMesh::Real> F_upd;
		this->eval_deform_gradient(grad_upd,grad_v,F_upd);

		libMesh::TensorValue<libMesh::Real> Fdot_upd;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    solid_context.interior_rate(F_upd(0,alpha), sqp, Fdot_upd(0,alpha));
		    solid_context.interior_rate(F_upd(1,alpha), sqp, Fdot_upd(1,alpha));
		  }

		libMesh::TensorValue<libMesh::Real> F_umd;
		this->eval_deform_gradient(grad_umd,grad_v,F_umd);

		libMesh::TensorValue<libMesh::Real> Fdot_umd;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    solid_context.interior_rate(F_umd(0,alpha), sqp, Fdot_umd(0,alpha));
		    solid_context.interior_rate(F_umd(1,alpha), sqp, Fdot_umd(1,alpha));
		  }

		libMesh::TensorValue<libMesh::Real> F_vpd;
		this->eval_deform_gradient(grad_u,grad_vpd,F_vpd);

		libMesh::TensorValue<libMesh::Real> Fdot_vpd;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    solid_context.interior_rate(F_vpd(0,alpha), sqp, Fdot_vpd(0,alpha));
		    solid_context.interior_rate(F_vpd(1,alpha), sqp, Fdot_vpd(1,alpha));
		  }

		libMesh::TensorValue<libMesh::Real> F_vmd;
		this->eval_deform_gradient(grad_u,grad_vmd,F_vmd);

		libMesh::TensorValue<libMesh::Real> Fdot_vmd;
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    solid_context.interior_rate(F_vmd(0,alpha), sqp, Fdot_vmd(0,alpha));
		    solid_context.interior_rate(F_vmd(1,alpha), sqp, Fdot_vmd(1,alpha));
		  }

		// Finite differencing the Fdot terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha)
								 - ((Fdot_upd(0,alpha)-Fdot_umd(0,alpha))/(2*delta)))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha)
								 - ((Fdot_upd(1,alpha)-Fdot_umd(1,alpha))/(2*delta)))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha)
								 - ((Fdot_vpd(0,alpha)-Fdot_vmd(0,alpha))/(2*delta)))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha)
								 - ((Fdot_vpd(1,alpha)-Fdot_vmd(1,alpha))/(2*delta)))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;
		  }

		// Finite differencing the F terms
		libMesh::TensorValue<libMesh::Real> gVtF_us, gVtF_vs;

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    for( unsigned int beta = 0; beta < 2; beta++ )
		      {
			gVtF_us(0,alpha) += grad_Vx(beta)*((F_upd(beta,alpha)-F_umd(beta,alpha))/(2*delta));
			gVtF_us(1,alpha) += grad_Vy(beta)*((F_upd(beta,alpha)-F_umd(beta,alpha))/(2*delta));
			gVtF_vs(0,alpha) += grad_Vx(beta)*((F_vpd(beta,alpha)-F_vmd(beta,alpha))/(2*delta));
			gVtF_vs(1,alpha) += grad_Vy(beta)*((F_vpd(beta,alpha)-F_vmd(beta,alpha))/(2*delta));
		      }
		  }

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gVtF_us(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gVtF_us(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gVtF_vs(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gVtF_vs(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;
		  }
		*/
		// Computing V and grad_V terms solid derivative terms

		libMesh::Real Vx_upd, Vy_upd, Vx_umd, Vy_umd;
		//libMesh::Gradient grad_Vx_upd, grad_Vy_upd, grad_Vx_umd, grad_Vy_umd;

		u_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_upd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_upd);
		/*
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_upd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_upd);
		*/
		u_coeffs(j) -= 2*delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_umd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_umd);
		/*
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_umd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_umd);
		*/
		u_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		libMesh::Real Vx_vpd, Vy_vpd, Vx_vmd, Vy_vmd;
		//libMesh::Gradient grad_Vx_vpd, grad_Vy_vpd, grad_Vx_vmd, grad_Vy_vmd;

		v_coeffs(j) += delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_vpd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_vpd);
		/*
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_vpd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_vpd);
		*/
		v_coeffs(j) -= 2*delta;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		fluid_context.interior_value(this->_flow_vars.u(), 0, Vx_vmd);
		fluid_context.interior_value(this->_flow_vars.v(), 0, Vy_vmd);
		/*
		fluid_context.interior_gradient(this->_flow_vars.u(), 0, grad_Vx_vmd);
		fluid_context.interior_gradient(this->_flow_vars.v(), 0, grad_Vy_vmd);
		*/
		v_coeffs(j) += delta;
		/*
		//H1 Norm
		// Finite differencing the grad_V terms w.r.t solid
		libMesh::TensorValue<libMesh::Real> gradVtF_us, gradVtF_vs;

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    for( unsigned int beta = 0; beta < 2; beta++ )
		      {
			gradVtF_us(0,alpha) += ((grad_Vx_upd(beta)-grad_Vx_umd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_us(1,alpha) += ((grad_Vy_upd(beta)-grad_Vy_umd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_vs(0,alpha) += ((grad_Vx_vpd(beta)-grad_Vx_vmd(beta))/(2*delta))*F(beta,alpha);
			gradVtF_vs(1,alpha) += ((grad_Vy_vpd(beta)-grad_Vy_vmd(beta))/(2*delta))*F(beta,alpha);
		      }
		  }

		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_us(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_us(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_vs(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradVtF_vs(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;
		  }

		//Finite differencing the V terms w.r.t solid
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(((Vx_upd-Vx_umd)/(2*delta)) - udot))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(((Vy_upd-Vy_umd)/(2*delta)) - vdot))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(((Vx_vpd-Vx_vmd)/(2*delta)) - udot))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(((Vy_vpd-Vy_vmd)/(2*delta)) - vdot))*jac;
		  }
		*/
		//L2 Norm
		//Finite differencing the V terms w.r.t solid

		Kulm_us(i,j) += lambda_phi[i][sqp]*((Vx_upd-Vx_umd)/(2*delta))*jac;

		Kvlm_us(i,j) += lambda_phi[i][sqp]*((Vy_upd-Vy_umd)/(2*delta))*jac;

		Kulm_vs(i,j) += lambda_phi[i][sqp]*((Vx_vpd-Vx_vmd)/(2*delta))*jac;

		Kvlm_vs(i,j) += lambda_phi[i][sqp]*((Vy_vpd-Vy_vmd)/(2*delta))*jac;

		this->prepare_fluid_context(system,solid_context,solid_qpoints,sqp,fluid_elem_id,fluid_context);

		//Finite differencing the udot and vdot terms

		u_coeffs(j) += delta;

		libMesh::Real udot_upd, vdot_upd;
		solid_context.interior_rate(this->_disp_vars.u(), sqp, udot_upd);
		solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot_upd);

		u_coeffs(j) -= 2*delta;

		libMesh::Real udot_umd, vdot_umd;
		solid_context.interior_rate(this->_disp_vars.u(), sqp, udot_umd);
		solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot_umd);

		u_coeffs(j) += delta;

		v_coeffs(j) += delta;

		libMesh::Real udot_vpd, vdot_vpd;
		solid_context.interior_rate(this->_disp_vars.u(), sqp, udot_vpd);
		solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot_vpd);

		v_coeffs(j) -= 2*delta;

		libMesh::Real udot_vmd, vdot_vmd;
		solid_context.interior_rate(this->_disp_vars.u(), sqp, udot_vmd);
		solid_context.interior_rate(this->_disp_vars.v(), sqp, vdot_vmd);

		v_coeffs(j) += delta;
		/*
		//H1 Norm
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - ((udot_upd-udot_umd)/(2*delta))))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - ((vdot_upd-vdot_umd)/(2*delta))))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - ((udot_vpd-udot_vmd)/(2*delta))))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - ((vdot_vpd-vdot_vmd)/(2*delta))))*jac;
		  }
		*/
		//L2 Norm

		Kulm_us(i,j) -= lambda_phi[i][sqp]*((udot_upd-udot_umd)/(2*delta))*jac;

		Kvlm_us(i,j) -= lambda_phi[i][sqp]*((vdot_upd-vdot_umd)/(2*delta))*jac;

		Kulm_vs(i,j) -= lambda_phi[i][sqp]*((udot_vpd-udot_vmd)/(2*delta))*jac;

		Kvlm_vs(i,j) -= lambda_phi[i][sqp]*((vdot_vpd-vdot_vmd)/(2*delta))*jac;


		//lambda terms not dependent on solid. No need to differentiate

		/*
		// Computing lambda_phi and lambda_dphi derivative terms

		u_coeffs(j) += delta;

		const std::vector<std::vector<libMesh::Real> > lambda_phi_upd =
		  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > lambda_dphi_upd =
		//  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

		u_coeffs(j) -= 2*delta;

		const std::vector<std::vector<libMesh::Real> > lambda_phi_umd =
		  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > lambda_dphi_umd =
		//  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

		u_coeffs(j) += delta;

		v_coeffs(j) += delta;

		const std::vector<std::vector<libMesh::Real> > lambda_phi_vpd =
		  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > lambda_dphi_vpd =
		//  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

		v_coeffs(j) -= 2*delta;

		const std::vector<std::vector<libMesh::Real> > lambda_phi_vmd =
		  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_phi();
		//const std::vector<std::vector<libMesh::RealGradient> > lambda_dphi_vmd =
		//  solid_context.get_element_fe(this->_lambda_var.u(),2)->get_dphi();

		v_coeffs(j) += delta;

		//H1 Norm
		// Finite differencing the lambda_dphi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (((lambda_dphi_upd[i][sqp](alpha)-lambda_dphi_umd[i][sqp](alpha))/(2*delta))
				     *(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_us(i,j) += (((lambda_dphi_upd[i][sqp](alpha)-lambda_dphi_umd[i][sqp](alpha))/(2*delta))
				     *(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;

		    Kulm_vs(i,j) += (((lambda_dphi_vpd[i][sqp](alpha)-lambda_dphi_vmd[i][sqp](alpha))/(2*delta))
				     *(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + lambda_phi[i][sqp]*(Vx - udot))*jac;

		    Kvlm_vs(i,j) += (((lambda_dphi_vpd[i][sqp](alpha)-lambda_dphi_vmd[i][sqp](alpha))/(2*delta))
				     *(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + lambda_phi[i][sqp]*(Vy - vdot))*jac;
		  }


		//Finite differencing the lambda_phi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + ((lambda_phi_upd[i][sqp]-lambda_phi_umd[i][sqp])/(2*delta))*(Vx - udot))*jac;

		    Kvlm_us(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + ((lambda_phi_upd[i][sqp]-lambda_phi_umd[i][sqp])/(2*delta))*(Vy - vdot))*jac;

		    Kulm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(0,alpha) - Fdot(0,alpha))
				     + ((lambda_phi_vpd[i][sqp]-lambda_phi_vmd[i][sqp])/(2*delta))*(Vx - udot))*jac;

		    Kvlm_vs(i,j) += (lambda_dphi[i][sqp](alpha)*(gradV_times_F(1,alpha) - Fdot(1,alpha))
				     + ((lambda_phi_vpd[i][sqp]-lambda_phi_vmd[i][sqp])/(2*delta))*(Vy - vdot))*jac;
		  }

		//L2 Norm
		//Finite differencing the lambda_phi terms
		for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
		  {
		    Kulm_us(i,j) += ((lambda_phi_upd[i][sqp]-lambda_phi_umd[i][sqp])/(2*delta))*(Vx - udot)*jac;

		    Kvlm_us(i,j) += ((lambda_phi_upd[i][sqp]-lambda_phi_umd[i][sqp])/(2*delta))*(Vy - vdot)*jac;

		    Kulm_vs(i,j) += ((lambda_phi_vpd[i][sqp]-lambda_phi_vmd[i][sqp])/(2*delta))*(Vx - udot)*jac;

		    Kvlm_vs(i,j) += ((lambda_phi_vpd[i][sqp]-lambda_phi_vmd[i][sqp])/(2*delta))*(Vy - vdot)*jac;
		  }

		*/

	      } // solid dof loop

	  } // if compute jacobian

      } // lambda dof loop

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

    // The F^T comes from needing the derivative of the fluid
    // shape function w.r.t. solid coordinates
    //tau = P*Ftrans;
  }

  //instantiate IBM classes
  template class ImmersedBoundary<ElasticCable<HookesLaw1D> >;
  template class ImmersedBoundary<ElasticMembrane<HookesLaw> >;
  template class ImmersedBoundary<ElasticMembrane<IncompressiblePlaneStressHyperelasticity<MooneyRivlin > > >;

} // namespace GRINS
