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
#include "grins/cartesian_fictitious_domain_fluid_structure_interaction_base.h"

// GRINS
#include "grins/cartesian_hyperelasticity.h"
#include "grins/mooney_rivlin.h"
#include "grins/incompressible_hyperelasticity_weak_form.h"

// libMesh
#include "libmesh/elem.h"

namespace GRINS
{
  template<unsigned int Dim, bool UseOldDisplacement>
  CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::
  CartesianFictitiousDomainFluidStructureInteractionBase( const PhysicsName & physics_name,
                                                          const GetPot & input )
    : FictitiousDomainFluidStructureInteractionAbstract(physics_name,input),
      _strain_energy(nullptr)
  {
    std::string mat_string("Physics/"+physics_name+"/Solid/material");
    const std::string material = input(mat_string,"DIE!");

    _strain_energy.reset(new MooneyRivlin(input,material));

    this->check_variable_dim_consistency();
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::init_context( AssemblyContext & context )
  {
    context.get_element_fe(_disp_vars.u(),Dim)->get_JxW();
    context.get_element_fe(_disp_vars.u(),Dim)->get_phi();
    context.get_element_fe(_disp_vars.u(),Dim)->get_dphi();

    context.get_element_fe(_solid_press_var.p(),Dim)->get_JxW();
    context.get_element_fe(_solid_press_var.p(),Dim)->get_phi();

    context.get_element_fe( _lambda_var.u(),Dim)->get_dphi();
    context.get_element_fe( _lambda_var.u(),Dim)->get_phi();
    context.get_element_fe( _lambda_var.u(),Dim)->get_JxW();
  }


  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext& context )
  {
    const libMesh::Elem & current_elem = context.get_elem();

    // Only compute this if we are on a solid element
    if( this->is_solid_elem( current_elem.subdomain_id() ) &&
        this->_fluid_solid_overlap->has_overlapping_fluid_elem(current_elem.id()) )
      {
        // For clarity
        AssemblyContext & solid_context = context;
        AssemblyContext & fluid_context = *(this->_fluid_context);

        MultiphysicsSystem & system = solid_context.get_multiphysics_system();

        // Matrix for coupled Jacobian terms that are not available in the solid_context
        libMesh::DenseMatrix<libMesh::Number> Kf_s;
        libMesh::DenseSubMatrix<libMesh::Number> Kuf_us(Kf_s), Kuf_vs(Kf_s);
        libMesh::DenseSubMatrix<libMesh::Number> Kvf_us(Kf_s), Kvf_vs(Kf_s);

        std::unique_ptr<libMesh::DenseSubMatrix<libMesh::Number>> Kuf_ws, Kvf_ws, Kwf_us, Kwf_vs, Kwf_ws;
        if(Dim==3)
          {
            Kuf_ws.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_s) );
            Kvf_ws.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_s) );
            Kwf_us.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_s) );
            Kvf_ws.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_s) );
            Kwf_ws.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_s) );
          }

        libMesh::DenseMatrix<libMesh::Number> Klm_f;
        libMesh::DenseSubMatrix<libMesh::Number> Kulm_uf(Klm_f), Kulm_vf(Klm_f);
        libMesh::DenseSubMatrix<libMesh::Number> Kvlm_uf(Klm_f), Kvlm_vf(Klm_f);

        std::unique_ptr<libMesh::DenseSubMatrix<libMesh::Number>> Kulm_wf, Kvlm_wf, Kwlm_uf, Kwlm_vf, Kwlm_wf;
        if(Dim==3)
          {
            Kulm_wf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Klm_f) );
            Kvlm_wf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Klm_f) );
            Kwlm_uf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Klm_f) );
            Kwlm_vf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Klm_f) );
            Kwlm_wf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Klm_f) );
          }

        libMesh::DenseMatrix<libMesh::Number> Kf_lm;
        libMesh::DenseSubMatrix<libMesh::Number> Kuf_ulm(Kf_lm), Kuf_vlm(Kf_lm);
        libMesh::DenseSubMatrix<libMesh::Number> Kvf_ulm(Kf_lm), Kvf_vlm(Kf_lm);

        std::unique_ptr<libMesh::DenseSubMatrix<libMesh::Number>> Kuf_wlm, Kvf_wlm, Kwf_ulm, Kwf_vlm, Kwf_wlm;
        if(Dim==3)
          {
            Kuf_wlm.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_lm) );
            Kvf_wlm.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_lm) );
            Kwf_ulm.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_lm) );
            Kwf_vlm.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_lm) );
            Kwf_wlm.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Kf_lm) );
          }


        libMesh::DenseMatrix<libMesh::Number> Ks_pf;
        libMesh::DenseSubMatrix<libMesh::Number> Kus_pf(Ks_pf), Kvs_pf(Ks_pf);

        std::unique_ptr<libMesh::DenseSubMatrix<libMesh::Number>> Kws_pf;
        if(Dim==3)
          Kws_pf.reset( new libMesh::DenseSubMatrix<libMesh::Number>(Ks_pf) );

        // We need to grab the fluid elements that are overlapping with this solid elem.
        // Then, for that fluid element, extract the indices of the *solid* quadrature points
        // that are in that fluid element
        const std::set<libMesh::dof_id_type> & fluid_elem_ids =
          _fluid_solid_overlap->get_overlapping_fluid_elems(current_elem.id());

        int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
        int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();

        IncompressibleHyperelasticityWeakForm<MooneyRivlin> weak_form;

        // Loop over those fluid elements
        for( const auto & fluid_id : fluid_elem_ids )
          {
            const libMesh::Elem * fluid_elem = system.get_mesh().elem_ptr(fluid_id);

            // Quadrature points that overlap with the current fluid element
            const std::vector<unsigned int> & qp_indices =
              _fluid_solid_overlap->get_solid_qps(current_elem.id(), fluid_elem->id() );

            fluid_context.pre_fe_reinit(system,fluid_elem);

            int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();
            int n_fluid_press_dofs = fluid_context.get_dof_indices(this->_fluid_press_var.p()).size();

            // Prepare the space for Jacobians first if we're computing it
            if ( compute_jacobian )
              {
                if(Dim==2)
                  this->prepare_jacobians(n_fluid_dofs, n_solid_dofs, n_lambda_dofs, n_fluid_press_dofs,
                                          Kf_s,
                                          Kuf_us,Kuf_vs,
                                          Kvf_us,Kvf_vs,
                                          Klm_f,
                                          Kulm_uf,Kulm_vf,
                                          Kvlm_uf,Kvlm_vf,
                                          Kf_lm,
                                          Kuf_ulm,Kuf_vlm,
                                          Kvf_ulm,Kvf_vlm,
                                          Ks_pf,
                                          Kus_pf,Kvs_pf);
                else if(Dim==3)
                  this->prepare_jacobians(n_fluid_dofs, n_solid_dofs, n_lambda_dofs, n_fluid_press_dofs,
                                          Kf_s,
                                          Kuf_us,Kuf_vs,(*Kuf_ws),
                                          Kvf_us,Kvf_vs,(*Kvf_ws),
                                          (*Kwf_us),(*Kwf_vs),(*Kwf_ws),
                                          Klm_f,
                                          Kulm_uf,Kulm_vf,(*Kulm_wf),
                                          Kvlm_uf,Kvlm_vf,(*Kvlm_wf),
                                          (*Kwlm_uf),(*Kwlm_vf),(*Kwlm_wf),
                                          Kf_lm,
                                          Kuf_ulm,Kuf_vlm,(*Kuf_wlm),
                                          Kvf_ulm,Kvf_vlm,(*Kvf_wlm),
                                          (*Kwf_ulm),(*Kwf_vlm),(*Kwf_wlm),
                                          Ks_pf,
                                          Kus_pf,Kvs_pf,(*Kws_pf) );
              }


            this->compute_ibm_residuals(system,solid_context,fluid_context,qp_indices);

            this->assemble_coupled_residuals(system,fluid_context);

            if ( compute_jacobian )
              {
                 this->compute_numerical_jacobians( system,solid_context,fluid_context,
                                                    qp_indices,
                                                    Kuf_ulm,Kuf_vlm,Kvf_ulm,Kvf_vlm,
                                                    Kulm_uf,Kulm_vf,Kvlm_uf,Kvlm_vf,
                                                    Kuf_us,Kuf_vs,Kvf_us,Kvf_vs,
                                                    Kus_pf, Kvs_pf);

                 // This will assemble the coupled terms that went into the fluid residual
                 // in the fluid context and the coupled Jacobians into the global data system
                 // residual and Jacobian
                 this->assemble_coupled_jacobians<Dim>( system, solid_context, fluid_context,
                                                        n_fluid_dofs, n_solid_dofs,
                                                        n_lambda_dofs, n_fluid_press_dofs,
                                                        Kf_s, Klm_f, Kf_lm, Ks_pf );
              }

          } // loop over fluid elems

      } // if(solid_elem)
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::
  check_variable_dim_consistency() const
  {
    if( _flow_vars.dim() != Dim )
      libmesh_error_msg(this->dim_error_msg(_flow_vars.dim()));

    if( _disp_vars.dim() != Dim )
      libmesh_error_msg(this->dim_error_msg(_disp_vars.dim()));

    if( _lambda_var.dim() != Dim )
      libmesh_error_msg(this->dim_error_msg(_lambda_var.dim()));
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  std::string CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::
  dim_error_msg(unsigned int var_dim) const
  {
    std::stringstream vs, ds;
    vs << var_dim;
    ds << Dim;

    return std::string("ERROR: Expected variable to have dimension "+ds.str()+"\n, but found "+vs.str()+"!\n");
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  libMesh::Tensor CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::form_def_gradient
  ( const libMesh::Gradient & grad_u, const libMesh::Gradient & grad_v ) const
  {
    return libMesh::Tensor( 1.0+grad_u(0), grad_u(1), 0.0,
                            grad_v(0), 1.0+grad_v(1), 0.0,
                            0.0,           0.0,       1.0 );
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  libMesh::Tensor CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::form_def_gradient
  ( const libMesh::Gradient & grad_u, const libMesh::Gradient & grad_v, const libMesh::Gradient & grad_w ) const
  {
    return libMesh::Tensor( 1.0+grad_u(0), grad_u(1), grad_u(2),
                            grad_v(0), 1.0+grad_v(1), grad_v(2),
                            grad_w(0), grad_w(1), 1.0+grad_w(2) );
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  libMesh::Tensor
  CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::form_fluid_def_gradient
  (AssemblyContext & solid_context, const libMesh::Tensor & F, const unsigned int qp) const
  {
    if(UseOldDisplacement)
      {
        const MultiphysicsSystem & system = solid_context.get_multiphysics_system();

        libMesh::DenseVector<libMesh::Number> old_elem_solution(solid_context.get_elem_solution().size());
        libMesh::DenseVector<libMesh::Number> elem_solution_copy(solid_context.get_elem_solution().size());

        // This populates the old_elem_solution vector for the solid element
        // using the solution values from the prevous timestep
        solid_context.get_old_elem_solution(system,old_elem_solution);

        // Put in the old_elem_solution so we use the previous timestep values
        elem_solution_copy = solid_context.get_elem_solution();
        solid_context.get_elem_solution() = old_elem_solution;

        libMesh::Gradient grad_u, grad_v, grad_w;
        solid_context.interior_gradient(this->_disp_vars.u(), qp, grad_u);
        solid_context.interior_gradient(this->_disp_vars.v(), qp, grad_v);
        if(Dim==3)
          solid_context.interior_gradient(this->_disp_vars.w(), qp, grad_w);

        // Copy back
        solid_context.get_elem_solution() = elem_solution_copy;

        libMesh::Tensor Fold;
        if(Dim==2)
          Fold = this->form_def_gradient(grad_u,grad_v);
        else if(Dim==3)
          Fold = this->form_def_gradient(grad_u,grad_v,grad_w);

        return Fold;
      }
    else
      return F;
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_residuals
  (AssemblyContext & solid_context,const AssemblyContext & fluid_context,unsigned int qp,
   IncompressibleHyperelasticityWeakForm<MooneyRivlin> & weak_form,
   libMesh::DenseSubVector<libMesh::Number> & Fuf,
   libMesh::DenseSubVector<libMesh::Number> & Fvf,
   libMesh::DenseSubVector<libMesh::Number> & Fus,
   libMesh::DenseSubVector<libMesh::Number> & Fvs,
   libMesh::DenseSubVector<libMesh::Number> & Fulm,
   libMesh::DenseSubVector<libMesh::Number> & Fvlm,
   libMesh::DenseSubVector<libMesh::Number> & Fps)
  {


    // Variable indices
    unsigned int us_var = this->_disp_vars.u();
    unsigned int vs_var = this->_disp_vars.v();
    unsigned int ws_var = libMesh::invalid_uint;
    if(Dim==3)
      ws_var = this->_disp_vars.w();

    unsigned int ps_var = this->_solid_press_var.p();

    unsigned int uf_var = this->_flow_vars.u();
    unsigned int vf_var = this->_flow_vars.v();
    unsigned int wf_var = libMesh::invalid_uint;
    if(Dim==3)
      wf_var = this->_flow_vars.w();

    unsigned int pf_var = this->_fluid_press_var.p();

    unsigned int lx_var = this->_lambda_var.u();
    unsigned int ly_var = this->_lambda_var.v();
    unsigned int lz_var = libMesh::invalid_uint;
    if(Dim==3)
      lz_var = this->_lambda_var.w();

    int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    int n_lambda_dofs = solid_context.get_dof_indices(this->_lambda_var.u()).size();
    int n_solid_press_dofs = solid_context.get_dof_indices(this->_solid_press_var.p()).size();
    int n_fluid_dofs = fluid_context.get_dof_indices(this->_flow_vars.u()).size();

    libMesh::Real delta_rho = this->_rho_solid - this->_rho_fluid;

    const std::vector<std::vector<libMesh::Real> > & solid_phi =
          solid_context.get_element_fe(us_var,Dim)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & solid_dphi =
      solid_context.get_element_fe(us_var,Dim)->get_dphi();

    const std::vector<std::vector<libMesh::Real> > & lambda_phi =
          solid_context.get_element_fe(lx_var,Dim)->get_phi();

    const std::vector<std::vector<libMesh::Real> > & solid_press_phi =
          solid_context.get_element_fe(ps_var,Dim)->get_phi();

    const std::vector<libMesh::Real> & solid_JxW =
      solid_context.get_element_fe(us_var,Dim)->get_JxW();

    const std::vector<std::vector<libMesh::Real> > & fluid_phi =
              fluid_context.get_element_fe(uf_var,Dim)->get_phi();


    // Compute quantities from the solid context
    libMesh::Real udot, vdot, wdot;
    solid_context.interior_rate(us_var, qp, udot);
    solid_context.interior_rate(vs_var, qp, vdot);
    if(Dim==3)
      solid_context.interior_rate(ws_var, qp, wdot);

    libMesh::Gradient Uddot;
    this->compute_displacement_accel<Dim>(solid_context,qp,Uddot);



    libMesh::Gradient grad_u, grad_v, grad_w;
    solid_context.interior_gradient(us_var, qp, grad_u);
    solid_context.interior_gradient(vs_var, qp, grad_v);
    if(Dim==3)
      solid_context.interior_gradient(ws_var, qp, grad_w);

    libMesh::Gradient grad_udot, grad_vdot, grad_wdot;
    solid_context.interior_rate_gradient(us_var, qp, grad_udot);
    solid_context.interior_rate_gradient(vs_var, qp, grad_vdot);
    if(Dim==3)
      solid_context.interior_rate_gradient(ws_var, qp, grad_wdot);

    libMesh::TensorValue<libMesh::Real> Fdot( grad_udot(0), grad_udot(1), 0,
                                              grad_vdot(0), grad_vdot(1), 0,
                                              0, 0, 0 );
    if(Dim==3)
      {
        Fdot(0,2) = grad_udot(2);
        Fdot(1,2) = grad_vdot(2);
        Fdot(2,0) = grad_wdot(0);
        Fdot(2,1) = grad_wdot(1);
        Fdot(2,2) = grad_wdot(2);
      }



    libMesh::Real solid_press;
    solid_context.interior_value(ps_var, qp, solid_press);



    libMesh::Real lambda_x, lambda_y, lambda_z;
    solid_context.interior_value(lx_var, qp, lambda_x);
    solid_context.interior_value(ly_var, qp, lambda_y);
    if(Dim==3)
      solid_context.interior_value(lz_var, qp, lambda_z);




    libMesh::Gradient grad_lambda_x, grad_lambda_y, grad_lambda_z;
    solid_context.interior_gradient(lx_var, qp, grad_lambda_x);
    solid_context.interior_gradient(ly_var, qp, grad_lambda_y);
    if(Dim==3)
      solid_context.interior_gradient(lz_var, qp, grad_lambda_z);




    // Compute quantities from the fluid context
    // The qp is always "0" since there should be only one point
    // in this prepared fluid context
    libMesh::Real Vx, Vy, Vz;
    fluid_context.interior_value(uf_var, 0, Vx);
    fluid_context.interior_value(vf_var, 0, Vy);
    if(Dim==3)
      fluid_context.interior_value(wf_var, 0, Vz);



    libMesh::Gradient grad_Vx, grad_Vy, grad_Vz;
    fluid_context.interior_gradient(uf_var, 0, grad_Vx);
    fluid_context.interior_gradient(vf_var, 0, grad_Vy);
    if(Dim==3)
      fluid_context.interior_gradient(wf_var, 0, grad_Vz);

    libMesh::Number pfluid;
    fluid_context.interior_value(pf_var,0,pfluid);

    // Setup the deformation gradients
    // F is the usual deformation gradient and F_fluid
    // is the one we'll used in the derivatives of
    // the fluid-shape-function-evaluated-at-deformed-solid-qps
    // quantities. For semi-implicit algorithms, F_fluid will
    // be the deformation gradient at U_n while for fully implicit
    // algorithms F_fluid = F
    libMesh::Tensor F, F_fluid;
    if(Dim==2)
      F = this->form_def_gradient(grad_u,grad_v);
    else if(Dim==3)
      F = this->form_def_gradient(grad_u,grad_v,grad_w);

    F_fluid = this->form_fluid_def_gradient(solid_context,F,qp);




    // Builds all our hyperelastic quantities
    CartesianHyperlasticity<MooneyRivlin>
      stress_law(F, (*(this->_strain_energy)));

    //libMesh::Tensor P(stress_law.pk1_stress());



    libMesh::Number J = stress_law.get_J();

    const libMesh::Tensor & Cinv = stress_law.get_C_inverse();

    libMesh::Tensor FCinv(F*Cinv);


    libMesh::TensorValue<libMesh::Real> C(F.transpose()*F);
    libMesh::Number I1 = C.tr();
    libMesh::Real J23 = std::pow(J,-2/3);
    libMesh::Real mus = 5;
    libMesh::Tensor  P( mus*J23*(F-(1/3)*I1*FCinv) );


    libMesh::Tensor grad_lam( grad_lambda_x(0), grad_lambda_x(1), 0,
                              grad_lambda_y(0), grad_lambda_y(1), 0,
                              0, 0, 0);


    if(Dim==3)
      {
        grad_lam(0,2) = grad_lambda_x(2);
        grad_lam(1,2) = grad_lambda_y(2);
        grad_lam(2,0) = grad_lambda_z(0);
        grad_lam(2,1) = grad_lambda_z(1);
        grad_lam(2,2) = grad_lambda_z(2);
      }

    libMesh::Tensor grad_lam_timesFT( grad_lam*(F_fluid.transpose()) );

    libMesh::Tensor gradV( grad_Vx(0),  grad_Vx(1), 0,
                           grad_Vy(0), grad_Vy(1), 0,
                           0, 0, 0);
    if(Dim==3)
      {
        gradV(0,2) = grad_Vx(2);
        gradV(1,2) = grad_Vy(2);
        gradV(0,2) = grad_Vz(0);
        gradV(1,2) = grad_Vz(1);
        gradV(2,2) = grad_Vz(2);
      }

    libMesh::TensorValue<libMesh::Real> gradV_times_F( gradV*F_fluid );

    libMesh::Real jac = solid_JxW[qp];

    //============================================
    // Fluid residual terms
    //============================================
    for (int i=0; i != n_fluid_dofs; i++)
      {
        libmesh_assert_equal_to( fluid_phi[i].size(), 1 );

        libMesh::Real phiJ = fluid_phi[i][0]*jac;

        if(Dim==2)
          {
            // L2 Norm
            Fuf(i) -= lambda_x*phiJ;
            Fvf(i) -= lambda_y*phiJ;

            // H1 Term
            /*
            libMesh::Gradient fluid_term = grad_lam_timesFT*fluid_dphi[i][0];
            Fuf(i) -= fluid_term(0)*jac;
            Fvf(i) -= fluid_term(1)*jac;
            */
          }

        if(Dim==3)
          {
            // L2 Norm
            Fuf(i) -= lambda_x*phiJ;
            Fvf(i) -= lambda_y*phiJ;
            //(*Fwf)(i) -= lambda_z*phiJ;

            // H1 Term
            /*
              libMesh::Gradient fluid_term = grad_lam_timesFT*fluid_dphi[i][0];
              Fuf(i) -= fluid_term(0)*jac;
              Fvf(i) -= fluid_term(1)*jac;
              (*Fwf)(i) -= fluid_term(2)*jac;
            */
          }
      }

    //============================================
    // Solid residual terms
    //============================================
    for (int i=0; i != n_solid_dofs; i++)
      {
        libMesh::Real phiJ = solid_phi[i][qp]*jac;
        libMesh::RealGradient dphiJ(solid_dphi[i][qp]*jac);

        libMesh::Number deltap = solid_press - pfluid;
        if(Dim==2)
          {

            //weak_form.evaluate_internal_stress_residual
            //  (P,dphiJ,Fus(i),Fvs(i));
            for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
              {
                Fus(i) += P(0,alpha)*solid_dphi[i][qp](alpha)*jac;
                Fvs(i) += P(1,alpha)*solid_dphi[i][qp](alpha)*jac;
              }

            weak_form.evaluate_pressure_stress_residual
              (J,deltap,FCinv,dphiJ,Fus(i),Fvs(i));

            // Acceleration term
            Fus(i) += delta_rho*Uddot(0)*phiJ;
            Fvs(i) += delta_rho*Uddot(1)*phiJ;

            //L2 Norm
            Fus(i) -= lambda_x*phiJ;
            Fvs(i) -= lambda_y*phiJ;

            //H1 Norm
            /*
            for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
              {
                Fus(i) += grad_lambda_x(alpha)*dphiJ(alpha);
                Fvs(i) += grad_lambda_y(alpha)*dphiJ(alpha);
              }
            */

          }
        else if(Dim==3)
          {
            /*
            weak_form.evaluate_internal_stress_residual
              (P,dphiJ,Fus(i),Fvs(i),(*Fws)(i));
            weak_form.evaluate_pressure_stress_residual
              (J,solid_press,FCinv,dphiJ,Fus(i),Fvs(i),(*Fws)(i));
            */
            // Acceleration term
            Fus(i) += delta_rho*Uddot(0)*phiJ;
            Fvs(i) += delta_rho*Uddot(1)*phiJ;
            //(*Fws)(i) += delta_rho*Uddot(2)*phiJ;

            //L2 Norm
            Fus(i) -= lambda_x*phiJ;
            Fvs(i) -= lambda_y*phiJ;
            //(*Fws)(i) -= lambda_z*phiJ;

            //H1 Norm
            /*
              for( unsigned int alpha = 0; alpha < this->_disp_vars.dim(); alpha++ )
              {
              Fus(i) -= grad_lambda_x(alpha)*dphiJ(alpha);
              Fvs(i) -= grad_lambda_y(alpha)*dphiJ(alpha);
              (*Fws)(i) -= grad_lambda_z(alpha)*dphiJ(alpha);
              }
            */
          }

      } // end displacement dof loop

    //============================================
    // Lambda residual terms
    //============================================
    for( int i = 0; i < n_lambda_dofs; i++ )
      {
        libMesh::Real phiJ = lambda_phi[i][qp]*jac;
        if(Dim==2)
          {
            //L2 Norm
            Fulm(i) += (Vx - udot)*phiJ;
            Fvlm(i) += (Vy - vdot)*phiJ;

            // H1 term
            /*
            libMesh::Gradient fluid_term = (gradV_times_F-Fdot)*lambda_dphi[i][qp];
            Fulm(i) += fluid_term(0)*jac;
            Fvlm(i) += fluid_term(1)*jac;
            */
          }
        if(Dim==3)
          {
            //L2 Norm
            Fulm(i) += (Vx - udot)*phiJ;
            Fvlm(i) += (Vy - vdot)*phiJ;
            //(*Fwlm)(i) += (Vz - wdot)*phiJ;

            // H1 term
            /*
            libMesh::Gradient fluid_term = (gradV_times_F-Fdot)*lambda_dphi[i][qp];
            Fulm(i) += fluid_term(0)*jac;
            Fvlm(i) += fluid_term(1)*jac;
            (*Fwlm)(i) += fluid_term(2)*jac;
            */

          }

      } // End lambda dof loop

    //============================================
    // Solid pressure residual terms
    //============================================
    for( int i = 0; i < n_solid_press_dofs; i++ )
      {
        libMesh::Real phiJ = solid_press_phi[i][qp]*jac;
        weak_form.evaluate_pressure_constraint_residual(J,phiJ,Fps(i));

      } // end solid pressure dof loop

  }


  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_ibm_residuals
  (MultiphysicsSystem & system,
   AssemblyContext & solid_context,
   AssemblyContext & fluid_context,
   const std::vector<unsigned int> & qp_indices)
  {
    // Variable indices
    unsigned int us_var = this->_disp_vars.u();
    unsigned int vs_var = this->_disp_vars.v();
    unsigned int ws_var = libMesh::invalid_uint;
    if(Dim==3)
      ws_var = this->_disp_vars.w();

    unsigned int ps_var = this->_solid_press_var.p();

    unsigned int uf_var = this->_flow_vars.u();
    unsigned int vf_var = this->_flow_vars.v();
    unsigned int wf_var = libMesh::invalid_uint;
    if(Dim==3)
      wf_var = this->_flow_vars.w();

    unsigned int pf_var = this->_fluid_press_var.p();

    unsigned int lx_var = this->_lambda_var.u();
    unsigned int ly_var = this->_lambda_var.v();
    unsigned int lz_var = libMesh::invalid_uint;
    if(Dim==3)
      lz_var = this->_lambda_var.w();

    libMesh::DenseSubVector<libMesh::Number> & Fus = solid_context.get_elem_residual(us_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvs = solid_context.get_elem_residual(vs_var);
    libMesh::DenseSubVector<libMesh::Number> * Fws = nullptr;
    if(Dim==3)
      Fws = &solid_context.get_elem_residual(ws_var);

    libMesh::DenseSubVector<libMesh::Number> & Fulm = solid_context.get_elem_residual(lx_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvlm = solid_context.get_elem_residual(ly_var);
    libMesh::DenseSubVector<libMesh::Number> * Fwlm = nullptr;
    if(Dim==3)
      Fwlm = &solid_context.get_elem_residual(lz_var);

    libMesh::DenseSubVector<libMesh::Number> & Fps = solid_context.get_elem_residual(ps_var);

    // Fluid residual data structures
    libMesh::DenseSubVector<libMesh::Number> & Fuf = fluid_context.get_elem_residual(uf_var);
    libMesh::DenseSubVector<libMesh::Number> & Fvf = fluid_context.get_elem_residual(vf_var);
    libMesh::DenseSubVector<libMesh::Number> * Fwf = nullptr;
    if(Dim==3)
      Fwf = &fluid_context.get_elem_residual(wf_var);


    const std::vector<libMesh::Point> & solid_qpoints =
      solid_context.get_element_fe(us_var,Dim)->get_xyz();

    // Extract fluid shape functions
    fluid_context.get_element_fe(uf_var,Dim)->get_phi();
    fluid_context.get_element_fe(uf_var,Dim)->get_dphi();
    fluid_context.get_element_fe(pf_var,Dim)->get_phi();
    fluid_context.get_element_fe(pf_var,Dim)->get_dphi();

    const libMesh::Elem & fluid_elem = fluid_context.get_elem();

    IncompressibleHyperelasticityWeakForm<MooneyRivlin> weak_form;

    /*
    std::cout << "Solid Element:" << std::endl;
    solid_context.get_elem().print_info();

    std::cout << "Fluid Element:" << std::endl;
    fluid_context.get_elem().print_info();
    */

    // Loop over this set of quadrature points
    for( const auto & qp : qp_indices )
      {
        const libMesh::Point & x_qp = solid_qpoints[qp];

        // Compute displacement. This will be U_{n+1} or U_n depending
        // on the compile-time value of UseOldDisplacement
        libMesh::Gradient U;
        this->compute_displaced_point<Dim,UseOldDisplacement>(system,solid_context,qp,U);

        // Reinit fluid context for current solid quadrature point
        this->reinit_fluid_context<Dim>(x_qp, U, &fluid_elem, fluid_context);

        this->compute_residuals(solid_context,fluid_context,qp,weak_form,
                                Fuf,Fvf,Fus,Fvs,Fulm,Fvlm,Fps);

      } // end qp loop
  }

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_numerical_jacobians
  (MultiphysicsSystem & system,
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
   libMesh::DenseSubMatrix<libMesh::Number> & Kvs_pf)
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

    // Compute solid pressure derivs
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

    // Compute fluid pressure derivs
    {
      libMesh::DenseSubVector<libMesh::Number> & p_coeffs =
        fluid_context.get_elem_solution(this->_fluid_press_var.p());

      this->compute_fluid_press_derivs(system,quad_points,solid_context,fluid_context,delta,
                                       backwards_solid_residual,backwards_fluid_residual,
                                       p_coeffs,Kus_pf,Kvs_pf);

    }

    // Restore original residuals
    solid_context.get_elem_residual() = original_solid_residual;
    fluid_context.get_elem_residual() = original_fluid_residual;
  }


  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::finite_difference_residuals
  (MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   AssemblyContext & fluid_context,
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

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_lambda_derivs
  (MultiphysicsSystem & system,
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

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_fluid_derivs
  (MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   AssemblyContext & fluid_context,
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

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_solid_derivs
  (MultiphysicsSystem & system,
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

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_press_derivs
  (MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   AssemblyContext & fluid_context,
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

  template<unsigned int Dim, bool UseOldDisplacement>
  void CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::compute_fluid_press_derivs
  (MultiphysicsSystem & system,
   const std::vector<unsigned int> & quad_points,
   AssemblyContext & solid_context,
   AssemblyContext & fluid_context,
   const libMesh::Real delta,
   libMesh::DenseVector<libMesh::Number> & backwards_solid_residual,
   libMesh::DenseVector<libMesh::Number> & backwards_fluid_residual,
   libMesh::DenseSubVector<libMesh::Number> & press_coeff,
   libMesh::DenseSubMatrix<libMesh::Number> & Kus,
   libMesh::DenseSubMatrix<libMesh::Number> & Kvs)
  {
    unsigned int n_solid_dofs = solid_context.get_dof_indices(this->_disp_vars.u()).size();
    unsigned int n_press_dofs = fluid_context.get_dof_indices(this->_fluid_press_var.p()).size();

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


  // Instantiate
  template class CartesianFictitiousDomainFluidStructureInteractionBase<2,false>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<3,false>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<2,true>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<3,true>;

} // end namespace GRINS
