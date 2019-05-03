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

// libMesh
#include "libmesh/elem.h"

namespace GRINS
{
  template<unsigned int Dim, bool UseOldDisplacement>
  CartesianFictitiousDomainFluidStructureInteractionBase<Dim,UseOldDisplacement>::
  CartesianFictitiousDomainFluidStructureInteractionBase( const PhysicsName & physics_name,
                                                          const GetPot & input )
    : FictitiousDomainFluidStructureInteractionAbstract(physics_name,input)
  {
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

  // Instantiate
  template class CartesianFictitiousDomainFluidStructureInteractionBase<2,false>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<3,false>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<2,true>;
  template class CartesianFictitiousDomainFluidStructureInteractionBase<3,true>;

} // end namespace GRINS
