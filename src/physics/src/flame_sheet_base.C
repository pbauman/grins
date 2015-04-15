//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/flame_sheet_base.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  FlameSheetBase<Mixture,Evaluator>::FlameSheetBase( const PhysicsName& physics_name, const GetPot& input )
    : _flow_vars(input,physics_name),
      _temp_vars(input,physics_name),
      _S_FE_family( libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>(input("Physics/"+physics_name+"/S_FE_family", "LAGRANGE")) ),
      _S_FE_order( libMesh::Utility::string_to_enum<GRINSEnums::Order>(input("Physics/"+physics_name+"/S_FE_order", "FIRST")) ),
      _gas_mixture(input)
  {}

  template<typename Mixture, typename Evaluator>
  void FlameSheetBase<Mixture,Evaluator>::init_variables( libMesh::FEMSystem* system )
  {
    _flow_vars.init(system);
    _temp_vars.init(system);

    _S_var =  system->add_variable( "S", this->_S_FE_order, _S_FE_family);
  }

  template<typename Mixture, typename Evaluator>
  void FlameSheetBase<Mixture,Evaluator>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(_flow_vars.u_var());
    system->time_evolving(_flow_vars.v_var());

    if (system->get_mesh().mesh_dimension() == 3)
      system->time_evolving(_flow_vars.w_var());

    system->time_evolving(_temp_vars.T_var());
    system->time_evolving(_S_var);
  }

  template<typename Mixture, typename Evaluator>
  void FlameSheetBase<Mixture,Evaluator>::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_flow_vars.u_var())->get_JxW();
    context.get_element_fe(_flow_vars.u_var())->get_phi();
    context.get_element_fe(_flow_vars.u_var())->get_dphi();
    context.get_element_fe(_flow_vars.u_var())->get_xyz();

    context.get_element_fe(_temp_vars.T_var())->get_JxW();
    context.get_element_fe(_temp_vars.T_var())->get_phi();
    context.get_element_fe(_temp_vars.T_var())->get_dphi();
    context.get_element_fe(_temp_vars.T_var())->get_xyz();

    context.get_element_fe(_S_var)->get_JxW();
    context.get_element_fe(_S_var)->get_phi();
    context.get_element_fe(_S_var)->get_dphi();
    context.get_element_fe(_S_var)->get_xyz();

    context.get_element_fe(_flow_vars.p_var())->get_phi();
    context.get_element_fe(_flow_vars.p_var())->get_xyz();
  }

} // end namespace GRINS
