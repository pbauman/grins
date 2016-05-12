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


#include "grins/immersed_boundary.h"


// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/physics_naming.h"
#include "grins/inc_nav_stokes_macro.h"
#include "grins/materials_parsing.h"
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"

// libMesh
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/fem_context.h"

namespace GRINS
{

  template<typename SolidMechanicsAbstract, typename Mu>
  void ImmersedBoundary<SolidMechanicsAbstract, Mu>::init_variables( libMesh::FEMSystem* system )
  {
    this->_dim = system->get_mesh().mesh_dimension();

    this->_flow_vars.init(system);
    this->_press_var.init(system);

    this->_mu.init(system);

    return;
  }

  template<typename SolidMechanicsAbstract, typename Mu>
  void ImmersedBoundary<SolidMechanicsAbstract, Mu>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    const unsigned int dim = system->get_mesh().mesh_dimension();

    // Tell the system to march velocity forward in time, but
    // leave p as a constraint only
    system->time_evolving(_flow_vars.u());
    system->time_evolving(_flow_vars.v());

    if (dim == 3)
      system->time_evolving(_flow_vars.w());

    return;
  }

  template<typename SolidMechanicsAbstract, typename Mu>
  void ImmersedBoundary<SolidMechanicsAbstract, Mu>::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_flow_vars.u())->get_JxW();
    context.get_element_fe(_flow_vars.u())->get_phi();
    context.get_element_fe(_flow_vars.u())->get_dphi();
    context.get_element_fe(_flow_vars.u())->get_xyz();

    context.get_element_fe(_press_var.p())->get_phi();
    context.get_element_fe(_press_var.p())->get_xyz();

    context.get_side_fe(_flow_vars.u())->get_JxW();
    context.get_side_fe(_flow_vars.u())->get_phi();
    context.get_side_fe(_flow_vars.u())->get_dphi();
    context.get_side_fe(_flow_vars.u())->get_xyz();

    return;
  }

  template<typename SolidMechanicsAbstract, typename Mu>
  libMesh::Real ImmersedBoundary<SolidMechanicsAbstract, Mu>::get_viscosity_value(AssemblyContext& context, unsigned int qp) const
  {
    return this->_mu(context, qp);
  }


  

} // namespace GRINS
