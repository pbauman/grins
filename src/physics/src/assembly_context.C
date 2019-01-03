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
#include "grins/assembly_context.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/unsteady_solver.h"

namespace GRINS
{
  MultiphysicsSystem & AssemblyContext::get_multiphysics_system()
  {
    libMesh::System & base_system = const_cast<libMesh::System &>(this->get_system());

    MultiphysicsSystem & multiphysics_system =
      libMesh::cast_ref<MultiphysicsSystem &>( base_system );

    return multiphysics_system;
  }

  const MultiphysicsSystem & AssemblyContext::get_multiphysics_system() const
  {
    const MultiphysicsSystem & multiphysics_system =
      libMesh::cast_ref<const MultiphysicsSystem &>( this->get_system() );

    return multiphysics_system;
  }

  void AssemblyContext::get_old_elem_solution( const MultiphysicsSystem & system,
                                               libMesh::DenseVector<libMesh::Number> & old_elem_solution ) const
  {
    unsigned int n_dofs = this->get_elem_solution().size();
    libmesh_assert_equal_to( old_elem_solution.size(),n_dofs);

    // Extract old nonlinear solution from TimeSolver
    // Error out if this is not an UnsteadySolver
    const libMesh::TimeSolver & time_solver = system.get_time_solver();

    const libMesh::UnsteadySolver * unsteady_solver =
      dynamic_cast<const libMesh::UnsteadySolver*>(&time_solver);

    if( !unsteady_solver )
      libmesh_error_msg("ERROR: Can only call get_old_elem_solution when using an UnsteadySolver!");

    for (unsigned int i=0; i != n_dofs; ++i)
      old_elem_solution(i) =
        unsteady_solver->old_nonlinear_solution(this->get_dof_indices()[i]);
  }

} // end namespace GRINS
