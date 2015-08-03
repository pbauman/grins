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
#include "grins/newmark_solver.h"

// GRINS
#include "grins/grins_enums.h"
#include "grins/solver_context.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/newmark_solver.h"

// C++
#include <ctime>

namespace GRINS
{
  NewmarkSolver::NewmarkSolver( const GetPot& input )
    : Solver(input),
      _n_timesteps( input("unsteady-solver/n_timesteps", 1 ) ),
      _deltat( input("unsteady-solver/deltat", 0.0 ) )
  {
    if( !input.have_variable("unsteady-solver/n_timesteps") )
      libmesh_error_msg("ERROR: Must specify input parameter unsteady-solver/n_timesteps!");

    if( !input.have_variable("unsteady-solver/deltat") )
      libmesh_error_msg("ERROR: Must specify input parameter unsteady-solver/deltat!");
  }

  void NewmarkSolver::init_time_solver(MultiphysicsSystem* system)
  {
    system->time_solver = libMesh::AutoPtr<libMesh::TimeSolver>(new libMesh::NewmarkSolver(*system) );
  }

  void NewmarkSolver::solve( SolverContext& context )
  {
    libmesh_assert( context.system );

    context.system->deltat = this->_deltat;

    libMesh::NewmarkSolver * newmark = libMesh::cast_ptr<libMesh::NewmarkSolver *>(context.system->time_solver.get());

    newmark->compute_initial_accel();

    libMesh::Real sim_time;

    if( context.output_vis )
      {
	context.postprocessing->update_quantities( *(context.equation_system) );
	context.vis->output( context.equation_system );
      }

    std::time_t first_wall_time = std::time(NULL);

    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step < this->_n_timesteps; t_step++)
      {
        std::time_t latest_wall_time = std::time(NULL);

	std::cout << "==========================================================" << std::endl
		  << "   Beginning time step " << t_step  <<
                     ", t = " << context.system->time <<
                     ", dt = " << context.system->deltat <<
                     ", runtime = " << (latest_wall_time - first_wall_time) <<
                     std::endl
		  << "==========================================================" << std::endl;

        context.system->solve();

	sim_time = context.system->time;

	if( context.output_vis && !((t_step+1)%context.timesteps_per_vis) )
	  {
	    context.postprocessing->update_quantities( *(context.equation_system) );
	    context.vis->output( context.equation_system, t_step, sim_time );
	  }

	if( context.output_residual && !((t_step+1)%context.timesteps_per_vis) )
	  context.vis->output_residual( context.equation_system, context.system,
                                        t_step, sim_time );

        if ( context.print_perflog && context.timesteps_per_perflog
             && !((t_step+1)%context.timesteps_per_perflog) )
          libMesh::perflog.print_log();

        // Advance to the next timestep
	context.system->time_solver->advance_timestep();
      }

    std::time_t final_wall_time = std::time(NULL);
    std::cout << "==========================================================" << std::endl
	      << "   Ending time stepping, t = " << context.system->time <<
                 ", runtime = " << (final_wall_time - first_wall_time) <<
                 std::endl
              << "==========================================================" << std::endl;
  }

} // end namespace GRINS
