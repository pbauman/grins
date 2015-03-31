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

// GRINS
#include "grins/mollified_point_value.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"

void init_dummy_system( const GetPot& input,
                        libMesh::EquationSystems& es,
                        GRINS::MollifiedPointValue& qoi );

libMesh::Real compute_1d_quad( GRINS::MollifiedPointValue& qoi );
libMesh::Real compute_2d_quad( GRINS::MollifiedPointValue& qoi );
libMesh::Real compute_3d_quad( GRINS::MollifiedPointValue& qoi );

int main( int argc, char* argv[] )
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );

  // Initialize libMesh and any dependent libaries, like in example 2.
  libMesh::LibMeshInit init (argc, argv);

  GRINS::MollifiedPointValue qoi( "Test" );

  libMesh::Mesh mesh(init.comm());

  int return_flag = 0;

  // First do 1D
  {
    libMesh::MeshTools::Generation::build_line(mesh,10,-1.0, 1.0,libMesh::EDGE2);
    libMesh::EquationSystems es (mesh);
    init_dummy_system(input,es,qoi);

    // Check that we integrate 1D mollifier function (approximately) correctly
    libMesh::Real value = compute_1d_quad(qoi);
    libMesh::Real tol = 1.0e-6; // Seems to be what we're getting for 1D
    libMesh::Real error = value - 1.0;
    if( std::fabs(error) > tol )
      {
        std::cerr << "Mismatch on 1D quadrature test!" << std::endl
                  << "error = " << error << std::endl
                  << "tol = " << tol << std::endl;
        return_flag = 1;
      }
  }


  // Now 2D
  {
    libMesh::MeshTools::Generation::build_square(mesh,10,10,-1.0, 1.0,-1.0,1.0,libMesh::QUAD4);
    libMesh::EquationSystems es (mesh);
    init_dummy_system( input, es, qoi );

    // Check that we integrate 2D mollifier function (approximately) correctly
    libMesh::Real value = compute_2d_quad(qoi);
    libMesh::Real tol = 1.0e-4; // Seems to be what we're getting for 2D
    libMesh::Real error = value - 1.0;
    if( std::fabs(error) > tol )
      {
        std::cerr << "Mismatch on 2D quadrature test!" << std::endl
                  << "error = " << error << std::endl
                  << "tol = " << tol << std::endl;
        return_flag = 1;
      }
  }


  // Now 3D
  {
    libMesh::MeshTools::Generation::build_cube(mesh,10,10,10,-1.0, 1.0,-1.0,1.0,-1.0,1.0,libMesh::HEX8);
    libMesh::EquationSystems es (mesh);
    init_dummy_system( input, es, qoi );

    // Check that we integrate 3D mollifier function (approximately) correctly
    libMesh::Real value = compute_3d_quad(qoi);
    libMesh::Real tol = 1.0e-4; // Seems to be what we're getting for 3D
    libMesh::Real error = value - 1.0;
    if( std::fabs(error) > tol )
      {
        std::cerr << "Mismatch on 3D quadrature test!" << std::endl
                  << "error = " << error << std::endl
                  << "tol = " << tol << std::endl;
        return_flag = 1;
      }
  }



  return return_flag;
}

void init_dummy_system( const GetPot& input,
                        libMesh::EquationSystems& es,
                        GRINS::MollifiedPointValue& qoi )
{
  GRINS::MultiphysicsSystem& system = es.add_system<GRINS::MultiphysicsSystem> ("TestQoI");
  system.add_variable("u");

  qoi.init(input,system);

  return;
}

libMesh::Real compute_1d_quad( GRINS::MollifiedPointValue& qoi )
{
  libMesh::QGauss qrule(1,libMesh::Order::FORTYTHIRD);
  const std::vector<libMesh::Real>& weights = qrule.get_weights();
  const std::vector<libMesh::Point>& points = qrule.get_points();
  const unsigned int n_qpoints = weights.size();

  libMesh::Real eps = 0.25;
  libMesh::Real C = qoi.get_constant(eps);

  libMesh::Real value = 0.0;

  for( unsigned int q = 0; q < n_qpoints; q++ )
    {
      /* This integral is from 0 to 1 so we need to map the quadrature points.
         \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f( (b-a)/2*z + (a+b)/2 ) dz */
      libMesh::Real a = -eps;
      libMesh::Real b = eps;
      libMesh::Real w_qp = (b-a)/2.0*weights[q];
      libMesh::Real x_qp = (a+b)/2.0 + (b-a)/2.0*points[q](0);
      libMesh::Real r2 = x_qp*x_qp;

      value += w_qp*qoi.mollification_function(C,eps,r2);
    }

  return value;
}

libMesh::Real compute_2d_quad( GRINS::MollifiedPointValue& qoi )
{
  libMesh::QGauss qrule(2,libMesh::Order::FORTYTHIRD);
  qrule.init(libMesh::QUAD4);
  const std::vector<libMesh::Real>& weights = qrule.get_weights();
  const std::vector<libMesh::Point>& points = qrule.get_points();
  const unsigned int n_qpoints = weights.size();

  libMesh::Real eps = 0.25;
  libMesh::Real C = qoi.get_constant(eps);

  libMesh::Real value = 0.0;

  for( unsigned int q = 0; q < n_qpoints; q++ )
    {
      /* This integral is from 0 to 1 so we need to map the quadrature points.
         \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f( (b-a)/2*z + (a+b)/2 ) dz */
      libMesh::Real a = -eps;
      libMesh::Real b = eps;
      libMesh::Real w_qp = (b-a)/2.0*(b-a)/2.0*weights[q];
      libMesh::Real x_qp = (a+b)/2.0 + (b-a)/2.0*points[q](0);
      libMesh::Real y_qp = (a+b)/2.0 + (b-a)/2.0*points[q](1);
      libMesh::Real r2 = x_qp*x_qp + y_qp*y_qp;

      value += w_qp*qoi.mollification_function(C,eps,r2);
    }

  return value;
}

libMesh::Real compute_3d_quad( GRINS::MollifiedPointValue& qoi )
{
  libMesh::QGauss qrule(3,libMesh::Order::FORTYTHIRD);
  qrule.init(libMesh::HEX8);
  const std::vector<libMesh::Real>& weights = qrule.get_weights();
  const std::vector<libMesh::Point>& points = qrule.get_points();
  const unsigned int n_qpoints = weights.size();

  libMesh::Real eps = 0.25;
  libMesh::Real C = qoi.get_constant(eps);

  libMesh::Real value = 0.0;

  for( unsigned int q = 0; q < n_qpoints; q++ )
    {
      /* This integral is from 0 to 1 so we need to map the quadrature points.
         \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f( (b-a)/2*z + (a+b)/2 ) dz */
      libMesh::Real a = -eps;
      libMesh::Real b = eps;
      libMesh::Real w_qp = (b-a)/2.0*(b-a)/2.0*(b-a)/2.0*weights[q];
      libMesh::Real x_qp = (a+b)/2.0 + (b-a)/2.0*points[q](0);
      libMesh::Real y_qp = (a+b)/2.0 + (b-a)/2.0*points[q](1);
      libMesh::Real z_qp = (a+b)/2.0 + (b-a)/2.0*points[q](2);
      libMesh::Real r2 = x_qp*x_qp + y_qp*y_qp + z_qp*z_qp;

      value += w_qp*qoi.mollification_function(C,eps,r2);
    }

  return value;
}
