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
#include "grins/mollified_point_value.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature_gauss.h"

namespace GRINS
{
  MollifiedPointValue::MollifiedPointValue( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  MollifiedPointValue::~MollifiedPointValue()
  {
    return;
  }

  QoIBase* MollifiedPointValue::clone() const
  {
    return new MollifiedPointValue( *this );
  }

  void MollifiedPointValue::init( const GetPot& input, const MultiphysicsSystem& system )
  {
    // What variable are we computing the point value of?
    if( !input.have_variable("QoI/MollifiedPointValue/variable_name") )
      {
        libmesh_error_msg("ERROR: Must specify QoI/MollifiedPointValue/variable_name.");
      }
    std::string var_name = input("QoI/MollifiedPointValue/variable_name", "DIE!");
    _var = system.variable_number(var_name);

    // What point are we computing the value of _var?
    unsigned int dim = system.get_mesh().mesh_dimension();
    if( input.vector_variable_size("QoI/MollifiedPointValue/point") < dim )
      {
        libmesh_error_msg("ERROR: Must specify as many components as mesh_dimension() for QoI/MollifiedPointValue/point");
      }

    for( unsigned int i = 0; i < dim; i++ )
      {
        (_point)(i) = input("QoI/MollifiedPointValue/point", 0.0, i );
      }

    // Grab \kappa. Default to \kappa = 0.25
    _kappa = input("QoI/MollifiedPointValue/kappa", 0.25);

    /* Now compute coefficient C such that
       \int_{\Omega} k_{\epsilon}(x-x_0)\;dx = 1.
       We just use a higher order Gauss rule to do the integral since,
       using a similar argument to Prudhomme/Oden */
    _int_constant = this->compute_integration(dim);

    return;
  }

  void MollifiedPointValue::init_context( AssemblyContext& context )
  {
    return;
  }

  void MollifiedPointValue::interior_qoi( AssemblyContext& context,
                                          const unsigned int qoi_index )
  {
    return;
  }

  void MollifiedPointValue::interior_qoi_derivative( AssemblyContext& context,
                                                     const unsigned int qoi_index )
  {
    return;
  }

  libMesh::Real MollifiedPointValue::compute_integration( unsigned int dim ) const
  {
    libMesh::Real value = 0.0;

    // All the integrals are 1D. We use a 15th order Gauss rule `cause why not.
    libMesh::QGauss qrule(1,libMesh::Order::FIFTEENTH);

    const std::vector<libMesh::Real>& weights = qrule.get_weights();
    const std::vector<libMesh::Point>& points = qrule.get_points();
    const unsigned int n_qpoints = weights.size();
    libmesh_assert_equal_to( weights.size(), points.size() );

    switch(dim)
      {
        /* 1D:
           \int_{\Omega} k_{\epsilon}(x-x_0)\;dx =
           \int_{x_0 - \epsilon}^{x_0 + \epsilon} k_{\epsilon}(x-x_0)\;dx =
           C \epsilon \int_{-1}^1 \exp{ \frac{-1}{1-x^2} } \; dx
        */
      case(1):
        {
          for( unsigned int q = 0; q < n_qpoints; q++ )
            {
              // This integral is from -1 to 1, so we don't need to map the quadrature points.
              libMesh::Real w_qp = weights[q];
              libMesh::Real x_qp = points[q](0);

              value += w_qp*std::exp( (-1.0)/(1 - x_qp*x_qp) );
            }
        }
        break;

        /* 2D:
           \int_{\Omega} k_{\epsilon}(x-x_0)\;dx =
           \int_{0}^{2\pi} \int_{0}^{\epsilon} k_{\epsilon}(r)\; r dr d\theta =
           2\pi \epsilon^2 \int_0^1 k_{\epsilon}(r)\; r dr =
           C 2\pi \epsilon^2 \int_0^1 \exp{\frac{-1}{1-r^2}} r dr
        */
      case(2):
        {
          for( unsigned int q = 0; q < n_qpoints; q++ )
            {
              /* This integral is from 0 to 1 so we need to map the quadrature points.
                 \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f( (b-a)/2*z + (a+b)/2 ) dz
                 Here: a = 0, b = 1.
                 So scale weights by 0.5 and map points accordingly. */
              libMesh::Real w_qp = 0.5*weights[q];
              libMesh::Real x_qp = 0.5 + 0.5*points[q](0);

              value += w_qp*x_qp*std::exp( (-1.0)/(1 - x_qp*x_qp) );
            }
          value *= Constants::two_pi;
        }
        break;

        /* 3D:
           \int_{\Omega} k_{\epsilon}(x-x_0)\;dx =
           \int_{0}^{2\pi} \int_0^{\pi} \int_{0}^{\epsilon} k_{\epsilon}(r)\; r^2 \sin(\phi) dr d\phi d\theta =
           4\pi \epsilon^3 \int_0^1 k_{\epsilon}(r)\; r^2 dr =
           C 4\pi \epsilon^3 \int_0^1 \exp{\frac{-1}{1-r^2}} r^2 dr
        */
      case(3):
        {
          for( unsigned int q = 0; q < n_qpoints; q++ )
            {
              /* This integral is from 0 to 1 so we need to map the quadrature points.
                 \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f( (b-a)/2*z + (a+b)/2 ) dz
                 Here: a = 0, b = 1.
                 So scale weights by 0.5 and map points accordingly. */
              libMesh::Real w_qp = 0.5*weights[q];
              libMesh::Real x_qp = 0.5 + 0.5*points[q](0);
              libMesh::Real r2 = x_qp*x_qp;

              value += w_qp*r2*std::exp( (-1.0)/(1.0-r2) );
            }
          value *= 2.0*Constants::two_pi;
        }
        break;

      default:
        {
          // This shouldn't happen
          libmesh_error();
        }

      } // switch(dim)

    // We'd better have computed a non-zero integral
    libmesh_assert_greater(value,0.0);

    // We return the inverse since that is what is desired in the end.
    return 1.0/value;
  }

  libMesh::Real MollifiedPointValue::mollification_function( libMesh::Real C, libMesh::Real eps,
                                                             libMesh::Real norm_x_sq ) const
  {
    libMesh::Real value = 0.0;

    if( std::sqrt(norm_x_sq) < eps )
      {
        libMesh::Real eps2 = eps*eps;
        value = C*std::exp( (-eps2)/(eps2 - norm_x_sq) );
      }

    return value;
  }

} // end namespace GRINS
