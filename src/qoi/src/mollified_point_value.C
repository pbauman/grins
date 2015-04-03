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
    _dim = system.get_mesh().mesh_dimension();
    if( input.vector_variable_size("QoI/MollifiedPointValue/point") < _dim )
      {
        libmesh_error_msg("ERROR: Must specify as many components as mesh_dimension() for QoI/MollifiedPointValue/point");
      }

    for( unsigned int i = 0; i < _dim; i++ )
      {
        (_point)(i) = input("QoI/MollifiedPointValue/point", 0.0, i );
      }

    // Grab \kappa. Default to \kappa = 0.25
    _kappa = input("QoI/MollifiedPointValue/kappa", 0.25);

    /* Now compute coefficient needed for constant C such that
       \int_{\Omega} k_{\epsilon}(x-x_0)\;dx = 1. */
    _int_constant = this->compute_integration(_dim);

    return;
  }

  void MollifiedPointValue::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* element_fe;
    context.get_element_fe<libMesh::Real>(_var, element_fe);
    element_fe->get_JxW();
    element_fe->get_xyz();
    element_fe->get_phi();

    return;
  }

  void MollifiedPointValue::element_qoi( AssemblyContext& context,
                                         const unsigned int qoi_index )
  {
    if( context.get_elem().contains_point( _point ) )
      {
        /* We use hmin to ensure we keep the mollifying function
           within the element */
        libMesh::Real h = context.get_elem().hmin();

        // Back out eps based on the mesh size (and input kappa)
        libMesh::Real eps = this->compute_eps( h );

        // Now compute constant C for mollifying function
        libMesh::Real C = this->get_constant(eps);

        libMesh::FEBase* element_fe;
        context.get_element_fe<libMesh::Real>(_var, element_fe);

        const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();
        const std::vector<libMesh::Point>& x_qp = element_fe->get_xyz();

        unsigned int n_qpoints = context.get_element_qrule().n_points();

        libMesh::Number& qoi = context.get_qois()[qoi_index];

        for (unsigned int qp = 0; qp != n_qpoints; qp++)
          {
            // Get the solution value at the quadrature point
            libMesh::Real u;
            context.interior_value(_var,qp,u);

            libMesh::Point xmx0 = x_qp[qp] - _point;
            libMesh::Real norm_x_sq = xmx0*xmx0;

            libMesh::Real k_eps = this->mollification_function(C,eps,norm_x_sq);

            qoi += u*k_eps*JxW[qp];
          }

      } // contains point

    return;
  }

  void MollifiedPointValue::element_qoi_derivative( AssemblyContext& context,
                                                    const unsigned int qoi_index )
  {
    if( context.get_elem().contains_point( _point ) )
      {
        /* We use hmin to ensure we keep the mollifying function
           within the element */
        libMesh::Real h = context.get_elem().hmin();

        // Back out eps based on the mesh size (and input kappa)
        libMesh::Real eps = this->compute_eps( h );

        // Now compute constant C for mollifying function
        libMesh::Real C = this->get_constant(eps);

        libMesh::FEBase* element_fe;
        context.get_element_fe<libMesh::Real>(_var, element_fe);

        const std::vector<libMesh::Real> &JxW = element_fe->get_JxW();
        const std::vector<libMesh::Point>& x_qp = element_fe->get_xyz();

        unsigned int n_qpoints = context.get_element_qrule().n_points();

        const unsigned int n_dofs = context.get_dof_indices(_var).size();

        const std::vector<std::vector<libMesh::Real> >& u_phi = element_fe->get_phi();

        libMesh::DenseSubVector<libMesh::Number>& dQ_du =
              context.get_qoi_derivatives(qoi_index, _var);

        for (unsigned int qp = 0; qp != n_qpoints; qp++)
          {
            libMesh::Point xmx0 = x_qp[qp] - _point;
            libMesh::Real norm_x_sq = xmx0*xmx0;

            libMesh::Real k_eps = this->mollification_function(C,eps,norm_x_sq);

            for( unsigned int i = 0; i != n_dofs; i++ )
              {
                dQ_du(i) += u_phi[i][qp]*k_eps*JxW[qp];
              }
          }

      } // contains point

    return;
  }

  libMesh::Real MollifiedPointValue::compute_integration( unsigned int dim ) const
  {
    libMesh::Real value = 0.0;

    switch(dim)
      {
        /* 1D:
           \int_{\Omega} k_{\epsilon}(x-x_0)\;dx =
           \int_{x_0 - \epsilon}^{x_0 + \epsilon} k_{\epsilon}(x-x_0)\;dx =
           C \epsilon \int_{-1}^1 \exp{ \frac{-1}{1-x^2} } \; dx
        */
      case(1):
        {
          // Generated from MATLAB using the following commands:
          // f = @(x) exp( (-1.0)./(1.0-x.^2) )
          // q1d = quad(f, -1.0, 1.0, 1.0e-13, 1)
          value = 0.443993816168134;
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
          // Generated from MATLAB using the following commands:
          // f = @(x) x.*exp( (-1.0)./(1.0-x.^2) )
          // q2d = 2*pi*quad(f, 0.0, 1.0, 1.0e-13, 1)
          value = 0.466512393178512;
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
          // Generated from MATLAB using the following commands:
          // f = @(x) x.^2.*exp( (-1.0)./(1.0-x.^2) )
          // q3d = 4*pi*quad(f, 0.0, 1.0, 1.0e-13, 1)
          value = 0.441088887276947;
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

  libMesh::Real MollifiedPointValue::get_constant(libMesh::Real eps) const
  {
    libMesh::Real value = 0.0;
    switch(_dim)
      {
      case(1):
        value = this->_int_constant/eps;
        break;
      case(2):
        value = this->_int_constant/(eps*eps);
        break;
      case(3):
        value = this->_int_constant/(eps*eps*eps);
        break;
      default:
        // Wat
        libmesh_error();
      }

    // We'd better have a positive constant
    libmesh_assert_greater(value,0.0);
    return value;
  }

} // end namespace GRINS
