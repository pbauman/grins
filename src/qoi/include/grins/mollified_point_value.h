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


#ifndef GRINS_MOLLIFIED_POINT_VALUE_H
#define GRINS_MOLLIFIED_POINT_VALUE_H

// GRINS
#include "grins/qoi_base.h"

// libMesh
#include "libmesh/point.h"

namespace GRINS
{
  // GRINS forward declarations
  class AssemblyContext;
  class MultiphysicsSystem;

  //! Point value QoI using mollification
  /*! This QoI follows the work of Prudhomme, Oden, "On goal-oriented error
    estimation for elliptic problems: application to the control of pointwise
    errors", CMAME, v.176, 1999. In particular, we take

    \f[ Q(u) = \int_{\Omega} u(x_0) k_{\epsilon}(x-x_0)\;dx \f]

    where \f$ k_{\epsilon}(x) \f$ are mollification functions

    \f[ k_{\epsilon}(x) =
        \begin{cases}
           C \exp{-\frac{\epsilon^2}{\epsilon^2 - \|x\|^2}} & \|x\| < \epsilon \\
           0                                                & \|x\| >= \epsilon
        \end{cases} \f]. The constant \f$ C \f$ is determined by ensuring that
    \f$ \int_{\Omega} k_{\epsilon}(x-x_0)\;dx = 1 \f$

    We follow the suggestion of Prudhomme/Oden to fix the ratio \f$\kappa \le \frac{2\epsilon}{h}\f$
    where \f$ h \f$ is the element size. So, as the mesh is refined, so is \f$ \epsilon \f$.
    Thus, the user sets \f$ \kappa \f$ from which $\f \epsilon \f$ and \f$ C \f$ are computed.
    The default for \f$ \kappa \f$ is 0.25, as suggested by Prudhomme/Oden.
    Because of the form of the mollification function, we allow the user
    to specify a higher order quadrature for the integration of this QoI.
  */
  class MollifiedPointValue : public QoIBase
  {
  public:

    MollifiedPointValue( const std::string& qoi_name );

    virtual ~MollifiedPointValue();

    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void interior_qoi( AssemblyContext& context,
                               const unsigned int qoi_index );

    virtual void interior_qoi_derivative( AssemblyContext& context,
                                          const unsigned int qoi_index );

    virtual void init( const GetPot& input, const MultiphysicsSystem& system );

    virtual void init_context( AssemblyContext& context );

    //! Retrieve precomputed integration constant
    /*! Mainly for testing, unlikely this is needed by the user. */
    libMesh::Real get_int_constant() const;

  protected:

    //! Compute value of integration constant
    /*! Formula depends on the dimension. Used to compute \f$ C \f$.
        In parcticular, we return the INVERSE of
        \f$ \int_{-1}^1 \exp{ \frac{-1}{1-x^2} } \; dx \f$ in 1D and
        the analogs for 2D and 3D. */
    libMesh::Real compute_integration( unsigned int dim ) const;

    //! Mollifying function. User supplies \f$ C \f$ and \f$ \|\mathbf{x}\|^2.
    /*! \f[ k_{\epsilon}(\mathbf{x}) =
            \begin{cases}
            C \exp{\frac{-\epsilon^2}{\epsilon^2 - \|\mathbf{x}\|^2}} & \|\mathbf{x}\| < \epsilon \\
            0 & \|\mathbf{x}\| \ge \epsilon
            \end{cases}  \f]
    */
    libMesh::Real mollification_function( libMesh::Real C, libMesh::Real eps, libMesh::Real norm_x_sq ) const;

    //! Point at which we want to compute value of variable.
    libMesh::Point _point;

    //! Variable that we're computing.
    VariableIndex _var;

    //! Used to determine \f$ \epsilon \f$.
    /*! \f$ \kappa \le \frac{2\epsilon}{h} \f$, \f$ h \f$ is the element size.
        \f$ \kappa \f$ should be at most 1. The recommended, and default, value
        is 0.25. This balances accuracy (\f$ \epsilon \rightarrow 0 \f$) while
        keeping number of quadrature points limited. */
    libMesh::Real _kappa;

    //! This is the integral we can precompute to compute \f$ C \f$.
    /*! In particular, \f$ C = (_int_constant)*\epsilon^{-d} \f$, where
        \f$ d \f$ is the dimension. \f$ \epsilon \f$ is computed once the
        element size is known. */
    libMesh::Real _int_constant;

  };

  inline
  bool MollifiedPointValue::assemble_on_interior() const
  {
    return true;
  }

  inline
  bool MollifiedPointValue::assemble_on_sides() const
  {
    return false;
  }

  inline
  libMesh::Real MollifiedPointValue::get_int_constant() const
  {
    return _int_constant;
  }

} // end namespace GRINS

#endif // GRINS_MOLLIFIED_POINT_VALUE_H
