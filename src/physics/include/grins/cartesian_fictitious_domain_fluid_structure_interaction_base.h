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

#ifndef GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H
#define GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H

// GRINS
#include "grins/fictitious_domain_fluid_structure_interaction_abstract.h"

namespace GRINS
{
  template<unsigned int Dim, bool UseOldDisplacement>
  class CartesianFictitiousDomainFluidStructureInteractionBase :
    public FictitiousDomainFluidStructureInteractionAbstract
  {
  public:

    CartesianFictitiousDomainFluidStructureInteractionBase( const PhysicsName & physics_name,
                                                            const GetPot & input );

    CartesianFictitiousDomainFluidStructureInteractionBase() = delete;

    virtual ~CartesianFictitiousDomainFluidStructureInteractionBase() = default;

    //! Context initialization
    virtual void init_context( AssemblyContext & context ) override;

  protected:

    void check_variable_dim_consistency() const;

    std::string dim_error_msg(unsigned int var_dim) const;

    //! 2D deformation gradient for plane strain
    /* F(2,2) = 1 will be consistent with plane strain for everything that consumes F */
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v ) const;

    //! 3D deformation gradient
    libMesh::Tensor form_def_gradient( const libMesh::Gradient & grad_u,
                                       const libMesh::Gradient & grad_v,
                                       const libMesh::Gradient & grad_w ) const;

  };

} // end namespace GRINS

#endif // GRINS_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_BASE_H
