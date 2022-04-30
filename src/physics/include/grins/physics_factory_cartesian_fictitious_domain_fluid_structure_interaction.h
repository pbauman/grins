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

#ifndef GRINS_PHYSICS_FACTORY_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H
#define GRINS_PHYSICS_FACTORY_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H

// GRINS
#include "grins/materials_parsing.h"
#include "grins/physics_factory_with_core.h"
#include "grins/variable_warehouse.h"
#include "grins/multi_component_vector_variable.h"

namespace GRINS
{
  template<template<unsigned int> class DerivedPhysics>
  class PhysicsFactoryCartesianFictitiousDomainFluidStructureInteraction : public PhysicsFactoryWithCore
  {
  public:

    PhysicsFactoryCartesianFictitiousDomainFluidStructureInteraction( const std::string & physics_name,
                                                                      const std::string & core_physics_name )
      : PhysicsFactoryWithCore(physics_name,core_physics_name)
    {}

    PhysicsFactoryCartesianFictitiousDomainFluidStructureInteraction() = delete;

    virtual ~PhysicsFactoryCartesianFictitiousDomainFluidStructureInteraction() = default;

  protected:

    virtual std::unique_ptr<Physics> build_physics( const GetPot & input,
                                                    const std::string & physics_name );

  };

  template<template<unsigned int> class DerivedPhysics>
  inline
  std::unique_ptr<Physics>
  PhysicsFactoryCartesianFictitiousDomainFluidStructureInteraction<DerivedPhysics>::build_physics
  ( const GetPot & input, const std::string & physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::unique_ptr<Physics> new_physics;

    // Deduce the dimension from the number displacment variable components
    const DisplacementVariable & disp_var =
      GRINSPrivate::VariableWarehouse::get_variable_subclass<DisplacementVariable>
      (VariablesParsing::disp_variable_name(input,physics_name,VariablesParsing::PHYSICS));

    int dim = disp_var.dim();

    if(dim==2)
      new_physics = std::make_unique<DerivedPhysics<2>>(physics_name,input);
    else if(dim==3)
      new_physics = std::make_unique<DerivedPhysics<3>>(physics_name,input);
    else
      libmesh_error_msg("ERROR: Cartesian Fictitious Domain FSI is only valid for dimensions 2 or 3!\n");

    libmesh_assert(new_physics);

    return new_physics;
  }

} // end namespace GRINS

#endif // GRINS_PHYSICS_FACTORY_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H
