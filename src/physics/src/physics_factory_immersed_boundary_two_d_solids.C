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
#include "grins/physics_factory_immersed_boundary_two_d_solids.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/immersed_boundary.h"

// Solid Mechanics headers
#include "grins/elastic_membrane.h"
#include "grins/elastic_membrane_rayleigh_damping.h"

// Stress Strain law headers
#include "grins/hookes_law.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"
#include "grins/mooney_rivlin.h"
#include "grins/hyperelasticity.h"



namespace GRINS
{

  template<template<typename> class DerivedPhysics>
  libMesh::UniquePtr<Physics> PhysicsFactoryImmersedBoundaryTwoDSolids<DerivedPhysics>::build_physics( const GetPot& input,
                                                                                                       const std::string& physics_name )
  {
    std::string core_physics = this->find_core_physics_name(physics_name);

    std::string solid_mech = "none";
    PhysicsFactoryHelper::parse_immersed_boundary_components( input,
                                                              core_physics,
                                                              solid_mech);

    //Check to see if input solid mechanics are pre-existing
    std::map< std::string, FactoryAbstract< Physics > * > & existing_factories  = this->factory_map();
    std::map< std::string, FactoryAbstract< Physics > * >::iterator factory_it;

    factory_it = existing_factories.find(solid_mech);
    if (factory_it == existing_factories.end())
      {
        std::string error = "Error: Invalid solid_mechanics specified:  " + solid_mech+ '\n';
        error += "       Valid values are: ElasticMembrane, ElasticCable"; 
          libmesh_error_msg(error);
        }
    
    //Now parse the supplied Solid Mechanics for valid stress strain models
    std::string stress_strain_model = "none";
    std::string strain_energy = "none";
    PhysicsFactoryHelper::parse_stress_strain_model( input,
                                                     core_physics,
                                                     stress_strain_model,
                                                     strain_energy);


    libMesh::UniquePtr<Physics> new_physics;

    if (solid_mech == std::string("ElasticMembrane"))
        {
          if( stress_strain_model == std::string("hookes_law") )
            new_physics.reset( new DerivedPhysics<HookesLaw>
                               (physics_name,input, false /*is_compressible*/) );

          else if( stress_strain_model == std::string("incompressible_hyperelasticity") )
            {
              if( strain_energy == std::string("mooney_rivlin") )
                {
                  new_physics.reset( new DerivedPhysics<IncompressiblePlaneStressHyperelasticity<MooneyRivlin> >(physics_name,input,false /*is_compressible*/) );
                }
              else
                {
                  std::string error = "ERROR: Invalid strain_energy "+strain_energy+"!\n";
                  error += "       Valid values are: mooney_rivlin\n";
                  libmesh_error_msg(error);
                }
            }
          else
            {
              std::string error = "Error: Invalid stress-strain model: "+ stress_strain_model +"!\n";
              error += "       Valid values are: hookes_law\n";
              error += "                         incompressible_hyperelasticity\n";
              libmesh_error_msg(error);
            }
        }

    libmesh_assert(new_physics);
    return new_physics;
    
  }//end PhysicsFactoryImmersedBoundaryOneDSolids
  
  // Instantiate all the immersed boundary factories.
  PhysicsFactoryImmersedBoundaryTwoDSolids<ElasticMembrane> grins_factory_ibm_elastic_membrane
  (PhysicsNaming::elastic_membrane(),PhysicsNaming::immersed_boundary());
  
  PhysicsFactoryImmersedBoundaryTwoDSolids<ElasticMembraneRayleighDamping> grins_factory_ibm_elastic_membrane_rayleigh_damping
  (PhysicsNaming::elastic_membrane_rayleigh_damping(),PhysicsNaming::immersed_boundary());
    
} // end namespace GRINS
