//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/physics_factory_immersed_boundary.h"

// GRINS
#include "grins/physics_factory_helper.h"
#include "grins/immersed_boundary.h"

// Solid Mechanics headers
#include "grins/elastic_cable.h"
#include "grins/elastic_membrane.h"
#include "grins/elastic_cable_rayleigh_damping.h"

// Stress Strain law headers
#include "grins/hookes_law.h"
#include "grins/hookes_law_1d.h"
#include "grins/hyperelasticity.h"
#include "grins/mooney_rivlin.h"
#include "grins/incompressible_plane_stress_hyperelasticity.h"

namespace GRINS
{

  libMesh::UniquePtr<Physics> PhysicsFactoryImmersedBoundary::build_physics( const GetPot& input,
                                                                             const std::string& physics_name )
  {

    std::cout << "Factory : Building Ibm Physics "<< std::endl;
    std::string solid_mech_input = "none";
    PhysicsFactoryHelper::parse_immersed_boundary_components( input,
                                                              physics_name,
                                                              solid_mech_input);

    //Check to see if input solid mechanics are pre-existing
    std::map< std::string, FactoryAbstract< Physics > * > & existing_factories  = this->factory_map();
    std::map< std::string, FactoryAbstract< Physics > * >::iterator factory_it;

    factory_it = existing_factories.find(solid_mech_input);
    if (factory_it == existing_factories.end())
      {
        std::string error = "Error: Invalid solid_mechanics specified: " + solid_mech_input + '\n';
        error += "       Valid values are: ElasticMembrane, ElasticCable";
        libmesh_error_msg(error);
      }


    //Now parse the supplied Solid Mechanics for valid stress strain models
    std::string stress_strain_model = "none";
    std::string strain_energy = "none";
    PhysicsFactoryHelper::parse_stress_strain_model( input,
                                                     solid_mech_input,
                                                     stress_strain_model,
                                                     strain_energy);

    libMesh::UniquePtr<Physics> new_physics;

    if (solid_mech_input == std::string("ElasticCable"))
      {
        if( stress_strain_model == std::string("hookes_law") )
          {
            libMesh::UniquePtr<ElasticCable<HookesLaw1D>> solid_mech_ptr( new ElasticCable<HookesLaw1D>(solid_mech_input,input,false) );
            new_physics.reset( new ImmersedBoundary<ElasticCable<HookesLaw1D> >(physics_name, solid_mech_ptr,  input) );
          }
        else
          {
            std::string error = "Error: Invalid stress-strain model: "+stress_strain_model+" for " + solid_mech_input + "!\n";
            error += "       Valid values are: hookes_law\n";
            libmesh_error_msg(error);
          }
      }
    else if (solid_mech_input == std::string("ElasticMembrane"))
      {
        if( stress_strain_model == std::string("hookes_law") )
          {
            libMesh::UniquePtr<ElasticMembrane<HookesLaw>> solid_mech_ptr( new ElasticMembrane<HookesLaw>(solid_mech_input,input,false) );
            new_physics.reset( new ImmersedBoundary< ElasticMembrane<HookesLaw> >(physics_name, solid_mech_ptr, input) );
          }

        else if( stress_strain_model == std::string("incompressible_hyperelasticity") )
          {
            if( strain_energy == std::string("mooney_rivlin") )
              {
                libMesh::UniquePtr< ElasticMembrane< IncompressiblePlaneStressHyperelasticity< MooneyRivlin>>> solid_mech_ptr( new ElasticMembrane< IncompressiblePlaneStressHyperelasticity< MooneyRivlin>>(solid_mech_input,input,false));
                new_physics.reset( new ImmersedBoundary< ElasticMembrane< IncompressiblePlaneStressHyperelasticity<MooneyRivlin>>>(physics_name, solid_mech_ptr, input) );
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
            std::string error = "Error: Invalid stress-strain model: "+stress_strain_model+"!\n";
            error += "       Valid values are: hookes_law\n";
            error += "                         incompressible_hyperelasticity\n";
            libmesh_error_msg(error);
          }
      }

    std::cout << "Factory : Building Ibm Physics Complete "<< std::endl;

    libmesh_assert(new_physics);
    return new_physics;

  } //end class PhysicsFactoryImmersedBoundary

  // Instantiate the immersed boundary factory
  PhysicsFactoryImmersedBoundary grins_factory_immersed_boundary(PhysicsNaming::immersed_boundary());

} // end namespace GRINS
