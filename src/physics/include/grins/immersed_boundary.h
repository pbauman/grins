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


#ifndef GRINS_IMMERSED_BOUNDARY_H
#define GRINS_IMMERSED_BOUNDARY_H

//GRINS
#include "grins/physics.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/common.h"
#include "grins/solid_mechanics_abstract.h"

//libMesh
#include "libmesh/fem_context.h"


namespace GRINS
{

  //! Physics class for the Immersed Boundary Method
  /*!
    This physics class implements the classical FEM Immersed  Boundary Method.
    This is a templated class, the class SolidMech can be instantiated as a specific type
    (right now as: ElasticCable or ElasticMembrane)
   */

  template<typename SolidMech>
  class ImmersedBoundary : public Physics
  {
  public:

    ImmersedBoundary( const std::string & my_physics_name,
                      libMesh::UniquePtr<SolidMech> & solid_mech_ptr,
                      const GetPot& input );

    virtual ~ImmersedBoundary(){};

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Init variables and stuff we need
    virtual void init_variables( libMesh::FEMSystem* system);

    //! Context initialization
    virtual void init_context( AssemblyContext& context );

    //! Residual contributions from the solid in the flow
    virtual void element_time_derivative( bool compute_jacobian, AssemblyContext& context,
                                          CachedValues& cache );

    //! Cache mesh information needed for residual computation
    virtual void preassembly( MultiphysicsSystem & system );

  private:

    //! FE variables for the flow
    VelocityVariable & _flow_vars;

    //! FE variables for the solid
    DisplacementVariable & _disp_vars;

    //! Solid Mechanics from the ibm factory
    libMesh::UniquePtr<SolidMech> _solid_mech;

    //! The fluid mechanics associated with the IBM method from the input
    std::string _fluid_mechanics;

    //! The solid mechanics associated with the IBM method from the input
    std::string _solid_mechanics;

    //! The subdomain ids for the solid that is read from input
    std::set<libMesh::subdomain_id_type> _solid_subdomain_set;

    //! The subdomain ids for the fluid that are read from input
    std::set<libMesh::subdomain_id_type> _fluid_subdomain_set;

    //! The locator object used to find fluid elements
    libMesh::UniquePtr<libMesh::PointLocatorBase> _point_locator;

    //! Typedef to reduce code verbosity
    typedef std::map<libMesh::dof_id_type,std::vector<unsigned int> > SolidElemToQpMap;

    //! Typedef to reduce code verbosity
    /*! \todo Hash map for the first map? */
    typedef std::map<libMesh::dof_id_type,SolidElemToQpMap> FluidToSolidMap;

    /*! Map from fluid element id to solid element ids that touch the fluid element and the
        associated solid element quadrature points that lie in the fluid element. */
    FluidToSolidMap _fluid_id_to_solid_ids_qps;

    libMesh::UniquePtr<libMesh::FEMContext> _solid_context;

    //! Residual contributions to the fluid
    void element_time_derivative_fluid(AssemblyContext& context);

    //! Residual contributions to the solid
    void element_time_derivative_solid(AssemblyContext& context);

    void assemble_fluid_var_residual_contributions( bool compute_jacobian,
                                                    AssemblyContext & context );

    bool is_solid_elem( libMesh::subdomain_id_type elem_id );

    bool is_fluid_elem( libMesh::subdomain_id_type elem_id );

    ImmersedBoundary();
  };

  template<typename SolidMech>
  inline
  bool ImmersedBoundary<SolidMech>::is_solid_elem( libMesh::subdomain_id_type elem_id )
  {
    return _solid_subdomain_set.find(elem_id) != _solid_subdomain_set.end();
  }

  template<typename SolidMech>
  inline
  bool ImmersedBoundary<SolidMech>::is_fluid_elem( libMesh::subdomain_id_type elem_id )
  {
    return _fluid_subdomain_set.find(elem_id) != _fluid_subdomain_set.end();
  }

} //End namespace block
#endif //GRINS_IMMERSED_BOUNDARY_H
