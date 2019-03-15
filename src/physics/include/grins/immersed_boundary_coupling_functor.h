//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_IMMERSED_BOUNDARY_COUPLING_FUNCTOR_H
#define GRINS_IMMERSED_BOUNDARY_COUPLING_FUNCTOR_H

// GRINS
#include "grins/multi_component_vector_variable.h"
#include "grins/overlapping_fluid_solid_map.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/ghosting_functor.h"

namespace GRINS
{

  class ImmersedBoundaryCouplingFunctor : public libMesh::GhostingFunctor
  {
  public:

    ImmersedBoundaryCouplingFunctor( const MultiphysicsSystem & system,
                                     const libMesh::CouplingMatrix & coupling_matrix,
                                     const OverlappingFluidSolidMap & overlapping_map )
      : libMesh::GhostingFunctor(),
        _system(system),
        _overlapping_map(overlapping_map),
        _coupling_matrix(coupling_matrix)
    {}

    virtual ~ImmersedBoundaryCouplingFunctor(){};

    virtual void operator()
    ( const libMesh::MeshBase::const_element_iterator & range_begin,
      const libMesh::MeshBase::const_element_iterator & range_end,
      libMesh::processor_id_type p,
      std::unordered_map<const libMesh::Elem *,const libMesh::CouplingMatrix*> & coupled_elements ) override;

  private:

    const MultiphysicsSystem & _system;

    const OverlappingFluidSolidMap & _overlapping_map;

    const libMesh::CouplingMatrix & _coupling_matrix;

  };

} // end namespace GRINS

#endif // GRINS_IMMERSED_BOUNDARY_COUPLING_FUNCTOR_H
