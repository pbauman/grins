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

#ifndef GRINS_IMMERSED_BOUNDARY_AUGMENTED_SPARSITY_H
#define GRINS_IMMERSED_BOUNDARY_AUGMENTED_SPARSITY_H

// GRINS
#include "grins/overlapping_fluid_solid_map.h"

// libMesh
#include "libmesh/dof_map.h"

namespace GRINS
{
  class MultiphysicsSystem;
  class DisplacementVariable;
  class VelocityVariable;
  class OverlappingFluidSolidMap;

  class ImmersedBoundaryAugmentedSparsity : public libMesh::DofMap::AugmentSparsityPattern
  {
  public:

    ImmersedBoundaryAugmentedSparsity( MultiphysicsSystem & system,
                                       const DisplacementVariable & disp_vars,
                                       const VelocityVariable & flow_vars,
                                       const OverlappingFluidSolidMap & fluid_solid_overlap )
      : _system(system),
        _disp_vars(disp_vars),
        _flow_vars(flow_vars),
        _fluid_solid_overlap(fluid_solid_overlap)
    {}

    ~ImmersedBoundaryAugmentedSparsity(){};

    //! Primary function of this class. Used by libMesh to augment sparsity pattern.
    /*! In the present case, we update n_nz (number of local dof "neighbors") and n_oz
        (number of off-processor dof "neighbors") for the fluid/solid overlap. We'll need
        both fluid and solid neighbors as we have contributions to both fluid and solid
        residuals from solid and fluid variables, respectively. */
    virtual void augment_sparsity_pattern( libMesh::SparsityPattern::Graph &,
                                           std::vector<libMesh::dof_id_type> & n_nz,
                                           std::vector<libMesh::dof_id_type> & n_oz );

  private:

    MultiphysicsSystem & _system;

    const DisplacementVariable & _disp_vars;

    const VelocityVariable & _flow_vars;

    //! Map from local fluid element id to solid element info
    /*! The FluidToSolidMap typedef is defined in the ImmersedBoundary Physics.
        It maps the fluid element id to a map of solid element ids and the quadrature
        point indices. Here, we only need the solid element ids, but we just
        grab the whole data structure instead of redundantly storing the separated
        info. */
    const OverlappingFluidSolidMap & _fluid_solid_overlap;

    //! Set sparsity for fluid/solid elem id "pairs"
    /*! Given the fluid element id and the set of overlapping solid element ids,
        this will augment the sparsity pattern for both solid element dof
        "neighbors" of each fluid dof and vice-versa. */
    void set_sparsity( MultiphysicsSystem & system,
                       const libMesh::Elem & fluid_elem,
                       const std::vector<libMesh::dof_id_type> & solid_elem_ids,
                       std::vector<libMesh::dof_id_type> & n_nz,
                       std::vector<libMesh::dof_id_type> & n_oz );

    bool dof_is_local( const libMesh::DofMap & dof_map, libMesh::dof_id_type dof ) const
    {
      return (dof >= dof_map.first_dof()) && (dof < dof_map.end_dof());
    }
  };

} // end namespace GRINS

#endif // GRINS_IMMERSED_BOUNDARY_AUGMENTED_SPARSITY_H
