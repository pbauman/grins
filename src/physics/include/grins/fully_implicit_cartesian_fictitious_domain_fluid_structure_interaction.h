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

#ifndef GRINS_FULLY_IMPLICIT_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H
#define GRINS_FULLY_IMPLICIT_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H

// GRINS
#include "grins/cartesian_fictitious_domain_fluid_structure_interaction_base.h"

namespace GRINS
{
  template<unsigned int Dim>
  class FullyImplicitCartesianFictitiousDomainFluidStructureInteraction :
    public CartesianFictitiousDomainFluidStructureInteractionBase<Dim,false>
  {
  public:

    FullyImplicitCartesianFictitiousDomainFluidStructureInteraction( const PhysicsName & physics_name,
                                                                     const GetPot & input )
      : CartesianFictitiousDomainFluidStructureInteractionBase<Dim,false>(physics_name,input)
    {}

    FullyImplicitCartesianFictitiousDomainFluidStructureInteraction() = delete;

    virtual ~FullyImplicitCartesianFictitiousDomainFluidStructureInteraction() = default;

    //! Cache mesh information needed for residual computation
    virtual void preassembly( MultiphysicsSystem & system ) override
    { this->reinit_overlapping_data(system,false); }

  };

} // end namespace GRINS
#endif // GRINS_FULLY_IMPLICIT_CARTESIAN_FICTITIOUS_DOMAIN_FLUID_STRUCTURE_INTERACTION_H
