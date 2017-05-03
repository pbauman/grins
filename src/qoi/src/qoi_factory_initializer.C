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
#include "grins/qoi_factory_initializer.h"

// GRINS
#include "grins/qoi_factory_basic.h"
#include "grins/qoi_names.h"

// GRINS-QoIs
#include "grins/average_nusselt_number.h"
#include "grins/vorticity.h"
#include "grins/parsed_boundary_qoi.h"
#include "grins/parsed_interior_qoi.h"
#include "grins/weighted_flux_qoi.h"

namespace GRINS
{
  QoIFactoryInitializer::QoIFactoryInitializer()
  {
    static QoIFactoryBasic<AverageNusseltNumber>
      grins_factory_avg_nusselt_qoi(avg_nusselt);

    static QoIFactoryBasic<ParsedBoundaryQoI>
      grins_factory_parsed_boundary_qoi(parsed_boundary);

    static QoIFactoryBasic<ParsedInteriorQoI>
      grins_factory_parsed_interior_qoi(parsed_interior);

    static QoIFactoryBasic<Vorticity>
      grins_factory_vorticity_qoi(vorticity);

    static QoIFactoryBasic<WeightedFluxQoI>
      grins_factory_weighted_flux_qoi(weighted_flux);
  }
} // end namespace GRINS
