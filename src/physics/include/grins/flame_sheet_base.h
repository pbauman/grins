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

//GRINS
#include "grins/physics.h"
#include "grins/primitive_flow_fe_variables.h"
#include "grins/primitive_temp_fe_variables.h"

#ifndef GRINS_FLAME_SHEET_BASE_H
#define GRINS_FLAME_SHEET_BASE_H

namespace GRINS
{
  template<typename Mixture, typename Evaluator>
  class FlameSheetBase : public Physics
  {
  public:

    FlameSheetBase( const PhysicsName& physics_name, const GetPot& input );
    virtual ~FlameSheetBase(){};

    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( AssemblyContext& context );

    const Mixture& gas_mixture() const;

  protected:

    PrimitiveFlowFEVariables _flow_vars;

    PrimitiveTempFEVariables _temp_vars;

    //! FE family for flame sheet variable
    GRINSEnums::FEFamily _S_FE_family;

    //! FE order for flame sheet variable
    GRINSEnums::Order _S_FE_order;

    //! index for flame sheet variable
    VariableIndex _S_var;

    Mixture _gas_mixture;

  private:

    FlameSheetBase();

  };

  template<typename Mixture, typename Evaluator>
  inline
  const Mixture& FlameSheetBase<Mixture,Evaluator>::gas_mixture() const
  {
    return _gas_mixture;
  }

} // end namespace GRINS

#endif // GRINS_FLAME_SHEET_BASE_H
