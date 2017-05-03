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

#ifndef GRINS_QOI_FACTORY_ABSTRACT_H
#define GRINS_QOI_FACTORY_ABSTRACT_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/composite_qoi.h"

namespace GRINS
{
  class QoIFactoryAbstract : public FactoryWithGetPot<CompositeQoI>
  {
  public:
    QoIFactoryAbstract( const std::string & qoi_name )
      : FactoryWithGetPot<CompositeQoI>(qoi_name)
    {}

    virtual ~QoIFactoryAbstract() =0;

  protected:

    //! Subclasses construct each individual QoI using this method
    virtual libMesh::UniquePtr<QoIBase> build_qoi( const GetPot & input ) =0;

  private:

    virtual libMesh::UniquePtr<CompositeQoI> create();

    void echo_qoi_list( libMesh::UniquePtr<CompositeQoI> & comp_qoi );

    QoIFactoryAbstract();
  };

  inline
  QoIFactoryAbstract:: ~QoIFactoryAbstract(){}

} // end namespace GRINS

#endif // GRINS_QOI_FACTORY_ABSTRACT_H
