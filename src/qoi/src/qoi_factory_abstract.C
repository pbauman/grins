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
#include "grins/qoi_factory_abstract.h"

// GRINS
#include "grins/string_utils.h"

namespace GRINS
{
  libMesh::UniquePtr<CompositeQoI> QoIFactoryAbstract::create()
  {
    if( !this->_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building QoIs!");

    // Extract out the list of QoIs requested by the User
    std::string qoi_list = (*(this->_input))("QoI/enabled_qois", "none" );

    std::vector<std::string> qoi_names;

    if( qoi_list != std::string("none") )
      StringUtilities::split_string( qoi_list, std::string(" "), qoi_names );

    libMesh::UniquePtr<CompositeQoI> comp_qoi;

    for( std::vector<std::string>::const_iterator name = qoi_names.begin();
         name != qoi_names.end(); ++name )
      {
        libMesh::UniquePtr<QoIBase> qoi_ptr( this->build_qoi( *(this->_input) ) );
        comp_qoi->add_qoi(qoi_ptr);
      }

    // Report on which QoIs were built
    if( (*(this->_input))( "screen-options/echo_qoi", false ) )
      this->echo_qoi_list( comp_qoi );

    return comp_qoi;
  }

  void QoIFactoryAbstract::echo_qoi_list( libMesh::UniquePtr<CompositeQoI> & comp_qoi )
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    libMesh::out << "==========================================================" << std::endl
                 << "List of Enabled QoIs:" << std::endl;

    for( unsigned int q = 0; q < comp_qoi->n_qois(); q++ )
      libMesh::out << comp_qoi->get_qoi(q).name() << std::endl;

     libMesh::out <<  "==========================================================" << std::endl;
  }

  // Full specialization for the Factory<QoI>
  template<>
  std::map<std::string, FactoryAbstract<CompositeQoI>*>&
  FactoryAbstract<CompositeQoI>::factory_map()
  {
    static std::map<std::string, FactoryAbstract<CompositeQoI>*> _map;
    return _map;
  }

  // Definition of static members
  template<>
  const GetPot* FactoryWithGetPot<CompositeQoI>::_input = NULL;

} // end namespace GRINS
