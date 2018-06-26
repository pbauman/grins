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

#ifndef GRINS_THERMOCHEMICAL_MIXTURE_BASH_H
#define GRINS_THERMOCHEMICAL_MIXTURE_BASH_H

namespace GRINS
{
  //! Interface for thermochemical mixture objects
  /*! This class defines the interface for thermochemical mixtures.
      This includes functions related to chemistry, thermodynamics, and
      chemical kinetics of the mixture. DerivedType objects will worry
      about details of parsing, setup, etc. Note that this interface
      assumes quantities are returned in SI units: kg, mol, s, K, m.
      We use a CRTP pattern for static polymorphism.

      Note that this class is intended to be constructed before the forking
      of threads. Use the Evaluator classes following the forking of threads
      to ensure thread safety.
   */
  template <typename DerivedType>
  class ThermochemicalMixtureBase
  {
  public:

    ThermochemicalMixtureBase(){};

    virtual ~ThermochemicalMixtureBase(){};

    //! Species molar mass (molecular weight), [kg/mol]
    libMesh::Real M( unsigned int species ) const
    { return static_cast<DerivedType*>(this)->M_imp(species); }

    //! Species gas constant, [J/kg-K]
    /*! R_universal/M(species) */
    libMesh::Real R( unsigned int species ) const
    { return static_cast<DerivedType*>(this)->R_imp(species); }

    unsigned int n_species() const
    { return static_cast<DerivedType*>(this)->n_species_imp(); }

    unsigned int species_index( const std::string & species_name ) const
    { return static_cast<DerivedType*>(this)->species_index_imp(species_name); }

    std::string species_name( unsigned int species_index ) const
    { return static_cast<DerivedType*>(this)->species_name_imp(species_index); }

    const DerivedType & mixture() const
    { return *(static_cast<DerivedType*>(this)); }

  };

} // end namespace GRINS

#endif // GRINS_THERMOCHEMICAL_MIXTURE_BASH_H
