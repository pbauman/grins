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

    //! Mixture molar mass (molecular weight), [kg/mol]
    libMesh::Real M_mix( const std::vector<libMesh::Real> & mass_fractions ) const;

    //! Species gas constant, [J/kg-K]
    /*! R_universal/M(species) */
    libMesh::Real R( unsigned int species ) const
    { return static_cast<DerivedType*>(this)->R_imp(species); }

    //! Mixture gas constant, [J/kg-K]
    libMesh::Real R_mix( const std::vector<libMesh::Real> & mass_fractions ) const;

    //! Species mole fraction, unitless
    libMesh::Real X( unsigned int species, libMesh::Real M_mix, libMesh::Real mass_fraction ) const
    { return mass_fraction*M_mix/this->M(species); }

    //! Mole fraction for all species, unitless
    void X( libMesh::Real M_mix,
            const std::vector<libMesh::Real> & mass_fractions,
            std::vector<libMesh::Real> & mole_fractions ) const;

    unsigned int n_species() const
    { return static_cast<DerivedType*>(this)->n_species_imp(); }

    unsigned int species_index( const std::string & species_name ) const
    { return static_cast<DerivedType*>(this)->species_index_imp(species_name); }

    std::string species_name( unsigned int species_index ) const
    { return static_cast<DerivedType*>(this)->species_name_imp(species_index); }

    const DerivedType & mixture() const
    { return *(static_cast<DerivedType*>(this)); }

  };

  template <typename DerivedType>
  inline
  libMesh::Real ThermochemicalMixtureBase<DerivedType>::
  M_mix( const std::vector<libMesh::Real> & mass_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), this->n_species() );

    libMesh::Real M = 0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      M += mass_fractions[s]/(this->M(s));

    return 1.0/M;
  }

  template <typename DerivedType>
  inline
  libMesh::Real ThermochemicalMixtureBase<DerivedType>::
  R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), this->n_species() );

    libMesh::Real R = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      R += mass_fractions[s]*this->R(s);

    return R;
  }

  template <typename DerivedType>
  inline
  void ThermochemicalMixtureBase<DerivedType>::X( libMesh::Real M_mix,
                                                  const std::vector<libMesh::Real> & mass_fractions,
                                                  std::vector<libMesh::Real> & mole_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), this->n_species() );
    libmesh_assert_equal_to( mole_fractions.size(), mass_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      mole_fractions[s] = this->X(s, M_mix, mass_fractions[s]);
  }

} // end namespace GRINS

#endif // GRINS_THERMOCHEMICAL_MIXTURE_BASH_H
