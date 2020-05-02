/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    This boundary condition provides a slip condition on the form of the
    first order Beaver Josepsh Saffmann interface between a clear fluid and
    a porous medium
        \f[
            U_p= - \phi (I - nn) \cdot \nabla_\perp U
        \f]

    where
    \vartable
        U_c   | internal field velocity
        U_p   | patch velocity
        \phi  | slip coefficient
        I     | identity matrix
        n     | face normal vector

    \endvartable

    The gradient in the normal direction is applied manually
Usage
    \table
        Property     | Description             | Required    | Default value
        slipCoeff    | slipCoefficient         | yes         |
        gradN        | gradient normal         | yes         |
    \endtable

SourceFiles
    porousInterfaceBJSGrad.C

Author
    Fedrico Municchi, Nottingham (2020)

References
    Bottaro, A. (2019). Flow over natural or engineered surfaces: an adjoint
    homogenization perspective. J. Fluid Mech, 877, 1.
    https://doi.org/10.1017/jfm.2019.607


\*---------------------------------------------------------------------------*/

#ifndef porousInterfaceBJSGrad_H
#define porousInterfaceBJSGrad_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class porousInterfaceBJSGrad
:
    public fixedValueFvPatchField<vector>
{
// Private data

    //- Slip coefficient
    scalarField     slipCoeff_;

    //- Permeability coefficient
    scalarField     gradN_;

public:

    //- Runtime type information
    TypeName("porousInterfaceBJSGrad");


    // Constructors

        //- Construct from patch and internal field
        porousInterfaceBJSGrad
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousInterfaceBJSGrad
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        porousInterfaceBJSGrad
        (
            const porousInterfaceBJSGrad&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        porousInterfaceBJSGrad
        (
            const porousInterfaceBJSGrad&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new porousInterfaceBJSGrad(*this)
            );
        }

        //- Copy constructor setting internal field reference
        porousInterfaceBJSGrad
        (
            const porousInterfaceBJSGrad&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField >
            (
                new porousInterfaceBJSGrad(*this, iF)
            );
        }


    //- Destructor
    virtual ~porousInterfaceBJSGrad();


    // Member Functions

        //- Do not fix value! Allow the flux to be adjusted to meet continuity
        virtual bool fixesValue() const
        {
            return false;
        }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Return slip coefficient
        virtual scalarField& slipCoeff()
        {
            return slipCoeff_;
        }

        virtual const scalarField& slipCoeff() const
        {
            return slipCoeff_;
        }

        //- Return permeability coefficient
        virtual scalarField& gradN()
        {
            return gradN_;
        }

        virtual const scalarField& gradN() const
        {
            return gradN_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
