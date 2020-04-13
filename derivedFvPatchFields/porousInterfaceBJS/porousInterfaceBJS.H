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
            U_p= - \phi (I - nn) \cdot \nabla_\perp U + nn \cdot U_c
        \f]

    where
    \vartable
        U_c   | internal field velocity
        U_p   | patch velocity
        \phi  | slip coefficient
        I     | identity matrix
        n     | face normal vector

    \endvartable

    Notice that it is possible to switch the condotion for the face-normal
    velocity from flully permeable to impermeable by simply adjusting a
    parameter \alpha.

    Specifically:
    \alpha = 1 corresponds to a fully permeable medium
    \alpha = 0 corresponds to a completely impermeable medium

    This boundary condition is derived from fixedValue and it is explicit.
    Thus, the condition on the pressure field will be fixedFluxPressure in most
    circumstances.

Usage
    \table
        Property     | Description             | Required    | Default value
        slipCoeff    | slipCoefficient         | yes         |
        alpha        | permeability coefficient| yes         |
    \endtable

SourceFiles
    porousInterfaceBJS.C

Author
    Fedrico Municchi, Nottingham (2020)

References
    Bottaro, A. (2019). Flow over natural or engineered surfaces: an adjoint
    homogenization perspective. J. Fluid Mech, 877, 1.
    https://doi.org/10.1017/jfm.2019.607


\*---------------------------------------------------------------------------*/

#ifndef porousInterfaceBJS_H
#define porousInterfaceBJS_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class porousInterfaceBJS
:
    public fixedValueFvPatchField<vector>
{
// Private data

    //- Slip coefficient
    scalarField     slipCoeff_;

    //- Permeability coefficient
    scalarField     alpha_;

public:

    //- Runtime type information
    TypeName("porousInterfaceBJS");


    // Constructors

        //- Construct from patch and internal field
        porousInterfaceBJS
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousInterfaceBJS
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        porousInterfaceBJS
        (
            const porousInterfaceBJS&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        porousInterfaceBJS
        (
            const porousInterfaceBJS&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new porousInterfaceBJS(*this)
            );
        }

        //- Copy constructor setting internal field reference
        porousInterfaceBJS
        (
            const porousInterfaceBJS&,
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
                new porousInterfaceBJS(*this, iF)
            );
        }


    //- Destructor
    virtual ~porousInterfaceBJS();


    // Member Functions

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
        virtual scalarField& alpha()
        {
            return alpha_;
        }

        virtual const scalarField& alpha() const
        {
            return alpha_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
