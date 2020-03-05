/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "porousInterfaceBJS.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvCFD.H"
#include "symmTransformField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    porousInterfaceBJS
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porousInterfaceBJS::porousInterfaceBJS
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    slipCoeff_(p.size())
{
}


porousInterfaceBJS::porousInterfaceBJS
(
    const porousInterfaceBJS& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    slipCoeff_(mapper(ptf.slipCoeff_))
{
}


porousInterfaceBJS::porousInterfaceBJS
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    slipCoeff_("slipCoeff", dict, p.size())
{
}


porousInterfaceBJS::porousInterfaceBJS
(
    const porousInterfaceBJS& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    slipCoeff_(ptf.slipCoeff_)
{
}


porousInterfaceBJS::porousInterfaceBJS
(
    const porousInterfaceBJS& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    slipCoeff_(ptf.slipCoeff_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

porousInterfaceBJS::~porousInterfaceBJS()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void porousInterfaceBJS::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const vectorField nHat(this->patch().nf());
    const vectorField pif(this->patchInternalField());  

    //- Coefficients are updated based on a discretized gradient
    //
    //  slipVel = -slipCoeff * (I-nn) & (U_p - U_c) * Delta
    //
    //  Since we already know that (I-nn)&U_p is slipVel, we can write:
    //
    //  slipVel*(1.0 + slipCoeff*Delta) = slipCoeff*Delta*(I-nn)&U_c
    // 
    //  Resulting in the expression below.

    const vectorField slipVelocity
    (
        slipCoeff()
      * transform(
            I - sqr(nHat),
            pif * this->patch().deltaCoeffs()
        )
      / (
          scalar(1.0) + slipCoeff() * this->patch().deltaCoeffs()
        )
    );
    
    //- Here the normal velocity is simply the normal component of the
    //  internal field (i.e., zero normal gradient)

    const vectorField normalVelocity(pif&sqr(nHat));

    //- Finally, the boundary velocity is computed
    operator==(normalVelocity + slipVelocity);

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

