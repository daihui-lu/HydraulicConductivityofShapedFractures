/*-----------------------------------------------------------------------------/
Name:
    BJSlip

Description:

Author:
    Daihui Lu, Purdue University (2019)
/-----------------------------------------------------------------------------*/

type            codedFixedValue;
value           uniform (0 0 0);
name            BJSlip;
codeInclude
#{
    #include "fvCFD.H"
    #include <cmath>
    #include <iostream>
    #include "fieldTypes.H"
    #include "volFields.H"
    #include "surfaceFields.H"
    #include "symmTransformField.H"
    #include "fieldTypes.H"
//#include "singlePhaseTransportModel.H"
#};
code
#{

    const fvPatch& boundaryPatch = patch();
    const vectorField& Cf = boundaryPatch.Cf();
    vectorField& field = *this;

    const scalar slipFactor_ = 0.00001;
    const scalar alpha_ = 0; // depth gradient of the cell
    const scalar v_0 = 0; // velocity scale on the permeation direction

    vectorField nHat = this->patch().nf();
    vectorField gradient = this->snGrad();
    gradient = transform(I - sqr(nHat), gradient);

    const scalar epsilon = 0.00001; //aspect ratio of the cell
    	vectorField gradientDirection = gradient / (mag(gradient) + SMALL);

    // slip velocity at the boundary = slip coefficient * velocity gradient
    	vectorField ux = -slipFactor_* mag(gradient)*gradientDirection;
    scalarField uy(patch().size(),scalar(0));

    forAll(Cf,faceI)
    {
        // permeation velocity at the boundary
    	uy[faceI] =
        (
            v_0*scalar(4.)
            *pow( scalar(1.)+ (alpha_/epsilon)*Cf[faceI].x(), 2)
        )
        *
        (
            alpha_/
            (
                epsilon
                *
                (
                    scalar(6.)*slipFactor_
                    +
                    (
                        scalar(2.)
                        *
                        pow(scalar(1.)+(alpha_/epsilon)*Cf[faceI].x(), 2)
                    )
                )
            )
        );

    }
    vectorField vw = uy*nHat;
    //Eq. 2-125 Leal's book

    operator==(ux - vw);// plus or minus?? nHat goes outside or inside?
#};

codeOptions
#{
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude
#};

codeLibs
#{
    -lmeshTools \
	-lfiniteVolume
#};
