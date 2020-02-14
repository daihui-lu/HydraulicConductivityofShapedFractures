/*-----------------------------------------------------------------------------/
Name:
    BJSlip

Description:

Author:
    Daihui Lu, Purdue University (2019)
/-----------------------------------------------------------------------------*/
#include "../parameters.H"


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

#};
code
#{

    const fvPatch& boundaryPatch = patch();
    const vectorField& Cf = boundaryPatch.Cf();
    vectorField& field = *this;
   // Info <<"\n this =  " << field  << endl;
    const scalar slipFactor_ = 0.02;
    const scalar alpha_ = -0.0005; // depth gradient of the cell
   // const scalar epsilon_ = $epsilon; //aspect ratio of the cell
    const scalar U_ave = 1; // velocity scale on the permeation direction

    vectorField nHat = this->patch().nf();
    vectorField gradient = this->snGrad(); //dudn
   // vectorField gradient = snGrad().field; //dudn
   // Info <<"\n gradient =  " << gradient  << endl;
  //  const scalar a = pow(0,0);
   // Info <<"\n a =  " << a  << endl;
    gradient = transform(I - sqr(nHat), gradient); // take the tangential component

    //vectorField gradientDirection = gradient / (mag(gradient) + SMALL);

    // slip velocity at the boundary = slip coefficient * velocity gradient
    vectorField ux = -slipFactor_* gradient;
    //   Info <<"\n ux =  " << ux << endl;
    //Info << ux + slipFactor_*gradient  << endl;

    scalarField uy(patch().size(),scalar(0));

    forAll(Cf,faceI)
    {
        // permeation velocity at the boundary
    	uy[faceI] =

	    U_ave
	*
	(
	    Cf[faceI].z()
	)
	/
	(
	    3.*slipFactor_+(Cf[faceI].z())
	)
	*
	    alpha_;


    }
    vectorField vw = uy*nHat;
    //Eq. 2-125 Leal's book

    operator==(ux - vw);
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
