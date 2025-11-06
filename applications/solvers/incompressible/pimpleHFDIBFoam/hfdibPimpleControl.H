/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

Class
    Foam::hfdibPimpleControl

Description
    Pimple control class. Provides time-loop control methods which exit the
    simulation once convergence criteria have been reached. Also provides
    Pimple-loop control methods which exit the iteration once corrector
    convergence criteria have been met. Example usage:

    \verbatim
    hfdibPimpleControl pimple(mesh);

    while (pimple.run(runTime))
    {
        // pre-time-increment operations ...

        runTime ++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            // solve ...
        }

        // post-solve operations ...
    }
    \endverbatim

SourceFiles
    hfdibPimpleControlI.H
    hfdibPimpleControl.C

\*---------------------------------------------------------------------------*/

#ifndef hfdibPimpleControl_H
#define hfdibPimpleControl_H

#include "pimpleControl.H"
#include "openHFDIBDEM.H"

//- Declare that pimpleControl will be used
#define HFDIB_PIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class hfdibPimpleControl Declaration
\*---------------------------------------------------------------------------*/

class hfdibPimpleControl
:
    public pimpleControl
{
    // Private member functions

        //- No copy construct
        hfdibPimpleControl(const hfdibPimpleControl&) = delete;

        //- No copy assignment
        void operator=(const hfdibPimpleControl&) = delete;


protected:

    // Protected data

        // Solution controls

            //- Reference to openHFDIBDEM
            openHFDIBDEM& openHFDIBDEM_;

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("hfdibPimpleControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        hfdibPimpleControl
        (
            fvMesh& mesh,
            openHFDIBDEM& openHFDIBDEM,
            const word& dictName = "PIMPLE",
            const bool verbose = true
        );


    //- Destructor
    virtual ~hfdibPimpleControl() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
