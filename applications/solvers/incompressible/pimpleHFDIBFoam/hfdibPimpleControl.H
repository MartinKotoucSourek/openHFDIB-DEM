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

#include "pimpleNoLoopControl.H"
#include "hfdibPimpleLoop.H"
#include "singleRegionConvergenceControl.H"
#include "singleRegionCorrectorConvergenceControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hfdibPimpleControl Declaration
\*---------------------------------------------------------------------------*/

class hfdibPimpleControl
:
    public pimpleNoLoopControl,
    public hfdibPimpleLoop
{
public:

    // Static Data Members

        //- Run-time type information
        TypeName("hfdibPimpleControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        hfdibPimpleControl(fvMesh& mesh, openHFDIBDEM& openHFDIBDEM, const word& algorithmName="PIMPLE");


    //- Destructor
    virtual ~hfdibPimpleControl();


    // Member Functions

        // IO

            //- Read controls
            virtual bool read();

        // Access

            //- Flag to indicate whether to solve the turbulence
            inline bool turbCorr() const;

        // Evolution

            //- Pimple loop
            bool loop();

            //- Time run loop
            bool run(Time& time);

            //- Time loop loop
            bool loop(Time& time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::hfdibPimpleControl::turbCorr() const
{
    return !turbOnFinalIterOnly() || finalPimpleIter();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
