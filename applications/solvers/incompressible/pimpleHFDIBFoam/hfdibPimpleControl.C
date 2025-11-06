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

\*---------------------------------------------------------------------------*/

#include "hfdibPimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hfdibPimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool Foam::hfdibPimpleControl::criteriaSatisfied()
{
    if (openHFDIBDEM_.minPimpleLoops() > 1 && corr_ < openHFDIBDEM_.minPimpleLoops())
    {
        if (debug)
        {
            Info<< algorithmName_ << " loop: minimal number of PIMPLE loops "
                << "not reached yet (" << corr_ << " < "
                << openHFDIBDEM_.minPimpleLoops() << ")" << endl;
        }
        return false;
    }

    return pimpleControl::criteriaSatisfied();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hfdibPimpleControl::hfdibPimpleControl
(
    fvMesh& mesh,
    openHFDIBDEM& openHFDIBDEM,
    const word& dictName,
    const bool verbose
)
:
    pimpleControl(mesh, dictName, verbose),
    openHFDIBDEM_(openHFDIBDEM)
{}

// ************************************************************************* //
