/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    contactModel

Description
    methods to solve particle-wall dlvo contact

SourceFiles
    dlvoContact.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef dlvoWallContact_H
#define dlvoWallContact_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"

#include "forces.H"
#include "ibContactVars.H"
#include "ibContactClass.H"

#include <optional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Namespace contactModel functions Declaration
\*---------------------------------------------------------------------------*/

namespace contactModel
{

    forces solveDlvoWallContact_ArbShape(
        ibContactClass& ibContact,
        ibContactVars& ibVars,
        List<std::shared_ptr<boundBox>> bboxes,
        dimensionedScalar nuF,
        dimensionedScalar rhoF
    );

    forces solveDlvoWallContact_Sphere
    (
        ibContactClass& ibContact,
        ibContactVars& ibVars,
        List<std::shared_ptr<boundBox>> bboxes,
        dimensionedScalar nuF,
        dimensionedScalar rhoF
    );

    forces solveDlvoWallContact_Cluster
    (
        ibContactClass& ibContact,
        ibContactVars& ibVars,
        List<std::shared_ptr<boundBox>> bboxes,
        dimensionedScalar nuF,
        dimensionedScalar rhoF
    );

    forces solveDlvoWallContact(
        ibContactClass& ibContact,
        ibContactVars& ibVars,
        List<std::shared_ptr<boundBox>> bboxes,
        dimensionedScalar nuF,
        dimensionedScalar rhoF
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
