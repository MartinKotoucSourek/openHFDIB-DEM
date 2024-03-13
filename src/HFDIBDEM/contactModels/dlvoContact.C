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
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "dlvoContact.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_ArbShape(
    dlvoContactInfo& cInfo
)
{
    Info << "Not implemented yet" << endl;
    return {};
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_Sphere
(
    dlvoContactInfo& cInfo
)
{
    scalar cRadius = cInfo.getcClass().getGeomModel().getDC() / 2;
    scalar tRadius = cInfo.gettClass().getGeomModel().getDC() / 2;
    vector tCenter = cInfo.gettClass().getGeomModel().getCoM();

    vector centerDir = cInfo.getcClass().getGeomModel().getCoM()
                        - tCenter;

    scalar d = mag(centerDir);

    scalar surfDist = d - (cRadius + tRadius);
    surfDist = surfDist < SMALL ? SMALL : surfDist;

    Info << "DLVO: Testing DLVO___________-" << endl;
    Info << "DLVO: surfDist: " << surfDist << endl;

    scalar A = 1e-20;
    scalar F_WdV = -A*(cRadius*tRadius/(cRadius + tRadius))/(6*surfDist*surfDist);

    scalar eps_0 = 8.854e-12;
    scalar eps_r = 50;
    scalar zeta = 2e-2;
    scalar rec_Debye = 5e-9;

    scalar F_elec = 4*3.14*eps_0*eps_r*zeta*zeta*(cRadius*tRadius/(cRadius + tRadius))/(rec_Debye*exp(surfDist/rec_Debye)+rec_Debye);

    Info << "DLVO: F_WdV: " << F_WdV << endl;
    Info << "DLVO: F_elec: " << F_elec << endl;
    scalar F_dlvo = F_WdV + F_elec;
    Info << "DLVO: F_dlvo: " << F_dlvo << endl;

    vector cDirNorm = centerDir/mag(centerDir);
    vector F_c = F_dlvo * cDirNorm;
    vector F_t = - F_c;
    Info << "DLVO: F_c: " << F_c << endl;
    Info << "DLVO: F_t: " << F_t << endl;

    return {forces(F_c, vector::zero), forces(F_t, vector::zero)};
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_Cluster
(
    dlvoContactInfo& cInfo
)
{
    Info << "Not implemented yet" << endl;
    return {};
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact
(
    dlvoContactInfo& cInfo
)
{
    if
    (
        cInfo.getcClass().getGeomModel().getcType() == sphere
        &&
        cInfo.gettClass().getGeomModel().getcType() == sphere
    )
    {
        return solveDlvoContact_Sphere(
            cInfo
        );
    }
    else if
    (
        cInfo.getcClass().getGeomModel().getcType() == cluster
        ||
        cInfo.gettClass().getGeomModel().getcType() == cluster
    )
    {
        return solveDlvoContact_Cluster(
            cInfo
        );
    }
    else
    {
        return solveDlvoContact_ArbShape(
            cInfo
        );
    }
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
