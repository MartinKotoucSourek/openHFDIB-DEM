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
#include "dlvoWallContact.H"

#include "dlvoInfo.H"
#include "periodicBody.H"

// #include "parameters.H"
#include "wallMatInfo.H"
#include "wallPlaneInfo.H"
#include "dlvoContact.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
forces solveDlvoWallContact_ArbShape(
    ibContactClass& ibContact,
    ibContactVars& ibVars,
    List<std::shared_ptr<boundBox>> bboxes,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    Info << "Not implemented yet" << endl;
    return {};
}
//---------------------------------------------------------------------------//
forces solveDlvoWallContact_Sphere
(
    ibContactClass& ibContact,
    ibContactVars& ibVars,
    List<std::shared_ptr<boundBox>> bboxes,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    vector F = vector::zero;

    scalar A = dlvoInfo::getA();
    scalar recDebay = dlvoInfo::getRecDebay();

    scalar radius = ibContact.getGeomModel().getDC() / 2;

    std::vector<vector> centers;
    for (auto bBox : bboxes)
    {
        centers.push_back((bBox->min() + bBox->max()) / 2);
    }

    auto wallPatches = wallMatInfo::getWallPatches();
    for (auto const patchI : wallPatches)
    {
        switch (wallMatInfo::getWallDlvoInfo()[patchI])
        {
            case wallDlvoModel::NONE:
                break;
            case wallDlvoModel::CHARGED:
            {
                auto planeInfo = wallPlaneInfo::getWallPlaneInfo()[patchI];
                scalar nn = magSqr(planeInfo[0]);

                if (nn < VSMALL)
                {
                    Info << "Wall patch " << patchI << " has zero normal vector. Skipping DLVO contact." << endl;
                    break;
                }

                scalar surfDist = VGREAT;
                for (auto cC : centers)
                {
                    scalar t = ((cC - planeInfo[1]) & planeInfo[0]) / nn;
                    vector pPoint = cC - t * planeInfo[0];
                    surfDist = min(mag(cC - pPoint), surfDist);
                }

                surfDist = surfDist - radius;
                surfDist = surfDist < dlvoInfo::getMinSurfDist() ? dlvoInfo::getMinSurfDist() : surfDist;

                scalar F_VdW = A * radius / (6 * surfDist * surfDist);
                scalar F_elec = 0;
                try {
                    F_elec = -recDebay * radius * dlvoInfo::getFactorZ() * exp(-recDebay * surfDist);
                }
                catch (...) {
                    F_elec = 0;
                }

                scalar F_dlvo = dlvoInfo::useDLVO() ? (F_VdW + F_elec) : 0;

                vector cDirNorm = planeInfo[0]/mag(planeInfo[0]);

                F += F_dlvo * cDirNorm;

                break;
            }
            case wallDlvoModel::SYMMETRIC:
            {
                // std::shared_ptr<geomModel> iBcopy(ibContact.getGeomModel().getCopy());
                // iBcopy->bodyMovePoints(vector::zero);

                auto planeInfo = wallPlaneInfo::getWallPlaneInfo()[patchI];
                scalar nn = magSqr(planeInfo[0]);

                if (nn < VSMALL)
                {
                    Info << "Wall patch " << patchI << " has zero normal vector. Skipping DLVO contact." << endl;
                    break;
                }

                vector closestCenter = centers[0];
                vector planePoint = vector::zero;
                scalar surfDist = VGREAT;
                for (auto cC : centers)
                {
                    scalar t = ((cC - planeInfo[1]) & planeInfo[0]) / nn;
                    vector pPoint = cC - t * planeInfo[0];
                    scalar tmpDist = mag(cC - pPoint);
                    if (tmpDist < surfDist)
                    {
                        surfDist = tmpDist;
                        closestCenter = cC;
                        planePoint = pPoint;
                    }
                }

                std::shared_ptr<geomModel> iBcopy(ibContact.getGeomModel().getCopy());
                iBcopy->bodyMovePoints(2 * (planePoint - closestCenter));

                ibContactClass copyIbClass(
                    iBcopy,
                    ibContact.getMatInfo().getMaterial()
                );

                ibContactVars copyIbVars = ibVars;
                // symetry the velocity
                copyIbVars.Vel_ = ibVars.Vel_ - 2 * ((ibVars.Vel_ & planeInfo[0]) / nn) * planeInfo[0];
                copyIbVars.Axis_ = ibVars.Axis_ - 2 * ((ibVars.Axis_ & planeInfo[0]) / nn) * planeInfo[0];

                dlvoContactInfo tmpInfo(
                    ibContact,
                    copyIbClass,
                    ibVars,
                    copyIbVars,
                    bboxes,
                    iBcopy->getBBoxes()
                );

                Tuple2<forces,forces> tmpF = solveDlvoContact(tmpInfo, nuF, rhoF);

                F += tmpF.first().F;
            }
        }
    }

    return forces(F, vector::zero);
}
//---------------------------------------------------------------------------//
forces solveDlvoWallContact_Cluster
(
    ibContactClass& ibContact,
    ibContactVars& ibVars,
    List<std::shared_ptr<boundBox>> bboxes,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    forces returnF = forces(vector::zero, vector::zero);

    std::vector<std::shared_ptr<geomModel>> bodies;
    periodicBody& cluster = dynamic_cast<periodicBody&>(ibContact.getGeomModel());
    std::vector<std::shared_ptr<geomModel>> tempBodies = cluster.getClusterBodies();
    for(std::shared_ptr<geomModel>& cBody : tempBodies)
    {
        if (cBody->getM() == 0)
            continue;

        bodies.push_back(cBody);
    }

    scalar mass = 0;

    for(std::shared_ptr<geomModel>& model : bodies)
    {
        mass += model->getM();
    }

    if (mass == 0)
    {
        return {};
    }

    for(std::shared_ptr<geomModel>& model : bodies)
    {
        ibContactClass cIbClassI(
            model,
            ibContact.getMatInfo().getMaterial()
        );

        forces tmpF = solveDlvoWallContact(cIbClassI, ibVars, bboxes, nuF, rhoF);
        return tmpF;
    }

    return returnF;
}
//---------------------------------------------------------------------------//
forces solveDlvoWallContact
(
    ibContactClass& ibContact,
    ibContactVars& ibVars,
    List<std::shared_ptr<boundBox>> bboxes,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    if
    (
        ibContact.getGeomModel().getcType() == sphere
    )
    {
        return solveDlvoWallContact_Sphere(
            ibContact,
            ibVars,
            bboxes,
            nuF,
            rhoF
        );
    }
    else if
    (
        ibContact.getGeomModel().getcType() == cluster
    )
    {
        return solveDlvoWallContact_Cluster(
            ibContact,
            ibVars,
            bboxes,
            nuF,
            rhoF
        );
    }
    else
    {
        return solveDlvoWallContact_ArbShape(
            ibContact,
            ibVars,
            bboxes,
            nuF,
            rhoF
        );
    }
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
