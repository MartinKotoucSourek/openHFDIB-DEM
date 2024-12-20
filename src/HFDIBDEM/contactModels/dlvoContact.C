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

#include "dlvoInfo.H"
#include "periodicBody.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_ArbShape(
    dlvoContactInfo& cInfo,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    Info << "Not implemented yet" << endl;
    return {};
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_Sphere
(
    dlvoContactInfo& cInfo,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    scalar A = dlvoInfo::getA();
    scalar recDebay = dlvoInfo::getRecDebay();

    scalar cRadius = cInfo.getcClass().getGeomModel().getDC() / 2;
    scalar tRadius = cInfo.gettClass().getGeomModel().getDC() / 2;
    vector cCenter = cInfo.getcClass().getGeomModel().getCoM();
    vector tCenter = cInfo.gettClass().getGeomModel().getCoM();

    vector centerDir = tCenter - cCenter;
    scalar d = mag(centerDir);

    scalar surfDist = d - (cRadius + tRadius);
    surfDist = surfDist < dlvoInfo::getMinSurfDist() ? dlvoInfo::getMinSurfDist() : surfDist;
    Info << "surfDist: " << surfDist << endl;

    scalar F_VdW = A*(cRadius*tRadius/(cRadius + tRadius))/(6*surfDist*surfDist);
    scalar F_elec = 0;
    try {
        F_elec = -recDebay * (cRadius*tRadius/(cRadius + tRadius)) * dlvoInfo::getFactorZ() * exp(-recDebay * surfDist);
    }
    catch (...) {
        F_elec = 0;
    }

    scalar F_dlvo = F_VdW + F_elec;

    vector cDirNorm = centerDir/mag(centerDir);

    scalar lower_limit = 0.2 * dlvoInfo::getCharCellSize();
    scalar upper_limit = 1.5 * dlvoInfo::getCharCellSize();

    scalar limFunction = surfDist <= lower_limit ? 1 : surfDist >= upper_limit ? 0 : 0.5 * (1 + Foam::cos(3.14 * (surfDist - lower_limit) / (upper_limit - lower_limit)));

    vector cCntPointDir = cRadius * cDirNorm;
    vector tCntPointDir = - tRadius * cDirNorm;

    vector cPlanarVec =  cCntPointDir - cInfo.getcVars().Axis_*(cCntPointDir&cInfo.getcVars().Axis_);
    vector tPlanarVec =  tCntPointDir - cInfo.gettVars().Axis_*(tCntPointDir&cInfo.gettVars().Axis_);

    vector cCntPVel = (-(cPlanarVec^cInfo.getcVars().Axis_)*cInfo.getcVars().omega_ + cInfo.getcVars().Vel_);
    vector tCntPVel = (-(tPlanarVec^cInfo.gettVars().Axis_)*cInfo.gettVars().omega_ + cInfo.gettVars().Vel_);

    vector relVel = cCntPVel - tCntPVel;

    vector F_lubr = 6 * 3.14 * rhoF.value() * nuF.value() * pow((cRadius*tRadius/(cRadius + tRadius)), 2) * relVel / surfDist;
    vector F_lubr_relaxed = cInfo.getLastTangLubrForce() + dlvoInfo::getTanLubrRelax() * (F_lubr - cInfo.getLastTangLubrForce());
    cInfo.getLastTangLubrForce() = F_lubr_relaxed;

    vector T_c = limFunction * dlvoInfo::getTanLubrC() * (cCntPointDir ^ (-F_lubr_relaxed));
    vector T_t = limFunction * dlvoInfo::getTanLubrC() * (tCntPointDir ^ F_lubr_relaxed);

    Info << "F_dlvo: " << F_dlvo << " F_lubr_relaxed: " << F_lubr_relaxed << " F_lubr_relaxed & cDirNorm: " << (F_lubr_relaxed & cDirNorm) << endl;

    // vector F_c = F_dlvo * cDirNorm;
    vector F_c = (F_dlvo - (F_lubr_relaxed & cDirNorm)) * cDirNorm;
    vector F_t = - F_c;

    return {forces(F_c, T_c), forces(F_t, T_t)};
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact_Cluster
(
    dlvoContactInfo& cInfo,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
)
{
    Tuple2<forces,forces> returnF = {forces(vector::zero, vector::zero), forces(vector::zero, vector::zero)};

    std::vector<std::shared_ptr<geomModel>> cBodies;
    std::vector<std::shared_ptr<geomModel>> tBodies;

    bool isCCluster = cInfo.getcClass().getGeomModel().isCluster();

    if(isCCluster)
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(cInfo.getcClass().getGeomModel());
        std::vector<std::shared_ptr<geomModel>> tempBodies = cCluster.getClusterBodies();

        for(std::shared_ptr<geomModel>& cBody : tempBodies)
        {
            if (cBody->getM() == 0)
                continue;

            cBodies.push_back(cBody);
        }
    }
    else
    {
        cBodies.push_back(cInfo.getcClass().getGeomModelPtr());
    }

    isCCluster = cBodies.size() != 1;

    if(cInfo.gettClass().getGeomModel().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(cInfo.gettClass().getGeomModel());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.getClusterBodies()[0]);
        }
        else
        {
            std::vector<std::shared_ptr<geomModel>> tempBodies = tCluster.getClusterBodies();

            for(std::shared_ptr<geomModel>& tBody : tempBodies)
            {
                if (tBody->getM() == 0)
                    continue;

                tBodies.push_back(tBody);
            }
        }
    }
    else
    {
        tBodies.push_back(cInfo.gettClass().getGeomModelPtr());
    }

    scalar cMass = 0;
    scalar tMass = 0;

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        cMass += cgModel->getM();
    }

    for(std::shared_ptr<geomModel>& tgModel : tBodies)
    {
        tMass += tgModel->getM();
    }

    if (cMass == 0 || tMass == 0)
    {
        return {};
    }

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        for(std::shared_ptr<geomModel>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                cInfo.getcClass().getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                cInfo.gettClass().getMatInfo().getMaterial()
            );

            dlvoContactInfo tmpDlvoInfoI(
                cIbClassI,
                tIbClassI,
                cInfo.getcVars(),
                cInfo.gettVars()
            );

            Tuple2<forces,forces> tmpF = solveDlvoContact(tmpDlvoInfoI, nuF, rhoF);

            // mass average of forces
            try
            {
                tmpF.first().F *= (cgModel->getM()/cMass) * (tgModel->getM()/tMass);
                tmpF.second().F *= (cgModel->getM()/cMass) * (tgModel->getM()/tMass);
            }
            catch(const std::exception& e)
            {
                Info << "DLVO error: " << e.what() << endl;
            }

            returnF.first() += tmpF.first();
            returnF.second() += tmpF.second();
        }
    }

    return returnF;
}
//---------------------------------------------------------------------------//
Tuple2<forces,forces> solveDlvoContact
(
    dlvoContactInfo& cInfo,
    dimensionedScalar nuF,
    dimensionedScalar rhoF
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
            cInfo,
            nuF,
            rhoF
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
            cInfo,
            nuF,
            rhoF
        );
    }
    else
    {
        return solveDlvoContact_ArbShape(
            cInfo,
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
