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
    scalar eps_0 = dlvoInfo::getEps0();
    scalar eps_r = dlvoInfo::getEpsR();
    scalar zeta = dlvoInfo::getZeta();
    scalar recK = dlvoInfo::getRecK();

    scalar cRadius = cInfo.getcClass().getGeomModel().getDC() / 2;
    scalar tRadius = cInfo.gettClass().getGeomModel().getDC() / 2;
    vector tCenter = cInfo.gettClass().getGeomModel().getCoM();

    vector centerDir = cInfo.getcClass().getGeomModel().getCoM()
                        - tCenter;

    scalar d = mag(centerDir);

    scalar surfDist = d - (cRadius + tRadius);
    surfDist = surfDist < dlvoInfo::getMinSurfDist() ? dlvoInfo::getMinSurfDist() : surfDist;

    scalar F_WdV = -A*(cRadius*tRadius/(cRadius + tRadius))/(6*surfDist*surfDist);
    scalar F_elec = 0;
    if (surfDist/recK < 100)
    {
        F_elec = 4*3.14*eps_0*eps_r*zeta*zeta*(cRadius*tRadius/(cRadius + tRadius))/(recK*exp(surfDist/recK)+recK);
    }

    scalar F_dlvo = F_WdV + F_elec;

    vector cDirNorm = centerDir/mag(centerDir);

    // Lubrication force
    scalar vn = -(cInfo.getcVars().Vel_ - cInfo.gettVars().Vel_) & cDirNorm;
    Info << "F_ vn: " << vn << endl;
    scalar F_lubr = 6*3.14*rhoF.value()*nuF.value()*pow(cRadius*tRadius/(cRadius + tRadius), 2)*vn/(surfDist);

    Info << "F_WdV: " << F_WdV << endl;
    Info << "F_elec: " << F_elec << endl;
    Info << "F_lubr: " << F_lubr << endl;

    vector F_c = (F_dlvo + F_lubr) * cDirNorm;
    Info << "F_c: " << F_c << endl;
    vector F_t = - F_c;

    return {forces(F_c, vector::zero), forces(F_t, vector::zero)};
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
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(cInfo.getcClass().getGeomModelPtr());
    }

    if(cInfo.gettClass().getGeomModel().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(cInfo.gettClass().getGeomModel());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.getClusterBodies()[0]);
        }
        else
        {
            tBodies = tCluster.getClusterBodies();
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
            tmpF.first().F *= (cgModel->getM()/cMass) * (tgModel->getM()/tMass);
            tmpF.second().F *= (cgModel->getM()/cMass) * (tgModel->getM()/tMass);

            returnF.first() += tmpF.first();
            returnF.second() += tmpF.second();
        }
    }

    Info << "Periodic F_c: " << returnF.first().F << endl;
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
