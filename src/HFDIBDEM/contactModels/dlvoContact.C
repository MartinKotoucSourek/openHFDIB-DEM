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

    scalar vdwMaxForce = min(dlvoInfo::getVdwMaxAcc() * cInfo.getcClass().getGeomModel().getM(), dlvoInfo::getVdwMaxAcc() * cInfo.gettClass().getGeomModel().getM());
    scalar eleMaxForce = min(dlvoInfo::getEleMaxAcc() * cInfo.getcClass().getGeomModel().getM(), dlvoInfo::getEleMaxAcc() * cInfo.gettClass().getGeomModel().getM());
    // scalar lubMaxForce = min(dlvoInfo::getLubMaxAcc() * cInfo.getcClass().getGeomModel().getM(), dlvoInfo::getLubMaxAcc() * cInfo.gettClass().getGeomModel().getM());

    // Info << "sphere dlvo contact cGetM: " << cInfo.getcClass().getGeomModel().getM() << " tGetM: " << cInfo.gettClass().getGeomModel().getM() << endl;

    scalar cRadius = cInfo.getcClass().getGeomModel().getDC() / 2;
    scalar tRadius = cInfo.gettClass().getGeomModel().getDC() / 2;
    vector cCenter = cInfo.getcClass().getGeomModel().getCoM();
    vector tCenter = cInfo.gettClass().getGeomModel().getCoM();

    vector centerDir = tCenter - cCenter;
    scalar d = mag(centerDir);

    scalar surfDist = d - (cRadius + tRadius);
    surfDist = surfDist < dlvoInfo::getMinSurfDist() ? dlvoInfo::getMinSurfDist() : surfDist;

    scalar F_VdW = A*(cRadius*tRadius/(cRadius + tRadius))/(6*surfDist*surfDist);
    scalar F_elec = 0;
    if (surfDist/recK < 100)
    {
        F_elec = -4*3.14*eps_0*eps_r*zeta*zeta*(cRadius*tRadius/(cRadius + tRadius))/(recK*exp(surfDist/recK)+recK);
    }

    scalar limitForce = std::abs(F_VdW) > std::abs(F_elec) ? vdwMaxForce : eleMaxForce;

    // if (std::abs(F_VdW) > vdwMaxForce)
    // {
    //     Info << "VdW force exceeds the maximum allowed force. Maximum: " << vdwMaxForce << " Current: " << F_VdW << endl;
    //     F_VdW = sign(F_VdW) * vdwMaxForce;
    // }
    // else
    // {
    //     Info << "Max VdW force: " << vdwMaxForce << " Current: " << F_VdW << endl;
    // }
    // if (std::abs(F_elec) > eleMaxForce)
    // {
    //     Info << "Electrostatic force exceeds the maximum allowed force. Maximum: " << eleMaxForce << " Current: " << F_elec << endl;
    //     F_elec = sign(F_elec) * eleMaxForce;
    // }
    // else
    // {
    //     Info << "Max elec force: " << eleMaxForce << " Current: " << F_elec << endl;
    // }

    // Info << "MKS F_VdW: " << F_VdW << " F_elec: " << F_elec << " surfDist: " << surfDist << " limitForce: " << limitForce << endl;

    scalar F_dlvo = F_VdW + F_elec;

    if (std::abs(F_dlvo) > limitForce)
    {
        // Info << "DLVO force exceeds the maximum allowed force. Maximum: " << limitForce << " Current: " << F_dlvo << endl;
        F_dlvo = sign(F_dlvo) * limitForce;
    }

    vector cDirNorm = centerDir/mag(centerDir);

    // Lubrication force

    // scalar beta = tRadius / cRadius;
    // scalar beta_m1 = 1 / beta;
    // scalar psi_m1 = 1 / (surfDist / ((tRadius + cRadius) / 2));
    // scalar log_psi_m1 = log(psi_m1);

    // scalar Y_B_11 = -12.56 * pow(cRadius, 2) * ((beta*(4 + beta)*log_psi_m1 / (5*pow(1+beta, 2))));
    // scalar Y_B_21 = -12.56 * pow(tRadius, 2) * ((beta_m1*(4 + beta_m1)*log_psi_m1 / (5*pow(1+beta_m1, 2))));

    // scalar m8_pi_pow_cRadius_3 = 25.12 * pow(cRadius, 3);
    // scalar Y_C_11 = m8_pi_pow_cRadius_3 * (2.0*beta*log_psi_m1 / (5*(1+beta)));
    // scalar Y_C_12 = m8_pi_pow_cRadius_3 * (pow(beta, 2)*log_psi_m1 / (10*(1+beta)));

    // scalar m8_pi_pow_tRadius_3 = 25.12 * pow(tRadius, 3);
    // scalar Y_C_21 = m8_pi_pow_tRadius_3 * (pow(beta_m1, 2)*log_psi_m1 / (10*(1+beta_m1)));
    // scalar Y_C_22 = m8_pi_pow_tRadius_3 * (2.0*beta_m1*log_psi_m1 / (5*(1+beta_m1)));

    scalar lower_limit = 0.5 * dlvoInfo::getCharCellSize();
    scalar upper_limit = 2.5 * dlvoInfo::getCharCellSize();

    scalar limFunction = surfDist <= lower_limit ? 1 : surfDist >= upper_limit ? 0 : 0.5 * (1 + Foam::cos(3.14 * (surfDist - lower_limit) / (upper_limit - lower_limit)));

    vector cCntPointDir = cRadius * cDirNorm;
    vector tCntPointDir = - tRadius * cDirNorm;

    vector cPlanarVec =  cCntPointDir - cInfo.getcVars().Axis_*(cCntPointDir&cInfo.getcVars().Axis_);
    vector tPlanarVec =  tCntPointDir - cInfo.gettVars().Axis_*(tCntPointDir&cInfo.gettVars().Axis_);

    vector cCntPVel = (-(cPlanarVec^cInfo.getcVars().Axis_)*cInfo.getcVars().omega_ + cInfo.getcVars().Vel_);
    vector tCntPVel = (-(tPlanarVec^cInfo.gettVars().Axis_)*cInfo.gettVars().omega_ + cInfo.gettVars().Vel_);

    vector relativeTanVel = cCntPVel - tCntPVel;

    vector F_t_lubr = 6 * 3.14 * rhoF.value() * nuF.value() * pow((cRadius*tRadius/(cRadius + tRadius)), 2) * relativeTanVel / surfDist;

    // Info << "limFunction: " << limFunction << " cInfo.getcVars().Axis_: " << cInfo.getcVars().Axis_ << " cInfo.getcVars().omega_: " << cInfo.getcVars().omega_ << endl;
    // Info << "t: " << " cInfo.gettVars().Axis_: " << cInfo.gettVars().Axis_ << " cInfo.gettVars().omega_: " << cInfo.gettVars().omega_ << endl;

    Info << "cCenter: " << cCenter << endl;
    Info << "tCenter: " << tCenter << endl;
    Info << "cCntPointDir: " << cCntPointDir << endl;
    Info << "tCntPointDir: " << tCntPointDir << endl;
    Info << "F_t_lubr: " << F_t_lubr << endl;

    vector T_c = limFunction * dlvoInfo::getTanLubrC() * (cCntPointDir ^ (-F_t_lubr));
    vector T_t = limFunction * dlvoInfo::getTanLubrC() * (tCntPointDir ^ F_t_lubr);

    Info << "T_c: " << T_c << " T_t: " << T_t << endl;

    // Info << "T_C: " << T_c << " T_t: " << T_t << endl;
    vector F_c = F_dlvo * cDirNorm;
    // Info << "F_c: " << F_c << endl;
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
