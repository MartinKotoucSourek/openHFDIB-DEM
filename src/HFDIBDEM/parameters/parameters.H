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
    DEM parameters

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef parameters_H
#define parameters_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialProperties.H"
#include "interAdhesion.H"
#include "wallMatInfo.H"
#include "wallPlaneInfo.H"
#include "virtualMeshLevel.H"
#include "cyclicPlaneInfo.H"
#include "contactModelInfo.H"
#include "solverInfo.H"
#include "dlvoInfo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Parameters Definition
\*---------------------------------------------------------------------------*/

HashTable<materialInfo,string,Hash<string>> materialProperties::materialProps_ =
    HashTable<materialInfo,string,Hash<string>>();

HashTable<scalar,string,Hash<string>> interAdhesion::matInterAdh_ =
    HashTable<scalar,string,Hash<string>>();

HashTable<materialInfo,string,Hash<string>> wallMatInfo::wallMatInfo_ =
    HashTable<materialInfo,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> wallPlaneInfo::wallPlaneInfo_ =
    HashTable<List<vector>,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> cyclicPlaneInfo::cyclicPlaneInfo_ =
    HashTable<List<vector>,string,Hash<string>>();

std::map<string, string> cyclicPlaneInfo::cyclicPatchPairs_ =
    std::map<string, string>();

scalar virtualMeshLevel::virtualMeshLevel_ = 1;
scalar virtualMeshLevel::levelOfDivision_ = 1;
scalar virtualMeshLevel::charCellSize_ = 1;
scalar contactModelInfo::LcCoeff_ = 4.0;
label  solverInfo::nSolidsTreshnold_ = 1000;
bool solverInfo::onlyDEM_ = false;
bool contactModelInfo::useChenRotationalModel_ = false;
bool contactModelInfo::useMindlinRotationalModel_ = false;
bool contactModelInfo::increasedDamping_ = false;

bool   dlvoInfo::active_ = false;
bool   dlvoInfo::useDLVO_ = true;
bool   dlvoInfo::useTangLubr_ = true;
bool   dlvoInfo::useTransLubr_ = true;
scalar dlvoInfo::A_ = 1e-20;
scalar dlvoInfo::eps0_ = 8.854e-12;
scalar dlvoInfo::epsR_ = 50;
scalar dlvoInfo::zeta_ = 2e-2;
scalar dlvoInfo::recDebay_ = 1.0 / 2.5e-8;
scalar dlvoInfo::factorZ = 0;
scalar dlvoInfo::cutOff_ = 1e-5;
scalar dlvoInfo::minSurfDist_ = 1e-10;
scalar dlvoInfo::tanLubrC_ = 0.1;
scalar dlvoInfo::tanLubrRelax_ = 1;
scalar dlvoInfo::charCellSize_ = 0.1;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //