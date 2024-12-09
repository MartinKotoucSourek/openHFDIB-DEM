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

Description
    class for storage of DLVO contact variables

SourceFiles
    dlvoContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef dlvoContactInfo_H
#define dlvoContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "prtSubContactInfo.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class prtContactInfo Declaration
\*---------------------------------------------------------------------------*/

class dlvoContactInfo
{
    private:

    ibContactClass& cIbContactClass_;
    ibContactClass& tIbContactClass_;

    ibContactVars cIbContactVars_;
    ibContactVars tIbContactVars_;

    Tuple2<label,label> contactPair_;

    vector lastTangLubrForce_;

    public:

    // Constructors ------------------------------------------------
    dlvoContactInfo
    (
        ibContactClass& cClass,
        ibContactClass& tClass,
        ibContactVars& cVars,
        ibContactVars& tVars
    );

    // Destructors -------------------------------------------------
    ~dlvoContactInfo();
    //-----------------------------Func-------------------------------//

    const Tuple2<label, label>& getCPair()
    {
        return contactPair_;
    }

    ibContactClass& getcClass()
    {
        return cIbContactClass_;
    }

    ibContactClass& gettClass()
    {
        return tIbContactClass_;
    }

    ibContactVars& getcVars()
    {
        return cIbContactVars_;
    }

    ibContactVars& gettVars()
    {
        return tIbContactVars_;
    }

    vector& getLastTangLubrForce()
    {
        return lastTangLubrForce_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
