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
#include "ibContactClass.H"

using namespace Foam;

//---------------------------------------------------------------------------//
ibContactClass::ibContactClass
(
    geomModel& geomModel,
    materialInfo& matInfo
)
:
geomModel_(geomModel),
isInWallContact_(false),
inContactWithStatic_(false),
timeStepsInContWStatic_(0),
matInfo_(matInfo)
{
    surfCells_.setSize(Pstream::nProcs());
    intCells_.setSize(Pstream::nProcs());
}
ibContactClass::~ibContactClass()
{
}
//---------------------------------------------------------------------------//
bool ibContactClass::cellNotInLists(label cell)
{
    forAll(surfCells_[Pstream::myProcNo()],cellI)
    {
        if(surfCells_[Pstream::myProcNo()][cellI] == cell)
        {
            return false;
        }
    }

    forAll(intCells_[Pstream::myProcNo()],cellI)
    {
        if(intCells_[Pstream::myProcNo()][cellI] == cell)
        {
            return false;
        }
    }

    return true;
}
//---------------------------------------------------------------------------//
void ibContactClass::appendLists
(
    const List<DynamicLabelList>& surfCells,
    const List<DynamicLabelList>& intCells
)
{
    forAll(surfCells[Pstream::myProcNo()],cellI)
    {
        if(cellNotInLists(surfCells[Pstream::myProcNo()][cellI]))
        {
            surfCells_[Pstream::myProcNo()].append(surfCells[Pstream::myProcNo()][cellI]);
        }
    }

    forAll(intCells[Pstream::myProcNo()],cellI)
    {
        if(cellNotInLists(intCells[Pstream::myProcNo()][cellI]))
        {
            intCells_[Pstream::myProcNo()].append(intCells[Pstream::myProcNo()][cellI]);
        }
    }
}
//---------------------------------------------------------------------------//
