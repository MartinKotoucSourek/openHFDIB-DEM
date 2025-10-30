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
    Static class for storing wall material information

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallMatInfo_H
#define wallMatInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallMatInfo Declaration
\*---------------------------------------------------------------------------*/

enum class wallDlvoModel
{
    NONE,
    CHARGED,
    SYMMETRIC
};

class wallMatInfo
{
    public:

    static const HashTable<materialInfo,string,Hash<string>>& getWallMatInfo()
    {
        return wallMatInfo_;
    }

    static const HashTable<wallDlvoModel,string,Hash<string>>& getWallDlvoInfo()
    {
        return wallDlvoInfo_;
    }

    static List<string> getWallPatches()
    {
        return wallMatInfo_.toc();
    }

    static void wallMatInfo_insert
    (
        const string& wallPatch,
        const materialInfo& matInfo,
        const wallDlvoModel& dlvoModel
    )
    {
        wallMatInfo_.insert(wallPatch,matInfo);
        wallDlvoInfo_.insert(wallPatch,dlvoModel);
    };

    //- Hash table for storing wall material information
    static HashTable<materialInfo,string,Hash<string>> wallMatInfo_;
    static HashTable<wallDlvoModel,string,Hash<string>> wallDlvoInfo_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
