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
    Static class to store dlvo coefficients

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef dlvoInfo_H
#define dlvoInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class interAdhesion Declaration
\*---------------------------------------------------------------------------*/


class dlvoInfo
{
    public:

    static const bool& isActive()
    {
        return active_;
    };
    static bool active_;

    static const scalar& getA()
    {
        return A_;
    };

    static const scalar& getEps0()
    {
        return eps0_;
    };

    static const scalar& getEpsR()
    {
        return epsR_;
    };

    static const scalar& getZeta()
    {
        return zeta_;
    };

    static const scalar& getRecK()
    {
        return recK_;
    };

    static const scalar& getCutOff()
    {
        return cutOff_;
    };

    static const scalar& getMinSurfDist()
    {
        return minSurfDist_;
    };

    static const scalar& getVdwMaxAcc()
    {
        return vdwMaxAcc_;
    };

    static const scalar& getEleMaxAcc()
    {
        return eleMaxAcc_;
    };

    static const scalar& getLubMaxAcc()
    {
        return lubMaxAcc_;
    };

    static scalar A_;
    static scalar eps0_;
    static scalar epsR_;
    static scalar zeta_;
    static scalar recK_;
    static scalar cutOff_;
    static scalar minSurfDist_;
    static scalar vdwMaxAcc_;
    static scalar eleMaxAcc_;
    static scalar lubMaxAcc_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
