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
    Static class for Brownina force

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef brownianForce_H
#define brownianForce_H

#include <math.h>
#include "ibContactClass.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class interAdhesion Declaration
\*---------------------------------------------------------------------------*/


class brownianForce
{
    public:

    vector computeBrownianForce(ibContactClass& ibClass, scalar dt)
    {
        const scalar dp = ibClass.getGeomModel().getDC();           // diameter [m]
        const scalar a  = 0.5*dp;                                   // radius  [m]
        // const scalar Tc = td.Tc();

        // const scalar alpha = 2.0*lambda_/dp;
        // const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

        // // Boltzmann constant
        const scalar kb = Foam::constant::physicoChemical::k.value();

        // const scalar s0 = 216*muc_*kb*T_/(sqr(Foam::constant::mathematical::pi)*rhoF_*pow5(dp)*sqr(ibClass.getGeomModel().getRhoS().value() / rhoF_)*cc);
        // Info << "Brownian force: s0 = " << s0 << endl;
        // scalar f = ibClass.getGeomModel().getM()*sqrt(Foam::constant::mathematical::pi*s0/dt);

        const scalar zeta = 6.0 * Foam::constant::mathematical::pi * muc_ * a;
        const scalar sigma = sqrt( 2.0 * kb * T_ * zeta / dt ) * stability_;

        // To generate a spherical distribution:

        // const scalar theta = rndGen_.scalar01()*Foam::constant::mathematical::twoPi;
        // const scalar u = 2*rndGen_.scalar01() - 1;

        // const scalar a = sqrt(1 - sqr(u));
        // const vector dir(a*cos(theta), a*sin(theta), u);

        vector f = vector
        (
            sigma * rndGen_.scalarNormal(),
            sigma * rndGen_.scalarNormal(),
            sigma * rndGen_.scalarNormal()
        );

        // Info << "Brownian force: " << f << endl;

        return f;
    };

    brownianForce
    (
        const scalar lambda,
        const scalar T,
        const scalar muc,
        const scalar rhoF,
        const scalar stability
    )
    :
    lambda_(lambda),
    T_(T),
    muc_(muc),
    rhoF_(rhoF),
    stability_(stability),
    rndGen_(0)
    {
        if (lambda_ <= 0 || T_ <= 0)
        {
            WarningInFunction
                << "Brownian force is not active due to non-positive parameters."
                << endl;
            return;
        }

        Info << "Brownian force is active with lambda = " << lambda_ << " and T = " << T_ << endl;
    };
    brownianForce(const brownianForce&) = delete;
    brownianForce& operator=(const brownianForce&) = delete;
    ~brownianForce() = default;

    private:

    scalar lambda_;
    scalar T_;
    scalar muc_;
    scalar rhoF_;
    scalar stability_;
    Random rndGen_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
