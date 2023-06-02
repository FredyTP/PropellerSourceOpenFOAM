#include "bemTargetCoeff.H"
#include "bladeElementModel.H"
#include "unitConversion.H"
namespace Foam 
{
defineTypeNameAndDebug(bemTargetCoeff,0);
addToRunTimeSelectionTable(bemControl,bemTargetCoeff, dictionary);

bemTargetCoeff::bemTargetCoeff(const dictionary &dict, const bladeElementModel& bem)
:
    bemControl(dict),
    bem_(bem),
    calcFrequency_(-1),
    target_(Zero),
    theta_(Zero),
    nIter_(50),
    tol_(1e-8),
    relax_(1.0),
    dTheta_(degToRad(0.1))
{
    read(dict);
}

void bemTargetCoeff::read(const dictionary &dict)
{
    dict.readEntry("thrust", target_[0]);
    dict.readEntry("pitch", target_[1]);
    dict.readEntry("roll", target_[2]);

    theta_[0] = degToRad(dict.get<scalar>("theta0Ini"));
    theta_[1] = degToRad(dict.get<scalar>("theta1cIni"));
    theta_[2] = degToRad(dict.get<scalar>("theta1sIni"));

    dict.readEntry("calcFrequency", calcFrequency_);

    dict.readIfPresent("nIter", nIter_);
    dict.readIfPresent("tol", tol_);
    dict.readIfPresent("relax", relax_);

    if (dict.readIfPresent("dTheta", dTheta_))
    {
        dTheta_ = degToRad(dTheta_);
    }

    angularVelocity_ = rotorControl::readAngularVelocity(dict);
}

void bemTargetCoeff::correctControl(const vectorField &U, const scalarField *rhoField)
{
    if (bem_.mesh().time().timeIndex() % calcFrequency_ == 0)
    {
        word calcType = "forces";

        Info<< type() << ":" << nl
            << "    solving for target trim " << calcType << nl;

        const scalar rhoRef = bem_.rhoRef();

        // iterate to find new pitch angles to achieve target force
        scalar err = GREAT;
        label iter = 0;
        tensor J(Zero);

        vector old = Zero;
        while ((err > tol_) && (iter < nIter_))
        {
            // cache initial theta vector
            vector theta0(theta_);

            // set initial values
            old = calcForces(U,rhoField);

            // construct Jacobian by perturbing the pitch angles
            // by +/-(dTheta_/2)
            for (label pitchI = 0; pitchI < 3; pitchI++)
            {
                theta_[pitchI] -= dTheta_/2.0;
                vector cf0 = calcForces(U, rhoField);

                theta_[pitchI] += dTheta_;
                vector cf1 = calcForces(U, rhoField);

                vector ddTheta = (cf1 - cf0)/dTheta_;

                J[pitchI + 0] = ddTheta[0];
                J[pitchI + 3] = ddTheta[1];
                J[pitchI + 6] = ddTheta[2];

                theta_ = theta0;
            }

            // calculate the change in pitch angle vector
            vector dt = inv(J) & (target_ - old);

            // update pitch angles
            vector thetaNew = theta_ + relax_*dt;

            // update error
            err = mag(thetaNew - theta_);
            // update for next iteration
            theta_ = thetaNew;
            iter++;
        }

        if (iter == nIter_)
        {
            Info<< "    solution not converged in " << iter
                << " iterations, final residual = " << err
                << "(" << tol_ << ")" << endl;
        }
        else
        {
            Info<< "    final residual = " << err << "(" << tol_
                << "), iterations = " << iter << endl;
        }

        Info<< "    current and target " << calcType << nl
            << "        thrust  = " << old[0]*rhoRef << ", " << target_[0] << nl
            << "        pitch   = " << old[1]*rhoRef << ", " << target_[1] << nl
            << "        roll    = " << old[2]*rhoRef << ", " << target_[2] << nl
            << "    new pitch angles [deg]:" << nl
            << "        theta0  = " << radToDeg(theta_[0]) << nl
            << "        theta1c = " << radToDeg(theta_[1]) << nl
            << "        theta1s = " << radToDeg(theta_[2]) << nl
            << endl;
    }
}


vector bemTargetCoeff::calcForces(const vectorField &U, const scalarField *rhoField) const
{
    vector forces;
    auto result = bem_.calculate(U,rhoField);
    forces[0] = result.thrust();
    forces[1] = result.torque.x();
    forces[2] = result.torque.y();
    return forces;
}


scalar bemTargetCoeff::getAzimuth(scalar azimuth0) const
{
    return azimuth0;
}

scalar bemTargetCoeff::getPitch(scalar azimuth) const
{
    return PitchFunction(azimuth,theta_[0],theta_[1],theta_[2]);
}

scalar bemTargetCoeff::getFlapping(scalar azimuth) const
{
    return 0;
}

scalar bemTargetCoeff::getAzimuthDot(scalar azimuth) const
{
    return angularVelocity_;
}
scalar bemTargetCoeff::getPitchDot(scalar azimuth) const
{
    return angularVelocity_*(-theta_[2]*cos(azimuth) + theta_[1]*sin(azimuth));
}
scalar bemTargetCoeff::getFlappingDot(scalar azimuth) const
{
    return 0;
}

scalar bemTargetCoeff::getAngularVelocity() const
{
    return angularVelocity_;
}


}




