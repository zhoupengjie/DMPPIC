/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "KinematicParcel.H"
#include "forceSuSp.H"
#include "integrationScheme.H"
#include "meshTools.H"
#include "cloudSolution.H"

#include "GlobalVariable.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::KinematicParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    tetIndices tetIs = this->currentTetIndices();

    td.rhoc() = td.rhoInterp().interpolate(this->coordinates(), tetIs);

    if (td.rhoc() < cloud.constProps().rhoMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed density in cell " << this->cell()
                << " to " << cloud.constProps().rhoMin() <<  nl << endl;
        }

        td.rhoc() = cloud.constProps().rhoMin();
    }

    td.Uc() = td.UInterp().interpolate(this->coordinates(), tetIs);

    td.muc() = td.muInterp().interpolate(this->coordinates(), tetIs);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::calcDispersion
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() = cloud.dispersion().update
    (
        dt,
        this->cell(),
        U_,
        td.Uc(),
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::calcUCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    this->UCorrect_ = Zero;

    this->UCorrect_ =
        cloud.dampingModel().velocityCorrection(p, dt);

    this->UCorrect_ +=
        cloud.packingModel().velocityCorrection(p, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() += cloud.UTrans()[this->cell()]/massCell(td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(td);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Su, dUTrans, Spu);

    this->U_ += this->UCorrect_;

    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;

        // Update momentum transfer coefficient
        cloud.UCoeff()[this->cell()] += np0*Spu;
    }
}


template<class ParcelType>
template<class TrackCloudType>
const Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
) const
{
    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const typename TrackCloudType::forceType& forces = cloud.forces();

    // Momentum source due to particle forces
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, mu);
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, mu);
    const scalar massEff = forces.massEff(p, ttd, mass);

    // Calculate Brownian force

    // static const scalar kB = 1.38064852e-23;  // Boltzmann constant
    // static const scalar T_kelvin = 298.15; // 25 celsius
    const scalar diffusion_coefficient = 1e-9;
    
    const scalar gamma = 6 * M_PI * (d_/2);

    const scalar amplitude = sqrt(2 * diffusion_coefficient * gamma * gamma / dt);

    vector FBrownian(
        rnd.GaussNormal<scalar>() * amplitude,
        rnd.GaussNormal<scalar>() * amplitude,
        rnd.GaussNormal<scalar>() * amplitude
    );

    /*
    // Proper splitting ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Sp()*td.Uc() + Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;
    const scalar bncp = Fncp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaUcp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, acp, bcp
        );
    const vector deltaUncp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, ancp, bncp
        );
    const vector deltaT = deltaUcp + deltaUncp;
    */

    // Shortcut splitting assuming no implicit non-coupled force ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaU = cloud.UIntegrator().delta(U_, dt, acp + ancp, bcp);
    const vector deltaUncp = ancp*dt;
    const vector deltaUcp = deltaU - deltaUncp;

    // The Brownian-induced velocity increment is computed as:
    // deltaU_B = (1/massEff)*sqrt( 2*(kBT)^2*dt / D ) * randomVector,
    // where each component of the random vector is drawn from N(0,1)
    // const scalar kBT = 1.380649e-23 * 298.15;
    // const scalar D = GV.duffusion_coefficient;
    // const scalar amplitude_std = sqrt( 2.0 * (kBT * kBT) * dt / D );

    // const vector deltaU_B = (1.0/massEff) * amplitude_std * vector(
    //     rnd.GaussNormal<scalar>(),
    //     rnd.GaussNormal<scalar>(),
    //     rnd.GaussNormal<scalar>()
    // );

    // print massEff, amplitude_std, deltaU_B;
    // Info << "massEff: " << massEff << 
    // " amplitude_std: " << amplitude_std << 
    // " amplitude: " << (1.0/massEff) * amplitude_std <<
    // " deltaU_B: " << deltaU_B << endl;
    // --- End Brownian contribution ---

    // Calculate the new velocity and the momentum transfer terms
    // vector Unew = U_ + deltaU + deltaU_B;
    vector Unew = U_ + deltaU;

    dUTrans -= massEff*deltaUcp;

    Spu = dt*Fcp.Sp();

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = cloud.pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    UCorrect_(p.UCorrect_)
{}


template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const KinematicParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    UCorrect_(p.UCorrect_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    auto& p = static_cast<typename TrackCloudType::parcelType&>(*this);
    auto& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    ttd.switchProcessor = false;
    ttd.keepParticle = true;

    const cloudSolution& solution = cloud.solution();
    const scalarField& cellLengthScale = cloud.cellLengthScale();
    const scalar maxCo = solution.maxCo();

    const scalar diffusion_std = sqrt(2.0 * GV.duffusion_coefficient  * trackTime * 1.0/0.29999); // compensate step fraction

    while (ttd.keepParticle && !ttd.switchProcessor && p.stepFraction() < 1)
    {
        // Cache the current position, cell and step-fraction
        const point start = p.position();
        const scalar sfrac = p.stepFraction();

        // Brownian motion displacement
        vector s_D
        (
            diffusion_std * distribution(generator),
            diffusion_std * distribution(generator),
            diffusion_std * distribution(generator)
        );

        // vector s_D(
        //     diffusion_std * rnd.GaussNormal<scalar>(),
        //     diffusion_std * rnd.GaussNormal<scalar>(),
        //     diffusion_std * rnd.GaussNormal<scalar>()
        // );


        // Total displacement over the time-step
        const vector s = trackTime * U_ + s_D;
        // const vector s = trackTime * U_;

        // Cell length scale
        const scalar l = cellLengthScale[p.cell()];

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        // Fraction of the displacement to track in this loop. This is limited
        // to ensure that the both the time and distance tracked is less than
        // maxCo times the total value.
        scalar f = 1 - p.stepFraction();
        f = min(f, maxCo);
        f = min(f, maxCo*l/max(SMALL*l, mag(s)));
        if (p.active())
        {
            // Track to the next face
            p.trackToFace(f*s - d, f);
        }
        else
        {
            // At present the only thing that sets active_ to false is a stick
            // wall interaction. We want the position of the particle to remain
            // the same relative to the face that it is on. The local
            // coordinates therefore do not change. We still advance in time and
            // perform the relevant interactions with the fixed particle.
            p.stepFraction() += f;
        }

        const scalar dt = (p.stepFraction() - sfrac)*trackTime;

        // Avoid problems with extremely small timesteps
        if (dt > ROOTVSMALL)
        {
            // Update cell based properties
            p.setCellValues(cloud, ttd);

            p.calcDispersion(cloud, ttd, dt);

            if (solution.cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(cloud, ttd, dt);
            }

            p.calcUCorrection(cloud, ttd, dt);

            p.calc(cloud, ttd, dt);
        }

        p.age() += dt;

        if (p.active() && p.onFace())
        {
            ttd.keepParticle = cloud.functions().postFace(p, ttd);
        }

        ttd.keepParticle = cloud.functions().postMove(p, dt, start, ttd);

        if (p.active() && p.onFace() && ttd.keepParticle)
        {
            p.hitFace(s, cloud, ttd);
        }
    }

    return ttd.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicParcel<ParcelType>::hitPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    auto& p = static_cast<typename TrackCloudType::parcelType&>(*this);
    auto& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    // Invoke post-processing model
    td.keepParticle = cloud.functions().postPatch(p, pp, ttd);

    if (isA<processorPolyPatch>(pp))
    {
        // Skip processor patches
        return false;
    }
    else if (cloud.surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // Surface film model consumes the interaction, i.e. all
        // interactions done
        return true;
    }
    else
    {
        // This does not take into account the wall interaction model
        // Just the polyPatch type. Then, a patch type which has 'rebound'
        // interaction model will count as escaped parcel while it is not
        if (!isA<wallPolyPatch>(pp) && !polyPatch::constraintType(pp.type()))
        {
            cloud.patchInteraction().addToEscapedParcels(nParticle_*mass());
        }

        // Invoke patch interaction model
        return cloud.patchInteraction().correct(p, pp, td.keepParticle);
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    TrackCloudType&,
    trackingData&
)
{
    // wall interactions are handled by the generic hitPatch method
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "KinematicParcelIO.C"

// ************************************************************************* //
