/**
 * \file    Orbit.thermal.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    November 2019
 *
 * \brief   Implements the Orbit class, toy model for particle orbits in LTX
 *          SOL w/ Pastukhov potential solved analytically.
 */


#include "Orbit.h"

Doub Orbit::getTe(Doub rr, Doub zz)
{
	return Te_;
}

Doub Orbit::getNe(Doub rr, Doub zz)
{
	return 1E19;
}

Doub Orbit::getTi(Doub rr, Doub zz)
{
	return Ti_;
}

Doub Orbit::getNi(Doub rr, Doub zz)
{
	return 1E19;
}

Doub Orbit::chandrasekharG(Doub x)
{
	Doub gauss = ((2 * x) / sqrt(PI)) * exp(-1 * pow(x, 2));
	return (erf(x) - gauss)/(2 * pow(x, 2));
}


void Orbit::xAndNu_ab(Particle& part, Particle& partB, Doub& xB, Doub& nu_ab)
{
	Doub dens, temp, ma, ea, mb, eb;
	ma = part.mass();
	ea = part.charge();
	mb = partB.mass();
	eb = partB.charge();

	Doub rr = sqrt(part.pos().x() * part.pos().x() + part.pos().y() * part.pos().y() );
	Doub zz = part.pos().z();
	if (partB.spec()){ // if an electron
		dens = getNe(rr, zz);
		temp = getTe(rr, zz) * QE; // convert to joules
	} else { // if background proton
		dens = getNi(rr, zz);
		temp = getTi(rr, zz) * QE;
	}
	// First return value
	nu_ab =  dens * pow(ea, 2) * pow(eb, 2) * COULOMBLOG / (4 * PI * pow(EPSILON0, 2) * pow(ma, 2));

	// Second return value
	Doub vtB = sqrt(2 * temp / mb);
	xB = part.speed() / vtB;

	return;
}


Doub Orbit::nu_s(Particle& part, Particle& partB, Doub xB, Doub nu_ab)
{
	Doub ma, mb, temp;
	Doub rr = sqrt(part.pos().x() * part.pos().x() + part.pos().y() * part.pos().y() );
	Doub zz = part.pos().z();

	if (partB.spec()){ // if an electron
		temp = getTe(rr, zz) * QE;
	} else { // if background proton
		temp = getTi(rr, zz) * QE;
	}

	ma = part.mass();
	mb = partB.mass();

	Doub g = chandrasekharG(xB);
	return g * (ma + mb) * nu_ab / (part.speed() * temp);
}

Doub Orbit::nu_D(Particle& part, Doub xB, Doub nu_ab)
{
	return nu_ab * (erf(xB) - chandrasekharG(xB)) / pow(part.speed(), 3);
}

Doub Orbit::nu_para(Particle& part, Doub xB, Doub nu_ab)
{
	return 2 * nu_ab * chandrasekharG(xB) / pow(part.speed(), 3);
}

void Orbit::collisions(Particle& part, Doub dt, std::default_random_engine& generator)
{
	Doub xe, nu_ae, xp, nu_ap;
	Particle electron;
	Particle proton;
	proton.setSpec(0);

	assert(proton.mass() != electron.mass() && proton.charge() == -1 * electron.charge());
	xAndNu_ab(part, electron, xe, nu_ae);
	xAndNu_ab(part, proton, xp, nu_ap);

	Doub v = part.speed();
	// slowing down on both species
	Doub nu_slow = nu_s(part, electron, xe, nu_ae) + nu_s(part, proton, xp, nu_ap);
	// perpendicular and parallel diffusion coeeficients on both electron and protons
	Doub D_D = pow(v, 2) * 0.5 * (nu_D(part, xp, nu_ap) + nu_D(part, xe, nu_ae));
	Doub D_para = pow(v, 2) * 0.5 * (nu_para(part, xp, nu_ap) + nu_para(part, xe, nu_ae));

	// update particle velocity: slowing down
	Vector vslow = part.vel() * nu_slow * dt;
	
	Matrix identity;
	identity.diagonal(1, 1, 1);

	Vector vel = part.vel();
	Matrix parallel = vel.tensor(vel) * pow(v, -2);
	Matrix perp     = identity - parallel;
	Matrix diffusion = parallel * sqrt(D_para) + perp * sqrt(D_D);

	std::normal_distribution<double> distribution(0, sqrt(dt)); // generate a Gaussian distributed velocity
    double wx = distribution(generator); // generate 3 normal distributed velocities.
    double wy =  distribution(generator);
    double wz = distribution(generator);
    Vector wiener(wx, wy, wz);

    Vector diffused = diffusion.dot(wiener);

    Vector updated = part.vel() - vslow + diffused;
	part.setVel(updated);
	return;
}
