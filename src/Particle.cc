/**
 * \file    Particle.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implements the Particle class
 * 
 */

#define ME 9.1E-31                        // electron mass, in kilogram
#define MI 1.6726219E-27                  // ion mass, in kilogram
#define QE 1.60217662E-19                 // elementary charge, in coulomb

#include "Particle.h"

Particle::Particle()
		:pos_(), vel_(), species_(true), lost_(false)
{
	charge_ = -1 * QE;
	mass_ = ME;
} 


Particle::Particle(const Vector& pos, const Vector& vel, bool spec)
		:pos_(pos.x(), pos.y(), pos.z()),
		 vel_(vel.x(), vel.y(), vel.z()),
		 species_(spec),lost_(false), generator_(int(time(NULL)))
{
	if (species_){
		charge_ = -1 * QE;
		mass_ = ME;
	} else {
		charge_ = 1 * QE;
		mass_ = MI;
	}
}

Particle::Particle(const Vector& pos, const Vector& vel, int mass, int charge)
		:mass_(mass * MI), charge_(charge * QE),
		 pos_(pos.x(), pos.y(), pos.z()),
		 vel_(vel.x(), vel.y(), vel.z()),
		 lost_(false), generator_(int(time(NULL)))
{

}

Vector Particle::pos() const
{
	return pos_;
}

Vector Particle::vel() const
{
	return vel_;
}

bool Particle::spec() const
{
	return species_;
}

double Particle::mass() const
{
	return mass_;
}

void Particle::setPos(const Vector& right)
{
	pos_.setX(right.x());
	pos_.setY(right.y());
	pos_.setZ(right.z());
	return;
}

void Particle::setVel(const Vector& right)
{
	vel_.setX(right.x());
	vel_.setY(right.y());
	vel_.setZ(right.z());
	return;
}

void Particle::setSpec(const bool right)
{
	species_ = right;
	return;
}

void Particle::lost()
{
	lost_ = true;
	return;
}

bool Particle::isLost()
{
	return lost_;
}

double Particle::mu(Vector& BField) const
{
	Vector vperp = vel().perp(BField);
	double perp = vperp.mod();
	double B = BField.mod();

	double mu = (mass() * perp * perp)/ (2 * B);

	return mu;
}

bool Particle::operator==(const Particle& right)
{
	bool posB = (pos() == right.pos());
	bool velB = (vel() == right.vel());
	bool specB = (spec() == right.spec());
	return (posB && velB && specB);
}

void Particle::move(Vector& E, Vector& B, const double dt)
{
	double qprime;
	
	qprime = dt * charge_ / (2.0 * mass_);

	Vector h = B * qprime;
	// std::cerr << "h:" << h << std::endl;

	double hMod = h.mod();
	// std::cerr << "hMod:" << hMod << std::endl;

	Vector s = (h * 2.0)/(1.0 + pow(hMod,2));
	// std::cerr << "s:" << s << std::endl;

	Vector u = vel_ + E * qprime;
	Vector uuh = u + (u.cross(h));
	// std::cerr << "uuh:" << uuh << std::endl;

	Vector uPrime = u + uuh.cross(s);
	// std::cerr << "u':" << uPrime << std::endl;

	vel_ = uPrime + E * qprime;
	pos_ = pos_ + vel_ * dt; // used updated velocity
	// std::cerr << pos_ << std::endl;
	return;
}

/* THIS IS WRONG - DIFFERENTIAL GEOMETRY WAS NOT HANDLED CORRECTLY

void Particle::moveCyl(Vector& E, Vector& B, const double dt)
{

	double oldR     = pos_.x();
	double oldTheta = pos_.y();
	double oldZ     = pos_.z();

	Vector newpos(oldR, 0, oldZ);

	setPos(newpos);
	move(E, B, dt);

	double newX = pos_.x();
	double newY = pos_.y();
	double newZ = pos_.z();

	// convert pos vector back to cylindrical coordinates

	double newR   = sqrt(pow(newX, 2) + pow(newY, 2));
	// changed from oldR to newX in the following line - 11/14/17

	double dTheta = atan( newY / newX );

	Vector newPos(newR, oldTheta + dTheta, newZ);
	// std::cerr << newPos << std::endl;

	// rotate velocity vector to align with current reference frame:
	double ux = vel_.x();
	double uy = vel_.y();
	double uz = vel_.z();

	double vx = cos(dTheta) * ux - sin(dTheta) * uy;
	double vy = sin(dTheta) * ux + cos(dTheta) * uy;


	Vector newVel(vx, vy, uz);

	setVel(newVel);
	setPos(newPos);

	return;
}

*/

void Particle::moveCyl(Vector& E, Vector& B, const double dt)
{
	/* The particle pusher is now entirely in cartesian coordinates. No
	artificial rotations of any kind. Toroidal motions should already be
	resolved in the 3D Boris algorithm. All fields (accelerations) are 
	converted to Cartesian before applied to the pusher.
	*/

	// Now assume that particle position and velocity vectors are in xyz

	// Magnetic field is still in R, phi, Z
	// Electric field is solved in RZ plane, forcing rotational
	// symmetry, therefore no field in phi direction. We'll keep electic
	// field in R, phi, Z coordinates too.

   	Vector Enew;
   	Vector Bnew;

   	E.cyl2Cart(pos(), Enew);
   	B.cyl2Cart(pos(), Bnew);
   	
   	// The good old Boris push.
   	move(Enew, Bnew, dt);

   	return;
}

void Particle::scatter(Vector& B)
{
	std::bernoulli_distribution distribution(0.5);
	bool sign = distribution(generator_); // choose the orientation of the new vector by random;
	Vector newv = vel().turn(B, sign);
	setVel(newv);
	return;
}
