/**
 * \file    Particle.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Declares the Particle class
 * 
 */

#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED 1

#include <iostream>
#include <time.h> // for seeding the random engine.
#include <random>
#include "Vector.h"
#include "Constants.h"

/**
 * \brief  A class for individual particle objects
 *
 *
 */
class Particle
{
    public:

    	/**
    	 * \brief Default constructor. 
    	 * \details Calls the default constructor for Vector class to initialize both
    	 *          pos and vel to [0,0,0]. Species default to electron.
    	 */
    	Particle();   // default constructor
    	// ~Particle();


    	/**
    	 * \brief Constructor for Class Particle
    	 * \param pos  position of the particle
    	 * \param vel  velocity of the particle
    	 * \param spec species of the particle - 0 for hydrogen, 1 for electron
    	 */
    	Particle(const Vector& pos, const Vector& vel, bool spec);

        /**
         * \brief Another constructor for particles other than proton and electrons
         * \param pos    position of particle
         * \param vel    velocity of particle
         * \param mass   Mass number of particle (in proton masses)
         * \param charge Charge number of particle
         *
         */
        Particle(const Vector& pos, const Vector& vel, int mass, int charge);



    	/**
    	 * \brief Data getter for position
    	 */
    	Vector pos() const;

    	/**
    	 * \brief Data getter for velocity
    	 */
    	Vector vel() const;

    	/**
    	 * \brief Data getter for species
    	 */
    	bool spec() const;

        /**
         * \brief Data getter for mass
         */
        double mass() const;

    	double charge() const;

	double energy() const;

	void setPos(const Vector& right);
    	void setVel(const Vector& right);
    	void setSpec(const bool right);

        /**
         * \brief Set the lost_ flag to true if the particle is lost.
         */
        void lost();

        /**
         * \brief Return whether the particle is lost
         */
        bool isLost();

        /**
         * \brief Calculate magnetic moment at the given B field
         * \note  B Field value must be the field at current position
         */
        double mu(Vector& BField) const;

    	bool operator==(const Particle& right);

    	void move(Vector& E, Vector& B, const double dt);
        void moveCyl(Vector& E, Vector& B, const double dt);

        /**
         * \brief Turn the particle pitch angle by 90 degrees
         *
         */
        void scatter(Vector& B);

    private:
        double   mass_;
        double   charge_;

    	Vector pos_; // pos at time step k
    	Vector vel_; // vel at time step k-1/2
    	bool   species_;

        bool   lost_; // always start with false. Change to true when limiter is reached.
        std::default_random_engine generator_; // a random number generator, seeded at Particle construction;
};

#endif // PARTICLE_H_INCLUDED
