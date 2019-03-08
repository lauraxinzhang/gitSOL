/**
 * \file    Pusher.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March 2019
 *
 * \brief   Declares the Pusher class.
 *
 *  
 */

#ifndef PUSHER_H_INCLUDED
#define PUSHER_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <stdexcept>
#include <cassert>
#include <sstream>
#include <fstream>

#include "Particle.h"
#include "Vector.h"
#include "Orbit.h"
#include "Particle.h"


template <class T> // templated to be used with Orbit or Mirror class geometry
class Pusher{

	public:
		Pusher(T& geo);
		~Pusher();

		/**
		 * \brief Pushes a single particle, returns final position
		 * \param part The particle to be pushed
		 * \param dt time step
		 * \param iter Number of iterations
		 * \bool write Whether the trajectory should be written to file
		 */
		Vector pushSingle(Particle& part, double dt, int iter, bool write, std::ofstream& coord);

		/**
		 * \brief Single particle pusher in cylindrical geometry
		 * \note only part.move is changed for now. Maybe there's more?
		 */
		Vector pushSingleCyl(Particle& part, double dt, int iter, bool write, std::ofstream& coord);

		// functions for NBI simulations

		/**
		 * \brief burst out particles from spherical surface
		 * \param radius Radius of curvature for surface source
		 */
		void gridBurst(double radius, double ylim, int nsources, bool write);

		// void conicBurst(Vector& pos, Vector norm, double dtheta);

	private:
		T * geo_; // a pointer to the geometry class object

};

template <class T>
Pusher<T>::Pusher(T& geo)
		:geo_(nullptr)
{
	geo_ = &geo;
	return;
}

template <class T>
Pusher<T>::~Pusher()
{
	//delete geo_;
	return;
}

template <class T>
Vector Pusher<T>::pushSingle(Particle& part, double dt, int iter, bool write, std::ofstream& coord)
{
	//coord.open("coordinates.out");
	//coord << std::setprecision(10);
	for (int i = 0; i < iter; i++){
		Vector posNow = part.pos();
    	Vector BNow = (*geo_).getB(posNow);
    	Vector ENow = (*geo_).getE(posNow);

    	part.move(ENow, BNow, dt);

    	if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
    		std::cerr << "particle lost to limiter after" << i \
    		<< "iterations." << std::endl;
    		break;
    	}

    	if (write && (i % 1 == 0)){
    		coord << part.pos() << std::endl;
    	}
	}

	return part.pos();
}

template <class T>
Vector Pusher<T>::pushSingleCyl(Particle& part, double dt, int iter, bool write, std::ofstream& coord)
{
	for (int i = 0; i < iter; i++){
		Vector posNow = part.pos();
    	Vector BNow = (*geo_).getB(posNow);
    	Vector ENow = (*geo_).getE(posNow);

    	part.moveCyl(ENow, BNow, dt);

    	if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
    		std::cerr << "particle lost to limiter after" << i \
    		<< "iterations." << std::endl;
    		break;
    	}

    	if (write && (i % 100 == 0)){
    		coord << part.pos() << std::endl;
    	}
	}

	return part.pos();
}

template <class T>
void Pusher<T>::gridBurst(double radius, double ylim, int nsources, bool write)
{
	std::ofstream coord;
	coord.open("coordBurst.out");
	coord << std::setprecision(10);

	std::default_random_engine generator(int(time(NULL)));
	std::uniform_real_distribution<double> distribution(-1 * ylim, ylim);

	for (int isource = 0; isource < nsources; isource++){
		double yRand = distribution(generator);
		double zRand = distribution(generator);
		double xCalc = -1*sqrt(radius * radius - yRand * yRand - zRand * zRand) + radius;

		Vector posi(xCalc, yRand, zRand);
		Vector veli((-1*xCalc + radius),-1* yRand, -1*zRand);

		std::cerr << posi << std::endl;
		std::cerr << veli << std::endl;

		Particle part(posi, veli, 1, 0); // one particle per source for now

		pushSingle(part, 0.01, 1000, 1, coord);
	}

	return;
}

#endif // PUSHER_H_INCLUDED
