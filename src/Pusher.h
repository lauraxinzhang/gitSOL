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
		Pusher(const T& geo);
		~Pusher();

		/**
		 * \brief Pushes a single particle, returns final position
		 * \param part The particle to be pushed
		 * \param dt time step
		 * \param iter Number of iterations
		 * \bool write Whether the trajectory should be written to file
		 */
		Vector pushSingle(Particle& part, double dt, int iter, bool write);

		/**
		 * \brief Single particle pusher in cylindrical geometry
		 * \note only part.move is changed for now. Maybe there's more?
		 */
		Vector pushSingleCyl(Particle& part, double dt, int iter, bool write);

	private:
		T * geo_; // a pointer to the geometry class object

};

template <class T>
Pusher<T>::Pusher(const T& geo)
		:geo_(nullptr)
{
	geo_ = geo;
	return;
}

template <class T>
Pusher<T>::~Pusher()
{
	delete geo_;
	return;
}

template <class T>
Vector Pusher<T>::pushSingle(Particle& part, double dt, int iter, bool write)
{
	std::ofstream coord;
	coord.open("coordinates.out");
	coord << std::setprecision(10);

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

    	if (write && (i % 100 == 0)){
    		coord << part.pos() << std::endl;
    	}
	}

	return part.pos();
}

template <class T>
Vector Pusher<T>::pushSingleCyl(Particle& part, double dt, int iter, bool write)
{
	std::ofstream coord;
	coord.open("coordinates.out");
	coord << std::setprecision(10);

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



#endif // PUSHER_H_INCLUDED
