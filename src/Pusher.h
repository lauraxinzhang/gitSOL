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
#include <string>

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

		/**
		 * \brief Burst out cones of particles from spherical surface. Basically a better 
		 *        version of gridBurst
		 * \param radius Radius of curvature for spherical source
		 * \param ylim   half height of source plate
		 * \param dtheta Divergence of gaussian cones
		 * \param nsources Number of particle sources
		 * \param partPerS Number of particles per source
		 * \param write  Whether to write things to file
		 *
		 */
		void conicBurst(double radius, double ylim, double dtheta, int nsources, int partPerS, bool write);

		/**
		 * \brief Generate a random pos vector on sphere surface
		 * \note  Assumes that the canter of the sphere is at (radius, 0, 0)
		 */
		Vector sphere(double radius, double ylim, std::default_random_engine& generator);


		/**
		 * \brief Generate a normal vector on sphere surface, pointing towards center
		 * \note  Assumes that the canter of the sphere is at (radius, 0, 0)
		 */
		Vector sphereNormal(double radius, Vector pos);

		/**
		 * \brief  Find a diverged vector on the conic surface. 
		 *
		 * \note   Calculates normal vector from radius and pos, generate a random vector on the tangent plane,
		 *         then add to normal vector to find the diverged vector.
		 *
		 */
		Vector diverge(double radius, Vector& pos, double dtheta, std::default_random_engine& generator);

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
	return;
}

template <class T>
Vector Pusher<T>::pushSingle(Particle& part, double dt, int iter, bool write, std::ofstream& coord)
{
	//coord.open("coordinates.out");
	//coord << std::setprecision(10);
	int lastCrossed = 26;// start at the last sightline (first to cross)
	for (int i = 0; i < iter; i++){
	//std::cerr << "pushing iteration: " << i << std::endl;	
		Vector posNow = part.pos();
    	Vector BNow = (*geo_).getB(posNow);
    	Vector ENow = (*geo_).getE(posNow);

    	part.move(ENow, BNow, dt);
    	//int lastCrossed = 26;// start at the last sightline (first to cross)

    	if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
    		// std::cerr << "particle lost to limiter after" << i \
    		// << "iterations." << std::endl;
    		break;
    	}
    	if (write){
	    	lastCrossed = (*geo_).sightline(part, lastCrossed);
	    }
    	if (write && (i % 100 == 0)){
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
    		//std::cerr << "particle lost to limiter after" << i \
    		<< "iterations." << std::endl;
    		break;
    	}
    	if (write){
	    	(*geo_).sightline(part);
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
	// std::uniform_real_distribution<double> distribution(-1 * ylim, ylim);

	for (int isource = 0; isource < nsources; isource++){
		// Vector posi(xCalc, yRand, zRand);
		// Vector veli((-1*xCalc + radius), -1* yRand, -1*zRand);
		Vector posi = sphere(radius, ylim, generator);
		Vector veli = sphereNormal(radius, posi);

		//std::cerr << posi << std::endl;
		//std::cerr << veli << std::endl;
		Particle part(posi, veli, 1, 0); // one particle per source for now
		pushSingle(part, 0.01, 1000, write, coord);
	}

	return;
}

template <class T>
void Pusher<T>::conicBurst(double radius, double ylim, double dtheta, int nsources, int partPerS, bool write)
{
	std::ofstream conic;
	conic.open("conicBurst.out");
	conic << std::setprecision(10);

	std::default_random_engine generator(int(time(NULL))); // initialize outside loop to avoid overseeding
	
	// std::uniform_real_distribution<double> location(-1 * ylim, ylim); // for particle sources
	// std::normal_distribution<double> pitchAngle(0, dtheta);

	for (int iS = 0; iS < nsources; iS++){
		Vector posi = sphere(radius, ylim, generator);

		for (int n = 0; n < partPerS; n++){
			Vector veli = diverge(radius, posi, dtheta, generator);
			//std::cerr << posi << std::endl;
			//std::cerr << veli << std::endl;
			Particle part(posi, veli, 1, 0);
			pushSingle(part, 0.001, 2000, write, conic);
		}
	}
	return;
}

template <class T>
Vector Pusher<T>::sphere(double radius, double ylim, std::default_random_engine& generator)
{
	std::uniform_real_distribution<double> distribution(-1 * ylim, ylim);

	double yRand = distribution(generator);	
	double zRand = distribution(generator);
	while (yRand * yRand + zRand * zRand >= ylim * ylim){
		// if it's not in the circle, try again.
		yRand = distribution(generator);
        zRand = distribution(generator);
    }
	double xCalc = -1 * sqrt(radius * radius - yRand * yRand - zRand * zRand) + radius;

	Vector posi(xCalc, yRand, zRand);
	return posi;
}

template <class T>
Vector Pusher<T>::sphereNormal(double radius, Vector pos)
{
	double x = radius - pos.x();
	double y = -1 * pos.y();
	double z = -1 * pos.z();
	Vector result(x, y, z);
	return result;
}

template <class T>
Vector Pusher<T>::diverge(double radius, Vector& pos, double dtheta, std::default_random_engine& generator)
{
	std::normal_distribution<double> pitchAngle(0, dtheta);
	std::uniform_real_distribution<double> uni(0, 1);

	Vector posNorm = pos.normalize();
	double x0 = posNorm.x();
	double y0 = posNorm.y();
	double z0 = posNorm.z();

	double b = uni(generator);
	double c = uni(generator);
	double a = (radius * x0 - b * y0 - c * z0) / (x0 - radius);

	Vector tangent(a - x0, b - y0, c - z0); // a randomly generated vector tangent to sphere at pos
	Vector velNorm = sphereNormal(radius, pos).normalize(); // normal vector of sphere at pos

	double theta = pitchAngle(generator);
	Vector vperp = tangent.normalize() * tan(theta);

	Vector result = velNorm + vperp;

	return result.normalize();
}



#endif // PUSHER_H_INCLUDED
