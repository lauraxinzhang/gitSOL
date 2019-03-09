/**
 * \file    Pusher.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March 2019
 *
 * \brief   Implements the Pusher class.
 *
 *  
 */

#include "Pusher.h"

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


