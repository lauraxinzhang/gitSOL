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
Pusher::Pusher(const T& geo)
		:geo_(nullptr)
{
	geo_ = geo;
	return;
}

Pusher::~Pusher()
{
	delete geo_;
	return;
}

Vector Pusher::pushSingle(Particle& part, double dt, int iter, bool write)
{
	ofstream coord;
	coord.open("coordinates.out");
	coord << std::setprecision(10);

	for (int i = 0; i < iter; i++){
		Vector posNow = part.pos();
    	Vector BNow = getB(posNow);
    	Vector ENow = getE(posNow);

    	part.move(ENow, BNow, dt);

    	if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
    		std::cerr << "particle lost to limiter after" << step \
    		<< "iterations." << std::endl;
    		break;
    	}

    	if (write && (step % 100 == 0)){
    		coord << part.pos() << std::endl;
    	}
	}

	return part.pos();
}

Vector Pusher::pushSingleCyl(Particle& part, double dt, int iter, bool write)
{
	ofstream coord;
	coord.open("coordinates.out");
	coord << std::setprecision(10);

	for (int i = 0; i < iter; i++){
		Vector posNow = part.pos();
    	Vector BNow = getB(posNow);
    	Vector ENow = getE(posNow);

    	part.moveCyl(ENow, BNow, dt);

    	if ((*geo_).isLimiter(part.pos())){ // TODO write this method in Mirror
    		std::cerr << "particle lost to limiter after" << step \
    		<< "iterations." << std::endl;
    		break;
    	}

    	if (write && (step % 100 == 0)){
    		coord << part.pos() << std::endl;
    	}
	}

	return part.pos();
}
