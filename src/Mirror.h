/**
 * \file    Mirror.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    Feb 2019
 *
 * \brief   Declares the Mirror class, a straight mirror field background for iteratively
 *          solving the electrostatic potential.
 * 
 * \note    Uses the same Particle and Vector classes to perform basic operations
 */

#ifndef MIRROR_H_INCLUDED
#define MIRROR_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <stdexcept>

// Programmer class includes
#include "Vector.h"
#include "Particle.h"

// Physical contants. Don't change unless you find the value to be wrong.
// TODO: if both Orbit and Mirror are compiled, are these repeated defs going to give problems?

#define PI 3.14159265358979323846		 // pi, no explanation needed...
#define CC 299792458.0                   // speed of light
// #define QE 4.8E-10                    // electron charge, in Fr

#define ME 9.1E-31                       // electron mass, in kilogram
#define MI 1.6726219E-27                 // ion mass, in kilogram
#define QE 1.60217662E-19                // elementary charge, in coulomb
#define EVTOJOULE 1.60217662E-19         // energy conversion factor

class Mirror
{
	public:
		Mirror();

	private:
};


#endif // MIRROR_H_INCLUDED