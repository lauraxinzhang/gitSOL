/**
 * \file    Constants.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March, 2019
 *
 * \brief   Declares the structs for root finding.
 *
 * \NOTE    See header file for detailed explainations for each 
 *          function \n
 *
 *          Uses Numerical Recipe nr3.h data containers for Vector and Matrix
 *          functionalities; uses Vector.h for 3D vectors with arithmetic operations.
 *      
 *          Uses nr3 routines for interpolations.     
 *      
 *          Uses double instead of float to avoid floating number errors.
 */

#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED 1

#include "nr3.h"
#include "interp_1d.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "dfridr.h"
#include "roots.h"

// Physical contants. Don't change unless you find the value to be wrong.
#define PI 3.14159265358979323846		 // pi, no explanation needed...
#define CC 299792458.0                   // speed of light
// #define QE 4.8E-10                    // electron charge, in Fr

#define ME 9.1E-31                       // electron mass, in kilogram
#define MI 1.6726219E-27                 // ion mass, in kilogram
#define QE 1.60217662E-19                // elementary charge, in coulomb
#define EVTOJOULE 1.60217662E-19         // energy conversion factor
#define EPSILON0 8.85E-12                // F/m SI units.


// Numbers specific to LTX geometry
#define BMAGAXIS  2000                   // mod(B) = 0.2 T = 2000 Gauss at magnetic axis, characteristic field strength
#define NPERORBIT 20                     // steps per Lamor orbit, from Boris convergence test.
#define RMAJOR    0.4  

#define COULOMBLOG 19                    // coulomb logarithm for collisions

/** Choose which interpolater to use throughout the class */
// Use these two lines to compile with linear interpolations
//
typedef Linear_interp INTERP1D;
typedef Bilin_interp  INTERP2D;

// Use these two lines to compile with Spline interpolations 
// (very slow. use only in production mode when very pretty pictures are needed)
//typedef Spline_interp    INTERP1D;
//typedef Spline2D_interp  INTERP2D;


#endif
