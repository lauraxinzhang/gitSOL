/**
 * \file    Orbit.struct.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March, 2019
 *
 * \brief   Implements the Orbit class nested structs.
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

#include "Orbit.h"

Orbit::passingHelp::passingHelp(Doub Ti, Doub Te, Doub R)
	:Ti_(Ti), Te_(Te), R_(R), xGrid_(nullptr), RGrid_(nullptr), integralTable_(nullptr)
{
	readTable("./input/integral.csv", 100, 50); // 100 points in x, 50 in R
	return;
}

void Orbit::passingHelp::readTable(std::string input, int nx, int nR)
{
	std::ifstream inputfile;
	inputfile.open(input);

	VecDoub * xGrid = new VecDoub(nx);
	VecDoub * RGrid = new VecDoub(nR);
	MatDoub * integral  = new MatDoub(nx, nR);
	Doub current, x, R;

	for (int i = 0; i < nR; i++){
		inputfile >> R;
		(*RGrid)[i] = R;
		for (int j = 0; j < nx; j++){
			// assumes each line of input is R, x, I
			inputfile >> x >> current;
			(*integral)[i][j] = current;
			// (*RGrid)[i] = R;
			(*xGrid)[j] = x;
		
			if (j != nx - 1){
				Doub Rnext;
				inputfile >> Rnext;
				assert(Rnext == R); // making sure we're indeed still on the same R
			}
		}
	}

	xGrid_ = xGrid;
	RGrid_ = RGrid;
	integralTable_ = integral;
	return;
}

Doub Orbit::passingHelp::operator() (Doub phiTilde)
{
	INTERP2D function((*xGrid_), (*RGrid_), (*integralTable_));

	Doub Ie = function.interp(-1 * phiTilde, R_);
	Doub Ii = function.interp( (Te_ / Ti_) * phiTilde, R_);

	Doub diff = Ie - (ME / MI) * (Ti_ / Te_) * Ii;

	return diff;
}
