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
	readTable("./input/I_func_TableData.txt", 250, 50); // 250 points in x, 50 in R
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
			(*integral)[j][i] = current;
			// (*RGrid)[i] = R;
			(*xGrid)[j] = x;
			Doub Rnext;
			if (j != nx - 1){
				
				inputfile >> Rnext;
			//	std::cerr << j << ','<< R<<',' << Rnext << std::endl;
				assert(Rnext == R); // making sure we're indeed still on the same R
				
			}
			// std::cerr << j << ','<< R<<',' << Rnext << std::endl;
		}
	}

	xGrid_ = xGrid;
	RGrid_ = RGrid;
	integralTable_ = integral;
	//std::cerr << "input parsing complete" << std::endl;
	return;
}

Doub Orbit::passingHelp::operator() (Doub phiTilde)
{
	INTERP2D function((*xGrid_), (*RGrid_), (*integralTable_));

	Doub xe = -1 * phiTilde;
	Doub xi = (Te_ / Ti_) * phiTilde;
	assert(abs(xe) <= 25 && abs(xi) <= 25); // make sure the interpolation is not out of range
	
	Doub Ie = function.interp(xe, R_);
	Doub Ii = function.interp(xi, R_);

	Doub diff = Ie - (ME / MI) * (Ti_ / Te_) * Ii;

	std::cerr << "input: " << phiTilde << " Ie " << Ie << " Ii " << Ii << " current diff: " << diff << std::endl;
	return diff;
}
