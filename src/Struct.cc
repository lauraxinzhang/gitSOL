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

#include "Struct.h"

PassingHelp::PassingHelp(Doub Ti, Doub Te, Doub R)
	:Ti_(Ti), Te_(Te), R_(R), xGrid_(nullptr), RGrid_(nullptr), integralTable_(nullptr)
{
	readTable("./input/I_func_TableData.txt", 250, 50); // 250 points in x, 50 in R
	return;
}

PassingHelp::~PassingHelp()
{
	delete xGrid_;
	delete RGrid_;
	delete integralTable_;
}

void PassingHelp::readTable(std::string input, int nx, int nR)
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

Doub PassingHelp::operator() (Doub phiTilde)
{
	INTERP2D function((*xGrid_), (*RGrid_), (*integralTable_));

	Doub xe = -1 * phiTilde;
	Doub xi = (Te_ / Ti_) * phiTilde;
	assert(abs(xe) <= 25 && abs(xi) <= 25); // make sure the interpolation is not out of range
	assert(R_ <= 5);

	Doub Ie = function.interp(xe, R_);
	Doub Ii = function.interp(xi, R_);

	Doub diff = Ie - pow((ME / MI) * (Ti_ / Te_), 0.5) * Ii;

	//std::cerr << "input: " << phiTilde << " Ie " << Ie << " Ii " << Ii << " current diff: " << diff << std::endl;
	return diff;
}

Psir::Psir(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub z)
		:function_(rG, zG, f), z_(z)
{
	/** interpolation objects need to be initialized under the colon 
	 initializer. Otherwise the data members will be messed up. */
	// function_ = function;
}

Doub Psir::operator() (Doub r)
{
	return function_.interp(r, z_);
}

Psiz::Psiz(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub r)
	:function_(rG, zG, f), r_(r)
{
	// function_ = function;
}

Doub Psiz::operator() (Doub z)
{
	return function_.interp(r_, z);
}

psiLimiter::psiLimiter(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, \
	const VecDoub& rL, const VecDoub& zL)
	: psifull_(rG, zG, f), zOfR_(rL, zL), RHS_(0)
{
	// nothing to do here?
}

void psiLimiter::setRHS(Doub r, Doub z)
{
	RHS_ = psifull_.interp(r, z);
	return;
}

// For a given r on the limiter, return the difference in flux 
// between limiter and given RHS
Doub psiLimiter::operator() (const Doub& r)
{
	Doub zNow( zOfR_.interp(r) );

	Doub psiLimiter( psifull_.interp(r, zNow) );

	return psiLimiter - RHS_;
}

Doub psiLimiter::getZ(Doub& r)
{
	return zOfR_.interp(r);
}

pastukhovHelp::pastukhovHelp(Doub Ti, Doub Te, Doub R)
	:Ti_(Ti), Te_(Te)
{
	Doub constant   = sqrt(2.0/PI)*sqrt(MI/ME) * pow( (Ti/Te), (3.0/2.0) );
	// MI and ME defined in particle.cc
	Doub gOfR       = (2.0*R+1.0)/(2.0*R) * log(4.0*R+2.0);
	Doub LHS        =  constant * log10(R)/gOfR;

	LHS_ = LHS;
}

Doub pastukhovHelp::operator()(const Doub& x)
{
	Doub I = 1 + ( sqrt(PI) * exp(x) * erfc( sqrt(x) ) / (2 * sqrt(x)) );
	Doub fOfX = x * exp(x)/I;

	return fOfX - LHS_;
}

eFieldHelp::eFieldHelp(VecDoub lList, VecDoub potList)
		   :helper_(lList, potList){}

Doub eFieldHelp::operator() (Doub l)
{
	return helper_.interp(l);
}

histogram::histogram(Doub min, Doub max, Doub numBin)
		  : min_(min), max_(max), bins_(nullptr)
{
	Doub gridsize = (max_ - min_) / numBin;
	gridsize_ = gridsize;
	VecDoub * bins = new VecDoub(numBin);
	bins_ = bins;
	return;
}

void histogram::addToBin(Doub val)
{
	int index = (int) ( (val - min_ ) / gridsize_);
	(*bins_)[index]++;
	return;
}

