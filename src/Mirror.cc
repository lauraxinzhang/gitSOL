/**
 * \file    Mirror.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    Feb 2019
 *
 * \brief   Implements the Mirror class, a straight mirror field background for iteratively
 *          solving the electrostatic potential.
 *
 *  
 */

#include "Mirror.h"

Mirror::Mirror(double Ti, double Te, double Buniform, double R, double L, int nx, )
		: Ti_(Ti), Te_(Te), Buniform_(Buniform), R_(R), L_(L),
		  nx_(nx), 
		  xGrid_(nullptr), xShift_(nullptr), potential_(nullptr), EField_(nullptr)
{
	// make sure nx is odd
	assert(nx % 2 == 1);

	// define computation box according to physical paramters
	xlim_ = 3 * L;
	ylim_ = L;
	zlim_ = L;

	// fill in x grid
	double dx = 2 * xlim_ / (nx_ - 1); // extends to a xlim on each side.
	VecDoub * newgrid = new VecDoub(nx);
	(*newgrid)[0] = -1 * xlim_; 
	for (int i = 0; i < nx_; i++){
		(*newgrid)[i] = (*newgrid)[i - 1] + dx;
	}
	xGrid_ = newgrid;
	assert(xGrid_ != nullptr && xGrid_->size() == nx_);

	setGridShift();
	assert(xShift_ != nullptr);

	//fill in potential grid
	setPotential(R_);
	// ensure that the potential grid is a valid vector
	assert(potential_!= nullptr);

	return;
}

Mirror::~Mirror()
{
	// delete all dynamically allocated objects.
	delete xGrid_;
	delete xShift_;
	delete potential_;
	return;
}

void Mirror::setGridShift()
{
	if(xShift_ == nullptr){ // if shifted grids haven't been created yet
		// create shifted x grids (right shift)
		VecDoub * xShift = new VecDoub( xGrid_ -> begin(), --xGrid_ -> end() );
		assert(xShift->size() == nx_ - 1);

		double dx = 2 * xlim_ / (nx_ - 1); // extends to a xlim on each side.
		for (int i = 0; i < nx_ - 1; ++i){
			(*xShift)[i] += dx / 2; // shift by half a grid size
		}
		xShift_ = xShift;
	}
	return;
}

void Mirror::setPotential()
{
	std::cerr << "warning: potential is set to be 0 everywhere." << std::endl;
	VecDoub * newgrid = new VecDoub(nx_);
	for (int i = 0; i < nx_; i++){
		(*newgrid)[i] = 0;
	}
	potential_ = newgrid;
	return;
}

void Mirror::setPotential(double Rratio)
{
	double phiMid = findPhiMid(Rratio);

	VecDoub * newgrid = new VecDoub(nx_);

	// initialize with a cosine shape:
	// phi = phiMid * cos(x * pi / 2 * xlim)

	for (int i = 0; i < nx_; i++){
		double x = (*xGrid_)[i] / xlim_;
		double phi = phiMid * cos(0.5 * PI * x);
	}

	assert((*newgrid)[0] == (*newgrid)[nx_ - 1]); // ensure potential is symetric
	potential_ = newgrid;
	return;
}


double Mirror::findPhiMid(double Rratio)
{
	PassingHelp help(Ti_, Te_, Rratio);
	Doub upper = 24.9 * Ti_ / Te_;
	Doub result = rtbis(help, 0, upper, 1E-9); // root bracketed between 0 and 10, required by input file.
	return result;
}

void Mirror::setEField()
{
	VecDoub * EField = new VecDoub(xShift_-> size());
	if (potential_ == nullptr){
		setPotential(R_);
	} else {
		for (int i = 0; i < xShift_-> size(); i++){
			double dphi = (*potential_)[i + 1] - (*potential_)[i];
			double dx   = 2 * xlim_ / (nx_ - 1);
			double E    = -1 * dphi / dx;
			(*EField)[i] = E;
		}
	}
	EField_ = EField;
	return;
}

Vector Mirror::getE(const Vector& pos)
{
	if (EField_ == nullptr){
		setEField();
	}
	INTERP1D electric((*xShift_), (*EField_));
	double field = electric.interp(pos.x());
	Vector result(field, 0, 0);
	return result;
}

Vector Mirror::getB(const Vector& pos)
{
	double x = pos.x() - (xlim_ / 2); // shift to coordinates of magnetic field description.
	double y = pos.y();

	double denom = R_ + pow( (x / L_), 4);
	double Bx = Buniform_ * ( 1 + pow( (x / L_), 4) ) / denom;
	double By = Buniform_ * 4 * ( 1 - R_ ) * pow( (x / L_), 3) * ( y / L_ ) / pow(denom, 2);
	Vector result(Bx, By, 0);
	return result;
}

Vector Mirror::getB()
{
	Vector result(Buniform_, 0, 0);
	return result;
}

bool Mirror::isLimiter(const Vector& pos)
{
	bool result = false;
	// hits the limit of computation box
	if (abs(pos.x()) >= xlim_ || abs(pos.y()) >= ylim_ || abs(pos.z()) >= zlim_){
		result = true;
	}
	return result;
}


