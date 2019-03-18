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

Mirror::Mirror(double xlim, double ylim, double zlim, int nx, double Buniform)
		: xlim_(xlim), ylim_(ylim), zlim_(zlim), nx_(nx), 
		  xGrid_(nullptr), potential_(nullptr),
		  Buniform_(Buniform)
{
	// make sure nx is odd
	assert(nx % 2 == 1);
	// fill in x grid
	double dx = xlim_ / (nx_ - 1);
	VecDoub * newgrid = new VecDoub(nx);
	(*newgrid)[0] = 0;
	for (int i = 0; i < nx_; i++){
		(*newgrid)[i] = (*newgrid)[i - 1] + dx;
	}
	xGrid_ = newgrid;

	//fill in potential grid
	setPotential();
	// ensure that the potential grid is a valid vector
	assert(potential_!= nullptr);

	return;
}

Mirror::~Mirror()
{
	// delete all dynamically allocated objects.
	delete xGrid_;
	delete potential_;
	return;
}

void Mirror::setPotential()
{
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
	int imid = (nx_ - 1) / 2;

	(*newgrid)[imid] = phiMid;
	double dphi = (phiMid / 2) / (nx_ - 1);

	for (int i = imid; i >= 0; i--){
		// linear profile for now
		(*newgrid)[i]          = (*newgrid)[i + 1] - dphi; // points to the left of iMid
		(*newgrid)[nx_ - 1 - i] = (*newgrid)[i - 1] - dphi; // point to the right of iMid
	}
	assert((*newgrid)[0] == (*newgrid)[nx_ - 1]); // ensure potential is symetric
	potential_ = newgrid;
	return;
}


double Mirror::findPhiMid(double Rratio)
{
	//TODO: implement this root finding goes here
	double result(0); //place holder

	return result;
}

void Mirror::setEField()
{
	if (potential_ == nullptr){
		setPotential();
	} else {
		// do something to setup a grid of electric field. Define it on a shifted grid?
	}
	return;
}


Vector Mirror::getE(const Vector& pos)
{
	//TODO implement this
	Vector result;
	return result;
}

Vector Mirror::getB(const Vector& pos)
{
	double x = pos.x();
	Vector result;
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
	// case 1 : hits the limit of computation box
	if (pos.x() >= xlim_ || abs(pos.y()) >= ylim_ || abs(pos.z()) >= zlim_){
		result = true;
	}
	// case 2: Beam scraper
	return result;
}

int Mirror::sightline(Particle& part, int lastCrossed)
{
	std::string prefix  = "./SLoutput/sl";
	std::string suffix  = ".out";
	ofstream output;

	if (sightlines_.size() == 0){
		std::string path = "./input/viewlineBL.csv";
		setSightlines(path, 26);
	}

	if ( abs( part.pos().z() ) < 0.008 ){ // pos is in the plane of the sightlines
//		std::cerr << "in the right plane" << std::endl;
		for (int i = lastCrossed + 1; i < 26; i++){
		// i is index for rows, goes from 0 to 25
			double yinter = sightlines_[3 * i + 1]/100.0;
			double slope  = sightlines_[3 * i + 2];

			double yexpect = slope * part.pos().x() + yinter;
			//std::cerr << "i = " << i;
			double deltaY = abs( yexpect - part.pos().y() );
			//std::cerr << ", deltaY = " << deltaY << std::endl;
			if ( deltaY < 0.008 ){
				// pos is on sightline numbered i+1
				std::cerr << "crossing sightline #" << i+1 << std::endl;

				Vector sl(1, slope + yinter, 0); // define vector direction for current sightline
				Vector vpara = part.vel().parallel(sl); // find parallel velocity

				double parallel = vpara.dot(sl.normalize());
				output.open(prefix + std::to_string(i+1) + suffix, ios::app);
				output << parallel << std::endl;
				//output.close();
				//output.clear();
				return i;
			}
		}
	}
	return lastCrossed;
}

void Mirror::setSightlines(const std::string& inputSL, int rows)
{
	ifstream input;
	input.open(inputSL);
	if (!input.is_open()) {
    	std::cerr << "Unable to open sightline file" << std::endl; 
    }

    double val;
    for (int row = 1; row <= 3 * rows; row++){
    	input >> val;
    //	std::cerr << "reading val = "<< val << std::endl;
	sightlines_.push_back(val);
    }
    //for (int i = 0; i < 26; i++){
//	std::cerr << sightlines_[3*i + 1] << "," << sightlines_[3 * i + 2] << std::endl;
  //  }
    return;
}

// void Mirror::writeSightLines(int slIndex, Vector& vel)
// {
// 	std::string prefix  = "./SLoutput/sl";
// 	std::string suffix  = ".out";

// 	ofstream output;
// 	output.open(prefix + strd::to_string(slIndex) + suffix);

// 	double slope  = sightlines_[2 * i];
// 	double yinter = sightlines_[2 * i + 1];

// 	return;
// }

void Mirror::run()
{
	// double mass = MI * (1 - spec) + ME * spec;

	// //initialize an hydrogen ion with energy and direction input by user;
	// double vx = sqrt(energy * ex * EVTOJOULE / mass); // thermal velocity
	// double vy = sqrt(energy * ey * EVTOJOULE / mass);
	// double vz = sqrt(energy * ez * EVTOJOULE / mass);

 //    Vector veli(vx, vy, vz);
 //    Vector posi(0.4 * xlim_ , 0, 0);

 //    Particle part(posi, veli, spec);
	//
	return;
}
