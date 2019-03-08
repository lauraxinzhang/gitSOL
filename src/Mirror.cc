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


Mirror::Mirror(double xlim, double ylim, double zlim, int nx, double Buniform)
		: Ti_(Ti), Te_(Te),
		  xlim_(xlim), ylim_(ylim), zlim_(zlim), nx_(nx), 
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
	VecDoub * newgrid = new VecDoub(nx);
	for (int i = 0; i < nx_; i++){
		(*newgrid)[i] = 0;
	}
	potential_ = newgrid;
	return;
}

void Mirror::setPotential(double Rratio)
{
	double phiMid = findPhiMid(Rratio);

	VecDoub * newgrid = new VecDoub(nx);
	int imid = (nx_ - 1) / 2;

	(*newgrid)[imid] = phiMid;
	dphi = (phiMid / 2) / (nx_ - 1);

	for (int i = imid; i >= 0; i--){
		// linear profile for now
		(*newgrid)[i]          = (*newgrid)[i + 1] - dphi; // points to the left of iMid
		(*newgrid)[nx - 1 - i] = (*newgrid)[i - 1] - dphi; // point to the right of iMid
	}
	assert((*newgrid)[0] == (*newgrid)[nx_ - 1]); // ensure potential is symetric
	potential_ = newgrid;
	return;
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


double Mirror::findPhiMid(double Rratio)
{
	//TODO: implement this
	double result(0); //place holder

	return result;
}

Vector Mirror::getE(const Vector& pos)
{
	//TODO implement this
	Vector result();
	return result;
}

Vector Mirror::getB(const Vector& pos)
{
	double x = pos.x();
	Vector result();
	return result;
}

Vector Mirror::getB()
{
	Vector result(Buniform_, 0, 0);
	return result;
}

void Mirror::run(bool spec, double vx, double vy, double vz)
{
	// double mass = MI * (1 - spec) + ME * spec;

	// //initialize an hydrogen ion with energy and direction input by user;
	// double vx = sqrt(energy * ex * EVTOJOULE / mass); // thermal velocity
	// double vy = sqrt(energy * ey * EVTOJOULE / mass);
	// double vz = sqrt(energy * ez * EVTOJOULE / mass);

 //    Vector veli(vx, vy, vz);
 //    Vector posi(0.4 * xlim_ , 0, 0);

 //    Particle part(posi, veli, spec);
	return;
}