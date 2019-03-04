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


Mirror::Mirror(double xlim, double ylim, double zlim, double dx, double Buniform)
		: Ti_(Ti), Te_(Te),
		  xlim_(xlim), ylim_(ylim), zlim_(zlim), dx_(dx), 
		  xGrid_(nullptr), potential_(nullptr),
		  Buniform_(Buniform)
{
	// nothing else to do.
	return;
}

Mirror::~Mirror()
{
	// delete all dynamically allocated objects.
	delete xGrid_;
	delete potential_;
	return;
}

Vector Mirror::getE(const Vector& pos)
{
	//TODO implement this
	Vector result();
	return result;
}

Vector Mirror::getB(const Vector& pos)
{
	Vector result();
	return result;
}

Vector Mirror::getB()
{
	Vector result(Buniform_);
	return result;
}