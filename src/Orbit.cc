/**
 * \file    Orbit.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implements the Orbit class, toy model for particle orbits in LTX
 *          SOL w/ Pastukhov potential solved analytically.
 *
 * \NOTE    See header file for detailed explainations for each 
 *          function \n
 *
 *          Uses Numerical Recipe nr3.h data containers for Vector and Matrix
 *          functionalities; uses Vector.h for 3D vectors with arithmetic operations.
 *      
 *          Uses nr3 routines for interpolations.     
 *
 *          Reads data from g-eqdsk files.
 *      
 *          Uses double instead of float to avoid floating number errors.
 */


#include "Orbit.h"

/*
	To use linear interpolations:

	 (See also: test case)

	Int n = length of vector to be interpolated.
	VecDoub xx(n), yy(n);
	....
	Linear_interp myfunc(xx, yy); // Set up the interpolation object.

	// To get an interpolated value:
	Doub x, y;
	y = myfuncinterp(x);

	// Int, Doub, and other capitalized data types are defined in nr3.h
*/


/*
	To use 2D bilinear interpolations:

	Int m = ..., n = ...;
	MatDoub yy(m, n);       // function values defined on grid points
	VecDoub xx1(m), xx2(n);	// coordinates of grid points
	...
	Bilin_interp myfunc(xx1, xx2, yy);

	To get an interpolated value:
	Doub x1, x2, y;
	y = myfunc.interp(x1, x2);
*/

Orbit::Orbit(const std::string& field, const std::string& limiter)
{
	// Read and configure field here, 
	// read limiter coordinates in helper function.

	ifstream myfile;
    myfile.open(field); // open .dat file

    if (!myfile.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    } else {
    	std::cerr << "Uniform Grid file successfully opened" << std::endl;
    }

    std::string a;

    while (myfile >> a){
    	// truncates all the crap in front until a number is reached. 
    	// 
    	if (std::isdigit(a[0])){
    		int dum = std::stoi(a);
    		if (dum == 0){ break; }
    	}	    
	}

	// LTX G-S data file: R, Z, Br, Bt, Bz, P, Psi
	double nw, nh;
	myfile >> std::setprecision(10)>> nw >> nh;
	nw_ = nw;
	nh_ = nh;

	// initialize data containers
	VecDoub * rGrid = new VecDoub(nw);
    VecDoub * zGrid = new VecDoub(nh);

    MatDoub * Br   = new MatDoub(nw_, nh_); 
	MatDoub * Bz   = new MatDoub(nw_, nh_);
	MatDoub * Btor = new MatDoub(nw_, nh_);
	MatDoub * Bmod = new MatDoub(nw_, nh_);

    MatDoub * psiRZ = new MatDoub(nw_, nh_);
    MatDoub * pRZ   = new MatDoub(nw_, nh_); 

    // // dummy declarations to avoid memory issues.
    // VecDoub * siGrid = new VecDoub(1);
    // siGrid_ = siGrid;


	// fill in data containers with field, pressure, and flux values;
	for (int ir = 0; ir < nw; ir++){
		double rNow;
		myfile >> std::setprecision(10) >> rNow;

		(*rGrid)[ir] = rNow;

		for (int iz = 0; iz < nh; iz++){
			double zNow;
			double BrNow;
			double BtNow;
			double BzNow;
			double PNow;
			double PsiNow;

			myfile >> std::setprecision(10) >> zNow >> BrNow >> BtNow >> BzNow >> PNow >> PsiNow;

			// transfer everything to the matrix
			// ignore pressure for now.
			(*zGrid)[iz] = zNow;

			(*Br)[ir][iz] = BrNow;
			(*Btor)[ir][iz] = BtNow;
			(*Bz)[ir][iz] = BzNow;

			Vector Bvec(BrNow, BtNow, BzNow);
			Doub mod = Bvec.mod();
			(*Bmod)[ir][iz] = mod;

			(*psiRZ)[ir][iz] = PsiNow;
			(*pRZ)[ir][iz]   = PNow;

			if (iz != nh - 1){
				// pop the R value from the next line
				double rNext;
				myfile >> std::setprecision(10) >> rNext;
				// Check and make sure we're still in the same r.
				assert(rNext == rNow);
			}
		}
	}

	rGrid_ = rGrid;
	zGrid_ = zGrid;

	Br_    = Br;
	Bz_    = Bz;
	Btor_  = Btor;
	Bmod_  = Bmod;

	psiRZ_ = psiRZ;
	pRZ_ = pRZ;

	rdim_ = (*rGrid_)[nw - 1] - (*rGrid_)[0];
	zdim_ = (*zGrid_)[nh - 1] - (*zGrid_)[0];
	rleft_ = (*rGrid_)[0];
	zmid_ = 0;

	//DONE with reading .dat file.

	// Read from limiter file now.

	ifstream limit;
	limit.open(limiter);
	readLimiter(limit);

	if (limit.is_open()){
		limit.close();
	}

	return;

}	

void Orbit::readLimiter(std::ifstream &input)
{
	if (!input.is_open()) {
    	std::cerr << "Unable to open limiter file" << std::endl; 
    } else {
    	std::cerr << "limiter location file successfully opened" << std::endl;
    }

    std::string a;

    while (input >> a){
    	// truncates all the crap in front until a number is reached. 
    	// 
    	if (std::isdigit(a[0])){
    		int dum = std::stoi(a);
    		if (dum == 0){ break; }
    	}	    
	}


	VecDoub * rLimit = new VecDoub(70); // coord file has 70 lines right now, hard coded in.
	VecDoub * zLimit = new VecDoub(70);

	Doub rnow, znow;
	for (int row = 0; row < 70; ++row){
		input >> rnow >> znow;
		// std::cerr << rnow << "," << znow << std::endl;
		(*rLimit)[row] = rnow * 0.0254;
		(*zLimit)[row] = znow * 0.0254;
	}
	rLimit_ = rLimit;
	zLimit_ = zLimit;

	rllmtr_ = (*rLimit)[0];
	rrlmtr_ = (*rLimit)[69];

	zllmtr_ = (*zLimit)[0];
	zrlmtr_ = (*zLimit)[69];
	return;
}



Orbit::~Orbit()
{
	// delete rlim_;
	// delete zlim_;
	delete Br_;
	delete Bz_;
	delete Btor_;
	delete Bmod_;

	delete psiRZ_;
	delete pRZ_;


	delete rGrid_;
	delete zGrid_;

	delete rLimit_;
	delete zLimit_;
	// delete siGrid_;
}



Vector Orbit::getB(const Vector& pos)
{
	// Doub rr = pos.x();

	Doub zz = pos.z();

	Doub rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y());


	Bilin_interp Br( (*rGrid_), (*zGrid_), (*Br_));
	Bilin_interp Bz( (*rGrid_), (*zGrid_), (*Bz_));
	Bilin_interp Btor( (*rGrid_), (*zGrid_), (*Btor_));

	Doub BBr = Br.interp(rr, zz);
	Doub BBz = Bz.interp(rr, zz);
	Doub BBtor = Btor.interp(rr, zz);

	Vector gotB(BBr, BBtor, BBz);
	return gotB;

}

Doub Orbit::getModB(Doub rr, Doub zz)
{
	Bilin_interp Br( (*rGrid_), (*zGrid_), (*Br_));
	Bilin_interp Bz( (*rGrid_), (*zGrid_), (*Bz_));
	Bilin_interp Btor( (*rGrid_), (*zGrid_), (*Btor_));

	Doub BBr = Br.interp(rr, zz);
	Doub BBz = Br.interp(rr, zz);
	Doub BBtor = Btor.interp(rr, zz);
	Vector gotB(BBr, BBtor, BBz);

	Doub result = gotB.mod();

	return result;
}

bool Orbit::isLimiter(Doub rr, Doub zz)
{

	psiLimiter limit( (*rGrid_), (*zGrid_), (*psiRZ_), (*rLimit_), (*zLimit_) );

	bool result = false;
	if (rr <= rllmtr_ || rr >= rrlmtr_){
		// our r is outside plasma region
		result = true;
	} else {

		Doub ztlmtr = limit.getZ(rr);

		if (zz >= abs(ztlmtr) || zz <= -abs(ztlmtr)){
			// our z is outside plasma region;
			result = true;
		}
	}
	return result;
}

Doub Orbit::getzLimiter(Doub rr)
{
	Linear_interp limiter((*rLimit_), (*zLimit_));

	return limiter.interp(rr);
}


void Orbit::mirrorRatio(MatDoub& ratio)
{
	psiLimiter limit( (*rGrid_), (*zGrid_), (*psiRZ_), (*rLimit_), (*zLimit_) );

	for (int iz = 0; iz < nh_; ++iz){
		for (int ir = 0; ir < nw_; ++ir){

    		Doub rr = (*rGrid_)[ir];
    		Doub zz = (*zGrid_)[iz];

    		limit.setRHS(rr, zz);

    		VecDoub xb1, xb2;
    		int nroot(0);

    		
    		if ((*pRZ_)[ir][iz] != 0){
    			// we're within last closed surface
    			ratio[ir][iz] = NAN;
    		} else if (isLimiter(rr, zz)){
    			ratio[ir][iz] = NAN;
			} else {
				// std::cerr << "here 3" << std::endl;
				// Now we solve for mirror ratio
    			// bracketing the solutions with the range of limiter
	    		zbrak(limit, rllmtr_, rrlmtr_, 20, xb1, xb2, nroot);

	    		if (nroot == 1 ) {
	    			// std::cerr << rr << "," << zz << std::endl;
	    			Doub r1 = xb1[0];
	    			Doub r2 = xb2[0];

	    			Doub rl = rtbis(limit, r1, r2, 1E-9);
	    			Doub zl = limit.getZ(rl);

	    			Doub Bzero = getModB(rr, zz);
	    			Doub Bmax = getModB(rl, zl);

	    			Doub ratioNow = Bmax/Bzero;
	    			ratio[ir][iz] = ratioNow;

	    		} else {
	    			// when no root or multiple roots are found
	    			ratio[ir][iz] = NAN;
	    		}		    	
	    	}   		
    	}
    }
	return;
}



Doub Orbit::pastukhov(Doub Ti, Doub Te, Doub R)
{
	if (R <= 1) {
		Doub zero(0);
		return zero;
	}
	// The bulk of the calculation preparation is done in struct
	pastukhovHelp help(Ti, Te, R);

    Doub foundX = rtbis(help, 0, 1000, 1E-9);

    return foundX;
}

void Orbit::fieldLines(std::ofstream &rBList, std::ofstream &zBList, std::ofstream &phiBList, const Vector& init, Doub dl, int iter)
{

    if (!rBList.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    }

    Vector now = init;

    // std::cerr << now << std::endl;

	for (int i = 0; i < iter; ++i){
		Doub rNow = now.x();
		Doub phiNow = now.y();
		Doub zNow = now.z();

		if (rNow >= rleft_ + rdim_ || zNow >= zmid_ + (zdim_/2) || zNow <= zmid_ - (zdim_/2)){

			break;
		}

		rBList   << std::setprecision(10) << rNow   << std::endl;
    	phiBList << std::setprecision(10) << phiNow << std::endl;
    	zBList   << std::setprecision(10) << zNow   << std::endl;

    	Vector bNow = getB(now);
		Vector bHat = bNow.normalize();

		Vector dB = bHat * dl;
		now = now + dB;
	}

	return;
}

void Orbit::eField(Doub Ti, Doub Te, Vector init, Doub dl, int iter, std::ofstream &output)
{
	if (!output.is_open()) {
    	std::cerr << "Unable to open output file" << std::endl; 
    }

    // Prepare mirror ratio matrix
    MatDoub ratio(nw_, nh_);
	mirrorRatio(ratio);

	// Interpolate mirror ratio onto entire RZ plane
	Bilin_interp mirrorRZ((*rGrid_), (*zGrid_), ratio);

    Vector  now  = init;
    Doub    lTot = 0;

    VecDoub lList(iter);
    VecDoub potList(iter);

    Doub maxIter(0);

	for (int i = 0; i < iter; ++i){
		Doub rr = now.x();
		Doub zz = now.z();

		std::cout << rr << "," << zz << std::endl;

		if (isLimiter(rr, zz)){
			std::cerr << "limiter reached after " << i << "iterations " << std::endl;
			break;
		} else {

			Doub ratioNow  = mirrorRZ.interp(rr, zz);

			if (std::isnan(ratioNow)){
				std::cerr << "mirror ratio undefined before limiter is reached. Break after " << i << "iterations " << std::endl;
				break;
			}

			Doub potential = pastukhov(Ti, Te, ratioNow);

			Vector bNow = getB(now);
			Vector bHat = bNow.normalize(); // find the direction of the magnetic field

			Vector dBl = bHat * dl;
			now = now + dBl;  // Each step's arc length is exactly dl

			lTot +=dl; // Add to arclength

			lList[i]   = lTot;
			potList[i] = potential;
			// std::cerr << lList[i] << " , " << potList[i]<< " , " << ratioNow << std::endl;
		}
		maxIter++;
	}

	// Done with calculating potential

	if (maxIter < iter){
	// move the elements over to a smaller container to not mess up the interpolation

		VecDoub lListNew(maxIter);
		VecDoub potListNew(maxIter);

		for (int i = 0; i < maxIter; ++i){
			lListNew[i] = lList[i];
			potListNew[i] = potList[i];
		}
		lList = lListNew;
		potList = potListNew;
	}

	
	int zeroHead(0);
	while(potList[zeroHead] == 0){
		++zeroHead;
	}
	potList.erase(potList.begin(), potList.begin() + zeroHead);
	

    VecDoub eList(potList.size());


	// Now calculate electric field via numerical differentiation.

	eFieldHelp phiOfL(lList, potList);
	Doub err;

	for (int i = 0; i < potList.size(); ++i){

		Doub lc = lList[i]; // get current arc location
		if (potList[i] == 0){
			eList[i] = NAN;
		} else {
	    	eList[i] = dfridr(phiOfL, lc, 0.2, err);
	    	// << eList[i] << ", " << err << std::endl;
	    	// std::cerr << lc << ", " << potList[i] << ", "<< phiOfL(lc) << ", " << eList[i] << ", " << err << std::endl;
	    	// And then write results to file.
	    }
    	output << std::setprecision(10) << lc << "," << potList[i] << "," << eList[i] << std::endl;
    }
	return;
}

void Orbit::writeModB(std::ofstream &output)
{
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Bmod_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writeFlux(std::ofstream &output)
{
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*psiRZ_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writeMirrorRatio(std::ofstream &output)
{
	MatDoub ratio(nw_, nh_);
	mirrorRatio(ratio);
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << ratio[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writePastukhov( Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output)
{
	MatDoub ratio(nw_, nh_);
	mirrorRatio(ratio);

	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		if (std::isnan(ratio[i][j])){
    			pas[i][j] = NAN;
    		} else {
	    		Doub mirrorR = ratio[i][j];
	    		Doub x = pastukhov(Ti, Te, mirrorR);
	    		pas[i][j] = x;
	    	}
	    	output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << pas[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::midPlanePhi( Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output)
{
	// MatDoub ratio(nw_, nh_);
	// mirrorRatio(ratio);

	Bilin_interp pasInterp((*rGrid_), (*zGrid_), pas);

	for(Doub r = rllmtr_; r < rrlmtr_; r += 0.002){

		Doub x = pasInterp.interp(r, 0);

		if (!std::isnan(x)) {
			// Doub x   = pastukhov(Ti, Te, rat);
			output << r << ',' << x << std::endl;
		}
	}
	return;
}


