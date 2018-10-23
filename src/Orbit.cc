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
	INTERP1D myfunc(xx, yy); // Set up the interpolation object.

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
	INTERP2D myfunc(xx1, xx2, yy);

	To get an interpolated value:
	Doub x1, x2, y;
	y = myfunc.interp(x1, x2);
*/

Orbit::Orbit(const std::string& field, const std::string& limiter)
	: Phi_(nullptr), Rratio_(nullptr)
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
	std::cerr << "nw:" << nw << "nh:" << nh << std::endl;

	// initialize data containers
	VecDoub * rGrid = new VecDoub(nw);
    VecDoub * zGrid = new VecDoub(nh);

    MatDoub * Br   = new MatDoub(nw_, nh_); 
	MatDoub * Bz   = new MatDoub(nw_, nh_);
	MatDoub * Btor = new MatDoub(nw_, nh_);
	MatDoub * Bmod = new MatDoub(nw_, nh_);

    MatDoub * psiRZ = new MatDoub(nw_, nh_);
    MatDoub * pRZ   = new MatDoub(nw_, nh_); 


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
	if (Rratio_!=nullptr){
		delete Rratio_;
	}
	if (Phi_ != nullptr){
		delete Phi_;
	}

	delete psiRZ_;
	delete pRZ_;


	delete rGrid_;
	delete zGrid_;

	delete rLimit_;
	delete zLimit_;
}

Vector Orbit::getB(const Vector& pos)
{
	// Doub rr = pos.x();

	Doub zz = pos.z();

	Doub rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y());


	INTERP2D Br( (*rGrid_), (*zGrid_), (*Br_));
	INTERP2D Bz( (*rGrid_), (*zGrid_), (*Bz_));
	INTERP2D Btor( (*rGrid_), (*zGrid_), (*Btor_));

	Doub BBr = Br.interp(rr, zz);
	Doub BBz = Bz.interp(rr, zz);
	Doub BBtor = Btor.interp(rr, zz);

	Vector gotB(BBr, BBtor, BBz);
	return gotB;

}

Doub Orbit::getModB(Doub rr, Doub zz)
{
	INTERP2D Br( (*rGrid_), (*zGrid_), (*Br_));
	INTERP2D Bz( (*rGrid_), (*zGrid_), (*Bz_));
	INTERP2D Btor( (*rGrid_), (*zGrid_), (*Btor_));

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
	INTERP1D limiter((*rLimit_), (*zLimit_));

	return limiter.interp(rr);
}

void Orbit::configMirror()
{
	assert(Rratio_ == nullptr);
	MatDoub * Rratio = new MatDoub(nw_, nh_);
	mirrorRatio( *Rratio );
	Rratio_ = Rratio;

	/* this following line doesn't work because Rratio is still a nullptr
	mirrorRatio(*Rratio_);
	*/
	return;
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
				// Now we solve for mirror ratio
    			// bracketing the solutions with the range of limiter
	    		zbrak(limit, rllmtr_, rrlmtr_, 20, xb1, xb2, nroot);

	    		if (nroot == 1 ) {
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

Doub Orbit::getMirrorRatio(Doub rr, Doub zz)
{
	if (Rratio_ == nullptr){
		configMirror();
	}
	INTERP2D mirror( (*rGrid_), (*zGrid_), (*Rratio_) );

	Doub result = mirror.interp(rr, zz);

	return result;
}

Doub Orbit::pastukhov(Doub Ti, Doub Te, Doub R)
{
	if (R <= 1) {
		Doub zero(0);
		return zero;
	}
	if (Rratio_ == nullptr){
		configMirror();
	}
	// The bulk of the calculation preparation is done in struct
	pastukhovHelp help(Ti, Te, R);
    Doub foundX = rtbis(help, 0, 1000, 1E-9);
    return foundX;
}

void Orbit::setPastukhov( Doub Ti, Doub Te, Doub multiplier)
{
	if (Rratio_ == nullptr){
		configMirror();
	}

	Phi_ = new MatDoub(nw_, nh_);
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		if (std::isnan((*Rratio_)[i][j])){
    			(*Phi_)[i][j] = NAN;
    		} else {
	    		Doub mirrorR = (*Rratio_)[i][j];
	    		Doub x = pastukhov(Ti, Te, mirrorR);
	    		(*Phi_)[i][j] = x * multiplier;
	    	}
    	}
    }
    return;
}

void Orbit::fieldLines(std::ofstream &rBList, std::ofstream &zBList, std::ofstream &phiBList, const Vector& init, Doub dl, int iter)
{
    if (!rBList.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    }

    Vector now = init;

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
	if (Rratio_ == nullptr){
		configMirror();
	}

	// Interpolate mirror ratio onto entire RZ plane
	INTERP2D mirrorRZ((*rGrid_), (*zGrid_), (*Rratio_));

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
	// I don't know what this is doing anymore. Why are there trailing zeros in front?	
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
	if (Rratio_ == nullptr){
		configMirror();
	}
	// MatDoub ratio(nw_, nh_);
	// mirrorRatio(ratio);
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Rratio_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::writePastukhov( Doub Ti, Doub Te, std::ofstream &output)
{
	std::cerr << "calling setPastukhov" << std::endl;
	setPastukhov(Ti, Te);
	std::cerr << "writing to file." << std::endl;

	if (!output.is_open()) {
    	std::cerr << "Unable to open output file" << std::endl; 
    	exit(0);
    }

	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
	    	output << (*rGrid_)[i] << ',' << (*zGrid_)[j] << ',' << (*Phi_)[i][j] << std::endl;
    	}
    }
    return;
}

void Orbit::midPlanePhi( Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output)
{
	INTERP2D pasInterp((*rGrid_), (*zGrid_), pas);

	for(Doub r = rllmtr_; r < rrlmtr_; r += 0.002){
		Doub x = pasInterp.interp(r, 0);
		if (!std::isnan(x)) {
			output << r << ',' << x << std::endl;
		}
	}
	return;
}


//----------------------------------------------------------------------
//------------------------ TESTS ---------------------------------------
//----------------------------------------------------------------------

/**
 * \brief Tests for base functions
 *
 */
void Orbit::test()
{
	std::cerr << "Entering tests:" << std::endl;
	bool testVector   = true;
	bool testPart     = true;
	bool testInterp   = true;

	// Flags for which tests to run
	// Vector Tests:
	bool vecInitAndAssign = true;
	bool plusAndMinus     = true;
	bool multAndDiv       = true;
	bool dotAndMod        = true;
	bool cross            = true;

	bool normalize        = true;
	bool parallelAndPerp  = true;

	// Particle Tests:
	bool partInitAndAssign = true;
	bool mover             = true;
	bool testMu                = true;

	// Interpolation Tests:
	if (testInterp){
		std::cerr << " - Testing Linear Interpolations" << std::endl;
		Int n = 10;
		VecDoub xx(n), yy(n);
		for (int i = 0; i < n ; ++i){
			xx[i] = i;
			yy[i] = 2.0 * i + 1;
		}

		Linear_interp line(xx, yy);

		Doub x, y;
		x = 5.5;
		y = line.interp(x);

		bool result = (y == 2*x + 1);

		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}



	// Vector Tests:
	if (testVector && vecInitAndAssign){

		std::cerr << " - Testing Vector init, assignment, and '==' operators" << std::endl;

		Vector v(2, 3, 5);
		Vector v2 = v;
		Vector v3;
		Vector v4(0,0,0);

		if (v == v2 && v3 == v4){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}

	if (testVector && plusAndMinus){
		std::cerr << " - Testing '+' and '-' operators" << std::endl;
		Vector v(2, 3, 5);
		Vector v2(1, 6, 0);

		Vector vPlus(3, 9, 5);
		Vector vMinus(1, -3, 5);

		if (v + v2 == vPlus && v - v2 == vMinus){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
			std::cerr << (v + v2) << "should be" << vPlus << std::endl;
			std::cerr << (v - v2) << "should be" << vMinus << std::endl;

		}

	}

	if (testVector && multAndDiv){
		std::cerr << " - Testing '*' and '/' operators" << std::endl;

		Vector v(2, 3, 5);
		double mult  = 0.7;
		double denom = 0.5;

		Vector vMult(2*0.7, 3*0.7, 5*0.7);
		Vector vDiv(2/0.5, 3/0.5, 5/0.5);

		if (v * mult == vMult && v / denom == vDiv){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
			std::cerr << (v * mult) << "should be" << vMult << std::endl;
			std::cerr << (v / denom) << "should be" << vDiv << std::endl;

		}

	}

	if (testVector && dotAndMod){
		std::cerr << " - Testing dot and mod()" << std::endl;

		Vector vl(2, 7, -8);
		Vector vr(3, 9, 1);

		double vDot = 2*3 + 7*9 + (-8)*1;
		double vlMod = sqrt(4 + 49 + 64);

		bool result = (vl.dot(vr)==vDot) && (vl.mod() == vlMod);

		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	
	}

	if (testVector && cross){
		std::cerr << " - Testing cross product" << std::endl;

		Vector vl(2, 7, -8);
		Vector vr(3, 9, 1);

		Vector vCross(79, -26, -3);
		bool result = (vl.cross(vr) == vCross);
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}
	
	if (testVector && normalize){
		std::cerr << " - Testing normalization" << std::endl;

		Vector vv(4, 3, 0);
		Vector vn(0.8, 0.6, 0);

		bool result = (vv.normalize() == vn);
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}

	if (testVector && parallelAndPerp){
		std::cerr << " - Testing parallel and perp" << std::endl;
		Vector vv(-4, 3, 0);
		Vector B(10, 0, 0);

		Vector vpara(-4, 0, 0);
		Vector vperp(0, 3, 0);

		bool para = (vv.parallel(B) == vpara);
		bool perp = (vv.perp(B) == vperp);
		bool result = para && perp;
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
			std::cerr << "para calculated: " << vv.parallel(B) <<", should be: " << vpara << std::endl;
			std::cerr << "perp calculated: " << vv.perp(B) <<", should be: " << vperp << std::endl;
		}
	}


	// Particle Tests:
	if (testPart && partInitAndAssign){
		std::cerr << " - Testing Particle Init and assign" << std::endl;
		
		Vector pos(0,0,0);
		Vector vel(0,0,0);
		bool species = true;

		Particle p;
		Particle p2(pos, vel, species);

		Vector pos2(2,80.7,0.4);
		Vector vel2(0,8,0);
		bool spec2 = false;

		Particle p3(pos2, vel2, spec2);
		Particle p4 = p3;

		Particle p5;
		p5.setPos(pos2);
		p5.setVel(vel2);
		p5.setSpec(spec2);

		bool result = (p == p2)&& (p4 == p3) && (p5 == p3);

		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}

	if (testPart && mover){
		std::cerr << " - Testing Particle mover" << std::endl;

		Vector E(0, 2.0, 0);
		Vector B(0, 0, 5.0);
		double dt = 0.1;

		Vector posk(0, 0, 0);
		Vector velk(2.0, 2.0, 0);
		Particle part(posk, velk, false);

		part.move(E, B, dt);

		Vector posf(-0.12, -0.19999999, 0.0);
		Vector velf(-1.2, -1.9999999, 0.0);

		Vector posDiff = part.pos() - posf;
		Vector velDiff = part.vel() - velf;
		bool posB = ( posDiff.mod() < 1E-6);
		bool velB = ( velDiff.mod() < 1E-6);
		bool result = posB && velB;

		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;

			std::cerr << part.vel() << "v.s." << velf << std::endl;

			std::cerr << part.pos() << "v.s." << posf << std::endl;
		}

	}

	if (testPart && testMu) {
		std::cerr << " - Testing Particle magnetic moment" << std::endl;

		Vector B(0.5, 0, 0); // B has mod 0.5.
		Vector pos(0, 0, 0); // position is not needed.
		Vector vel(0, 0, 1); // v is entirely perpendicular, has magnitude 1

		Particle part(pos, vel, false);
		double muCorrect = part.mass();

		bool result = (muCorrect == part.mu(B));
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}

	}

}

// A constructor for a simple field configuration
Orbit::Orbit()
	 :nw_(50), nh_(50), rdim_(2), zdim_(2), rleft_(1), zmid_(1)
	 // simag_(0),sibry_(0),dsi_(0)
{
	// Initialize all private data members
	// VecDoub * rlim = new VecDoub(1);
	// VecDoub * zlim = new VecDoub(1); // we don't need these

    VecDoub * rGrid = new VecDoub(nw_);
    VecDoub * zGrid = new VecDoub(nh_);
    // VecDoub * siGrid = new VecDoub(1); // we won't need siGrid.

    const Doub zero = 0;

    MatDoub * Br   = new MatDoub(nw_, nh_, zero); // Br and Bz should be filled with 0s by default
	MatDoub * Bz   = new MatDoub(nw_, nh_, zero);
	MatDoub * Btor = new MatDoub(nw_, nh_);

	// dummy declaration to avoid memory issues.
	MatDoub * psiRZ = new MatDoub(1, 1);
	psiRZ_ = psiRZ;
	MatDoub * pRZ   = new MatDoub(1, 1);
	pRZ_   = pRZ;

	VecDoub * rLimit = new VecDoub(1);
	VecDoub * zLimit = new VecDoub(1);

	rLimit_ = rLimit;
	zLimit_ = zLimit;


	MatDoub * Bmod = new MatDoub(nw_, nh_);
	Bmod_ = Bmod; // Temperary

	// Fill in values for rGrid and zGrid
	for (int i = 0; i < nw_; ++i){     // row index
    	(*rGrid)[i] = rleft_ + i * (rdim_/nw_);
    }

    for (int j = 0; j < nh_; ++j){ // column index
		(*zGrid)[j] = (zmid_ - zdim_/2.0) + j * (zdim_/nh_);
	}

    rGrid_ = rGrid;
    zGrid_ = zGrid;

    Doub bMax = 0.7;

    // // Fill Btor with a simple 1/r dependence w/ maximum 0.2 Tesla.
    // for (int i = 0; i < nw_; ++i){
    // 	for (int j = 0; j < nh_; ++j){
    // 		Doub r = (*rGrid)[i];
    // 		(*Btor)[i][j] = bMax/r;
    // 	}
    // }

    // // A linearly decreasing field
    // for (int i = 0; i < nw_; ++i){
    // 	for (int j = 0; j < nh_; ++j){
    // 		Doub r = (*rGrid)[i];
    // 		(*Btor)[i][j] = bMax - 0.1 * r;
    // 	}
    // } 

    // A constant field
    for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		Doub r = (*rGrid)[i];
    		(*Btor)[i][j] = bMax;
    	}
    } 

    // // A one sided magnetic mirror
    // for (int i = 0; i < nw_; ++i){
    // 	for (int j = 0; j < nh_; ++j){
    // 		Doub r = (*rGrid)[i];
    // 		(*Br)[i][j] = 4 / r;
    // 	}
    // }    
   

    Br_   = Br;
	Bz_   = Bz;
	Btor_ = Btor;

	return;
}

// testing function, only to be used with default constructor.
void Orbit::emptytest()
{
	// run with default constructor for testing
	// print a particle orbit and magnetic field data to file

    ofstream rList;
    rList.open ("./output/rList.out");

    ofstream zList;
    zList.open ("./output/zList.out");

    ofstream phiList;
    phiList.open("./output/phiList.out");

    Vector posi(0.2, 0, zmid_);

    //initialize an hydrogen ion with Ti = 100 eV;
	Doub vv = sqrt(2.0 * 400 * EVTOJOULE / MI);
	std::cerr << vv << std::endl;
    Vector veli(-vv/10, 0, vv);

    Particle part(posi, veli, false);

    Doub dt = 10E-11;
    // Vector EField(0, 0.0005, 0.05);
    Vector EField(0, 0, 0);


    for (int step = 0; step < 1000000; ++step)
    {
    	Vector posNow = part.pos();
    	if (part.pos().x() < 10E-6){
    		std::cerr << "hit the wall" << std::endl;
    		break;
    	}
    	Vector BNow = getB(posNow);

    	// std::cerr << "BField:" << BNow << std::endl;
    	part.moveCyl(EField, BNow, dt);
    	// part.move(EField, BNow, dt);

    	Doub rNow = part.pos().x();
    	Doub phiNow = part.pos().y();
    	Doub zNow = part.pos().z();

    	rList   << std::setprecision(10) << rNow   << std::endl;
    	phiList << std::setprecision(10) << phiNow << std::endl;
    	zList   << std::setprecision(10) << zNow   << std::endl;
    }

    rList.close();
    phiList.close();
    zList.close();
	
	/*
	Vector Bstart( orbit.rleft_, 0, orbit.zmid_ );
	Vector dv(0.1, 0, 0.1);


	ofstream rBList;
    rBList.open ("./output/rBList.bout");

    ofstream zBList;
    zBList.open ("./output/zBList.bout");

    ofstream phiBList;
    phiBList.open("./output/phiBList.bout");

	for (int i = 0; i < 5; i++){
		Vector Bnow = Bstart + dv * i;
		// std::cerr << Bnow << std::endl;
	    orbit.fieldLines(rBList, zBList, phiBList, Bnow, 0.0001, 100000);
	}

	rBList.close();
	zBList.close();
	phiBList.close();
	*/
}


