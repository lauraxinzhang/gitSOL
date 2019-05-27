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
*/

/*
	To use 2D interpolations:

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
	: Phi_(nullptr), Rratio_(nullptr), thetaRZ_(nullptr), Er_(nullptr), Ez_(nullptr),\
	  rShift_(nullptr), zShift_(nullptr)
{
	// Read and configure field here, 
	// read limiter coordinates in helper function.

	ifstream myfile;
    myfile.open(field); // open .dat file

    if (!myfile.is_open()) {
    	std::cerr << "Unable to open file" << std::endl; 
    } else {
    	// std::cerr << "Uniform Grid file successfully opened" << std::endl;
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
	// std::cerr << "nw:" << nw << "nh:" << nh << std::endl;

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

			myfile >> std::setprecision(10) >> zNow >> BrNow >> BtNow >> \
			BzNow >> PNow >> PsiNow;

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

	dr_ = (*rGrid_)[1] - (*rGrid_)[0];
	dz_ = (*zGrid_)[1] - (*zGrid_)[0];

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
    	// std::cerr << "Unable to open limiter file" << std::endl; 
    } else {
    	// std::cerr << "limiter location file successfully opened" << std::endl;
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
	if (Rratio_!=nullptr) delete Rratio_;
	if (Phi_ != nullptr) delete Phi_;

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
	MatDoub * thetaRZ = new MatDoub(nw_, nh_);

	mirrorRatio( *Rratio , *thetaRZ);
	Rratio_ = Rratio;
	thetaRZ_ = thetaRZ;

	/* this following line doesn't work because Rratio is still a nullptr
	mirrorRatio(*Rratio_);
	*/
	return;
}

void Orbit::mirrorRatio(MatDoub& ratio, MatDoub& thetaRZ)
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
    			thetaRZ[ir][iz] = NAN;
    		} else if (isLimiter(rr, zz)){
    			ratio[ir][iz] = NAN;
    			thetaRZ[ir][iz] = NAN;

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

	    			Doub thetaZero = theta(rr, zz);
	    			Doub thetaMax = theta(rl, zl);
	    			thetaRZ[ir][iz] = thetaZero / thetaMax;
	    		} else {
	    			// when no root or multiple roots are found
	    			ratio[ir][iz] = NAN;
    				thetaRZ[ir][iz] = NAN;

	    		}		    	
	    	}   		
    	}
    }
	return;
}

Doub Orbit::theta(Doub rr, Doub zz)
{
	Doub deltaR = rr - RMAJOR;
	Doub k = zz/deltaR;
	Doub result(0);
	if (deltaR < 0) {
		if (zz < 0){
			result = atan(k) - PI;
		} else {
			result = atan(k) + PI;
		}
	} else {
		result = atan(k);
	}
	return result;
}

Doub Orbit::getMirrorRatio(Doub rr, Doub zz)
{
	if (Rratio_ == nullptr) configMirror();
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
	if (Rratio_ == nullptr) configMirror();
	setTemp(Ti, Te);

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

Doub Orbit::passing(Doub Ti, Doub Te, Doub R)
{
	if (Rratio_ == nullptr){
		configMirror();
	}
	setTemp(Ti, Te);
	PassingHelp help(Ti, Te, R);
	Doub upper = 24.9 * Ti_ / Te_;
	Doub foundX = rtbis(help, 0, upper, 1E-9); // root bracketed between 0 and 10, required by input file.
	return foundX;
}

// void Orbit::setPassing(Doub Ti, Doub Te, Doub multiplier)
// {
// 	if (Rratio_ == nullptr) configMirror();
// 	setTemp(Ti, Te);

// 	Phi_ = new MatDoub(nw_, nh_);
// 	for (int i = 0; i < nw_; ++i){
//     	for (int j = 0; j < nh_; ++j){
//     		if (std::isnan((*Rratio_)[i][j])){
//     			(*Phi_)[i][j] = NAN;
//     		} else {
// 	    		Doub mirrorR = (*Rratio_)[i][j];
// 	    		Doub x = passing(Ti, Te, mirrorR);
// 	    		(*Phi_)[i][j] = x * multiplier;
// 	    	}
//     	}
//     }
// }

void Orbit::setPassing(Doub Ti, Doub Te, Doub multiplier)
{
	VecDoub phiList;
	VecDoub psiList;

	INTERP2D psiRZ(rGrid_, zGrid_, psiRZ_);

	for(Doub r = rllmtr_; r < rrlmtr_; r += 0.002){
		Doub R = getMirrorRatio(r, 0);
		if (!std::isnan(R)) {
			Doub x = passing(Ti, Te, R); // in normalized units e phi/Te
			Doub psi = psiRZ.interp(r, 0);
			phiList.push_back(x);
			psiList.push_back(psi);
		}
	}
	INTERP1D phiOfPsi(psiList, phiList);

	// now fill in the potential
	Phi_ = new MatDoub(nw_, nh_);
	for (int i = 0; i < nw_; ++i){
    	for (int j = 0; j < nh_; ++j){
    		if (std::isnan((*Rratio_)[i][j])){
    			(*Phi_)[i][j] = NAN;
    		} else {
	    		Doub thetaRatio = (*thetaRZ_)[i][j];

	    		Doub psiNow = (*psiRZ_)[i][j];
	    		Doub phiMid = phiOfPsi.interp(psiNow);

	    		Doub phi = phiMid * pow( cos(0.5 * PI * thetaRatio), 6);
	    		(*Phi_)[i][j] = phi * multiplier;
	    	}
    	}
    }
    return;
}

void Orbit::setTemp(Doub Ti, Doub Te)
{
	Ti_ = Ti;
	Te_ = Te; 
}

void Orbit::setEField()
{
	int zero(0);
	// std::cerr << "zero is: " << zero << std::endl;
	MatDoub * Er = new MatDoub(nw_ - 1, nh_, zero); // initialize to 0
	MatDoub * Ez = new MatDoub(nw_, nh_ - 1, zero);
	assert( (*Er)[1][1] == 0);

	Doub phiNow, phiRight, phiUp, dPhidR, dPhidZ;
	// calculate field
	for (int ir = 0; ir < nw_ - 1; ++ir){
		// Doub rNow = rShift[i];
    	for (int iz = 0; iz < nh_ - 1; ++iz){
    		
    		phiNow   = (*Phi_)[ir][iz] * Te_;
    		phiRight = (*Phi_)[ir + 1][iz] * Te_;
    		phiUp    = (*Phi_)[ir][iz + 1] * Te_; // return to standard unit.

    		if ( phiNow != 0 && !std::isnan(phiNow) && !std::isnan(phiRight) && phiRight != 0 ) { // r derivative is defined
	    		dPhidR = (phiRight - phiNow) / dr_;
	    		(*Er)[ir][iz] = -1 * dPhidR; // NEGATIVE GRADIENT
	    	} 
	    	if ( phiNow != 0 && !std::isnan(phiNow) && !std::isnan(phiUp) && phiUp != 0 ){ // z derivative is defined
	    		dPhidZ = (phiUp    - phiNow) / dz_;
	    		(*Ez)[ir][iz] = -1 * dPhidZ;
	    	} 
    	}
    }
    // top row
    for (int ir = 0; ir < nw_ - 1; ++ir){
		phiNow   = (*Phi_)[ir][nh_ - 1] * Te_;
    	phiRight = (*Phi_)[ir + 1][nh_ - 1] * Te_;
    	if ( phiNow != 0 && !std::isnan(phiNow) && !std::isnan(phiRight) && phiRight != 0 ) { // r derivative is defined
    		dPhidR = (phiRight - phiNow) / dr_;
    		(*Er)[ir][nh_ - 1] = -1 * dPhidR; // NEGATIVE GRADIENT
    	} 
    }
    // right-most column
    for (int iz = 0; iz < nh_ - 1; ++iz){
		phiNow   = (*Phi_)[nw_ - 1][iz] * Te_;
		phiUp    = (*Phi_)[nw_ - 1][iz + 1] * Te_; // return to standard unit.
		if ( phiNow != 0 && !std::isnan(phiNow) && !std::isnan(phiUp) && phiUp != 0 ){ // z derivative is defined
    		dPhidZ = (phiUp - phiNow) / dz_;
    		(*Ez)[nw_ - 1][iz] = -1 * dPhidZ;
    	} 
    }

    assert( (*Er)[5][5] == 0);
    Er_ = Er; // try this
    Ez_ = Ez;

    setGridShift();
    return;
}

void Orbit::setGridShift()
{
	if(rShift_ == nullptr){ // if shifted grids haven't been created yet
		// create shifted r-z grids (up - right shift)
		VecDoub * rShift = new VecDoub( rGrid_ -> begin(), --rGrid_ -> end() );
		VecDoub * zShift = new VecDoub( zGrid_ -> begin(), --zGrid_ -> end() );
		assert(rShift->size() == nw_ - 1);
		assert(zShift->size() == nh_ - 1);

		for (int i = 0; i < nw_ - 1; ++i){
			(*rShift)[i] += dr_ / 2; // shift by half a grid size
		}
		for (int j = 0; j < nh_ - 1; ++j){
			(*zShift)[j] += dz_ / 2;
		}
		rShift_ = rShift;
		zShift_ = zShift;
	}
	return;
}

Vector Orbit::getE(const Vector& pos)
{
	if (rShift_ == nullptr) setGridShift();
	if (Er_ == nullptr) setEField();

//	std::cerr << "before interp" << std::endl;
	INTERP2D fieldR((*rShift_),  (*zGrid_ ), *Er_);
	INTERP2D fieldZ((*rGrid_ ),  (*zShift_), *Ez_);

//	std::cerr << "after interp" << std::endl;
	Doub zz = pos.z();
	Doub rr = sqrt(pos.x() * pos.x() + pos.y() * pos.y());

	Doub Er(0), Ez(0); // start with default value 0
	// std::cerr << (*rShift_)[0] << std::endl;

	if ( rr >= (*rShift_).front() && rr <= (*rShift_).back()){
//		std::cerr << "calling interp in r" << std::endl;
		Er = fieldR.interp(rr, zz);
	}
	if ( zz >= (*zShift_).front() && zz <= (*zShift_).back()){
//		std::cerr << "calling interp in z" << std::endl;
		
		Ez = fieldZ.interp(rr, zz);
	}
	Vector gotE(Er, 0, Ez);
	return gotE;
}

Doub Orbit::timeStep(Doub vv)
{
	return 0; // place holder for maybe when this is needed.
}

void Orbit::particlePush(Doub dr, Doub energy, bool spec, Doub er, Doub ephi, Doub ez, Doub mult)
{
	int species = spec;
	std::string prefix  = "./output/";
	std::string suffix = "_E_" + std::to_string(energy) + "_dr_" + \
	std::to_string(dr) + "spec" + std::to_string(species) + "mult" + std::to_string(mult)+ ".out";

	std::string coordRZ = prefix + "coordRZ" + suffix;
	std::string coordXYZ = prefix + "coordXYZ" + suffix;
	std::string totenergy   = prefix + "totalEnergy" + suffix;
	std::string mag      = prefix + "Mu" + suffix;

	ofstream coordinatesRZ;
	coordinatesRZ.open(coordRZ);
	coordinatesRZ << std::setprecision(10);

	ofstream coordinatesXYZ;
	coordinatesXYZ.open(coordXYZ);
	coordinatesXYZ << std::setprecision(10);

	ofstream totalEnergy;
	totalEnergy.open(totenergy);
	totalEnergy << std::setprecision(16);

	ofstream magMoment;
	magMoment.open(mag);
	magMoment << std::setprecision(10);

	//prepare electric potential
	setPastukhov(40, 100, mult);
	setEField();

	// initialize particle position
	Vector posi(rrlmtr_ - dr, 0, zmid_);

	// mass toggle
	assert(species == 1 || species == 0);
	Doub mass = MI * (1 - species) + ME * species;

	//initialize an hydrogen ion with energy and direction input by user;
	Doub vr = sqrt(energy * er * EVTOJOULE / mass); // thermal velocity
	Doub vphi = sqrt(energy * ephi * EVTOJOULE / mass);
	Doub vz = sqrt(energy * ez * EVTOJOULE / mass);

    Vector veli(vr, vphi, vz);
    Particle part(posi, veli, spec);

    // A default electric field of 0;
    // Vector EField(0, 0, 0);

    // calculate time step
	Doub fLamor = ( 1520 * (1 - species) + 2.8E6 * species ) * BMAGAXIS; // another logical, constants from NRL p28
	Doub TLamor = 1/fLamor;
	Doub dt = TLamor / NPERORBIT;

	Doub Binit = getB(posi).mod();
	Doub rLamor = (2.38 * species + 102 * ( 1 - species ) ) * sqrt(energy) * 0.01 / BMAGAXIS; // in m
	Doub rinit = rrlmtr_ - dr;
	Doub mirrorMid = getMirrorRatio(rinit, zmid_);
	Doub mirrorLeft = getMirrorRatio(rinit - rLamor, zmid_);
	Doub mirrorRight = getMirrorRatio(rinit + rLamor, zmid_);

	std::cerr<< "Range of R:" << mirrorLeft << "," << mirrorMid << "," << mirrorRight << std::endl;

    // Doub dt = 10E-10; // keep this number for ions.

    int step = 0;
    for (step; step < 1500000; ++step) // basically run it till it's lost
    {
    	Vector posNow = part.pos();
    	Vector BNow = getB(posNow);
    	Vector ENow = getE(posNow);

    	part.moveCyl(ENow, BNow, dt);

    	Doub xNow = part.pos().x();
    	Doub yNow = part.pos().y();
    	Doub zNow = part.pos().z();

    	Doub rNow = sqrt( xNow * xNow + yNow * yNow );
    	Doub phiNow = atan(yNow / xNow);

    	if (isLimiter(rNow, zNow)){
    		std::cerr << "particle lost to limiter after" << step \
    		<< "iterations." << std::endl;
    		break;
    	}

    	if (step % 1000 == 0){ // output every 1000 steps
    	// if (true){ // always output
	    	coordinatesRZ  << rNow << "," << phiNow << "," << zNow << std::endl;
	    	coordinatesXYZ << xNow << "," << yNow << "," << zNow << std::endl;

	    	Vector Bupdated = getB(posNow);
	    	Vector Bcart;
	    	Bupdated.cyl2Cart(part.pos(), Bcart); // convert B to cartesian vector

	    	Vector vperp = part.vel().perp(Bcart);

	    	Doub energy = mass * part.vel().dot(part.vel()) / (EVTOJOULE);
    		Doub mu = part.mu(Bcart);
	    	totalEnergy    << energy << std::endl;
	    	magMoment      << mu << ',' << vperp.mod() << ',' << Bcart.mod() << std::endl;
	    	// std::cerr << step << std::endl;
	    }
    }

    std::cerr << "program terminated after" << step << "iterations." << std::endl;

    coordinatesRZ.close();

    coordinatesXYZ.close();

    totalEnergy.close();
    magMoment.close();

    return;
}

Doub Orbit::particleStats(Doub dr, Doub energy, bool spec, int nparts, \
	Doub Ti, Doub Te, Doub mult, int maxiter, bool write)
{
	std::list<Vector> initVel;
	std::list<Vector> finlVel;

	std::list<Vector> initVel3;
	std::list<Vector> finlVel3;

	std::list<Doub> paraVel;

	Doub rinit = rrlmtr_ - dr;
	std::cerr << "mirror ratio: " << getMirrorRatio(rinit, zmid_) \
	<< std::endl;

	Doub mass = MI * (1 - spec) + ME * spec; // logical statement, choosing between ion and electron mass.
	Doub vbar = sqrt(energy  *  EVTOJOULE / mass); // thermal velocity, <v^2> in distribution, sigma.

	Doub fLamor = ( 1520 * (1 - spec) + 2.8E6 * spec ) * BMAGAXIS; // another logical, constants from NRL p28
	Doub TLamor = 1/fLamor;
	Doub dt = TLamor / NPERORBIT;


	// //prepare electric potential
	setPastukhov(Ti, Te, mult);
	setEField();

	std::default_random_engine generator(int(time(NULL)));
    std::normal_distribution<double> distribution(0.0, vbar); // generate a Gaussian distributed velocity

	#pragma omp parallel
	{
		generator.seed( int(time(NULL)) ^ omp_get_thread_num() ); // seed the distribution generator

		int er, ephi, ez, diff;
		Doub vr, vphi, vz, xNow, yNow, zNow, rNow, phiNow;
		Vector veli, posi, posNow, BNow, ENow, vPara, vPerp, vGC;
		Particle part;

		part.setSpec(spec);

		// A default electric field of 0;
	    // Vector EField(0, 0, 0);

	    std::list<Vector> initVel_private; // A list of (vpara, vperp)
	    std::list<Vector> finlVel_private;

	    std::list<Vector> initVel3_private; // A list of 3D velocity (vr, vphi, vz)
	    std::list<Vector> finlVel3_private;


	    std::list<Doub>   paraVel_private; // A list of parallel velocity at exit

		#pragma omp for private(part, er, ephi, ez, vr, vphi, vz, veli, \
		posi, xNow, yNow, zNow, rNow, phiNow, diff, posNow, BNow, vGC)		
			for (int i=0; i < nparts; ++i ){
				// initialize particle position
				posi = Vector(rrlmtr_ - dr, 0, zmid_);

				//initialize an hydrogen ion with energy input by user, in random directions;
				vr = distribution(generator); // generate 3 normal distributed velocities.
				vphi =  distribution(generator);
				vz = distribution(generator);

			    veli = Vector(vr, vphi, vz);

			    BNow = getB(posi);
	    	    ENow = getE(posi);

			    vPara = veli.parallel(BNow);
			    vPerp = veli.perp(BNow);
			    vGC = Vector(vPara.mod(), vPerp.mod(), 0); // Guiding Center velocity in (vpara, vperp)

			    initVel_private.push_back(vGC);
			    initVel3_private.push_back(veli); 

			    part.setPos(posi);
			    part.setVel(veli);
			    part.setSpec(spec);
			    for (int step = 0; step < maxiter; ++step){ 
					posNow = part.pos();
			    	BNow = getB(posNow);
	    			ENow = getE(posNow);

			    	part.moveCyl(ENow, BNow, dt);
			    	// part.moveCyl(EField, BNow, dt);

			    	xNow = part.pos().x();
			    	yNow = part.pos().y();
			    	zNow = part.pos().z();
			    	rNow = sqrt( xNow * xNow + yNow * yNow );
			    	phiNow = atan(yNow / xNow);

			    	if (isLimiter(rNow, zNow)){
			    		part.lost();
					    finlVel_private.push_back(vGC); // a list of initial velocities that are lost
			    		finlVel3_private.push_back(veli);
			    		// collect the final paralell velocity here
					    vPara = part.vel().parallel(BNow);
					    paraVel_private.push_back(vPara.mod());
			    		break;
			    	}
			    }
			}
		#pragma omp critical
			// collect everything from all threads back into the main structure
			initVel.insert(initVel.end(), initVel_private.begin(), initVel_private.end());
			finlVel.insert(finlVel.end(), finlVel_private.begin(), finlVel_private.end());

			initVel3.insert(initVel3.end(), initVel3_private.begin(), initVel3_private.end());
			finlVel3.insert(finlVel3.end(), finlVel3_private.begin(), finlVel3_private.end());

			paraVel.insert(paraVel.end(), paraVel_private.begin(), paraVel_private.end());

	}

	std::cerr << "initial velocities: " << initVel.size() << std::endl;	
	std::cerr << std::endl << "ones that were lost: " << finlVel.size() << std::endl;

	int species = spec;
	std::string suffix = "_Ti_" + std::to_string(Ti).substr(0, 3) \
	+ "_Te_" + std::to_string(Te).substr(0, 3) \
	+ "_dr_" + std::to_string(dr).substr(0, 4) \
	+ "_mult_" + std::to_string(mult).substr(0, 3) \
	+ "_spec_" + std::to_string(species).substr(0, 1) \
	+ ".out";

	if (write){
		ofstream initial;
		initial.open("output/initial" + suffix);
		ofstream final;
		final.open("output/final" + suffix);
		initial << Ti << Te << std::endl;

		ofstream initial3;
		initial3.open("output/initial3" + suffix);
		ofstream final3;
		final3.open("output/final3" + suffix);

		while (!initVel.empty()){
			initial << initVel.front() << std::endl;
			initial3 << initVel3.front() << std::endl;

			initVel.pop_front();
			initVel3.pop_front();
		}	

		while (!finlVel.empty()){
			final << finlVel.front() << std::endl;
			final3 << finlVel3.front() << std::endl;

			finlVel.pop_front();
			finlVel3.pop_front();
		}
	} 

	Doub gammaOut = 0;
	while(!paraVel.empty()){
		gammaOut += paraVel.front();
		paraVel.pop_front();
	}	
	return gammaOut;
}


