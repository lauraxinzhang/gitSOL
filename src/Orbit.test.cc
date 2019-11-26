/**
 * \file    Orbit.test.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implementations for unit testing suite of the Orbit class.
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
	bool turning          = true;

	// Particle Tests:
	bool partInitAndAssign = true;
	bool mover             = true;
	bool testMu                = true;

	bool testPastukhov = false;
	bool testEField    = false;

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

		Vector v(2, 3, 5); // Constructor
		Vector v2 = v;     // assignment operator
		Vector v3;         // default constructor
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
		Vector vv(-4, 3, -1);
		Vector B(10, 0, 0);

		Vector vpara(-4, 0, 0);
		Vector vperp(0, 3, -1);

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

	if (testVector && turning){
		std::cerr << " - Testing vector turning" << std::endl;
		Vector axis(1, 0, 0);

		Vector vi(0, 4, 2); // a test vector in x-y plane;
		Vector vplus(0, -2, 4); // a test vector in x-z plane;
		Vector vminus(0, 2, -4);

		bool result = (vi.turn(axis, 1) == vplus) && (vi.turn(axis, 0) == vminus);
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
			std::cerr << "positive turn:" << vi.turn(axis, 1) << "should be: " << vplus<< std::endl;
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
		bool perp = (muCorrect == part.mu(B));

		Vector veldiag(1, 1, 0);
		part.setVel(veldiag);
		bool diag = (part.mu(B) == muCorrect); // adding a parallel shouldn't change mu

		Vector velpara(1, 0, 0);
		part.setVel(velpara);
		bool para = (part.mu(B) == 0);

		bool result = (perp && diag && para);
		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}

	}

	if (testPastukhov) {
		std::cerr << " - Testing empty pastukhov calculation" << std::endl;

		Vector pos(rrlmtr_ - 0.06, 0, zmid_);

		setPastukhov(0.2, 1, 0); // create a uniformly zero potential
		Vector zero(0, 0, 0);

		bool notnull = (Phi_ != nullptr);
		bool blank = ((*Phi_)[230][130] == 0 || (*Phi_)[230][130] == NAN);

		bool result = notnull && blank;

		if (result){
			std::cerr << " -- Passed" << std::endl;
		} else {
			std::cerr << "Something's wrong" << std::endl;
		}
	}

	if (testEField) {
		std::cerr << " - Testing electric field calculation" << std::endl;

		setPastukhov(0.2, 1, 0); // create a uniformly zero potential
		setEField();
		Vector pos(rrlmtr_ - 0.06, 0, zmid_);
		Vector eTest = getE(pos);

		bool result = (eTest == Vector(0, 0, 0));
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
