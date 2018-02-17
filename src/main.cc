// Programmer class includes
#include "Orbit.h"


//-----------------------------------------------------------------
//------------------------- Controllers  --------------------------
//-----------------------------------------------------------------



void printData(Orbit& orbit)
{
	/*
	Vector Bstart( orbit.rleft_+ orbit.rdim_/2 + 0.1, 0, orbit.zmid_  - orbit.zdim_/2 + 0.2);
	Vector dv(0.05, 0, 0);

	
	ofstream rBList;
    rBList.open ("rBList.bout");

    ofstream zBList;
    zBList.open ("zBList.bout");

    ofstream phiBList;
    phiBList.open("phiBList.bout");

	for (int i = 0; i < 8; i++){
		Vector Bnow = Bstart + dv * i;
		// std::cerr << Bnow << std::endl;
	    orbit.fieldLines(rBList, zBList, phiBList, Bnow, 0.01, 1000);
	    // orbit.fieldLines(rBList, zBList, phiBList, Bnow, -0.01, 1000);
	}

	rBList.close();
	zBList.close();
	phiBList.close();
	*/

	ofstream modB;
	modB.open("./output/modB.out");

	orbit.writeModB(modB);

	modB.close();

	ofstream psiout;
	psiout.open("./output/flux.out");
	orbit.writeFlux(psiout);
	psiout.close();

	ofstream RRout;
	RRout.open("./output/mirrorRatio.out");
	orbit.writeMirrorRatio(RRout);
	RRout.close();

	MatDoub pass(orbit.nw_, orbit.nh_);

	ofstream passout1; // LOLOLOLOL
	passout1.open("./output/pass02.out");
	// input Ti, Te here. keep Te 0.
	orbit.writePastukhov(0.2, 1, pass, passout1); // pass modified
	passout1.close();

	ofstream midplane1;
	midplane1.open("./output/midplane02.out");
	orbit.midPlanePhi(0.2, 1, pass, midplane1);
	midplane1.close();

	// MatDoub pass(orbit.nw_, orbit.nh_);

	ofstream passout2; // LOLOLOLOL
	passout2.open("./output/pass1.out");
	// input Ti, Te here. keep Te 0.
	orbit.writePastukhov(1, 1, pass, passout2); // pass modified
	passout2.close();

	ofstream midplane2;
	midplane2.open("./output/midplane1.out");
	orbit.midPlanePhi(1, 1, pass, midplane2);
	midplane2.close();

	// MatDoub pass(orbit.nw_, orbit.nh_);
	ofstream passout3; // LOLOLOLOL
	passout3.open("./output/pass2.out");
	// input Ti, Te here. keep Te 0.
	orbit.writePastukhov(2, 1, pass, passout3);
	passout3.close();

	ofstream midplane3;
	midplane3.open("./output/midplane2.out");
	orbit.midPlanePhi(2, 1, pass, midplane3);
	midplane3.close();


	ofstream eFieldLine2;
	eFieldLine2.open("./output/eField2.out");

	ofstream eFieldLine1;
	eFieldLine1.open("./output/eField1.out");

	ofstream eFieldLine02;
	eFieldLine02.open("./output/eField02.out");

	Doub rStart = orbit.rllmtr_ + 0.2;
	Doub zStart = abs(orbit.getzLimiter(rStart)) - 0.01;

	// Vector start( 0.575, 0, 0);
	Vector start(rStart, 0, zStart);

	orbit.eField(2, 1, start, 0.005, 15000, eFieldLine2);
	orbit.eField(1, 1, start, 0.005, 15000, eFieldLine1);
	orbit.eField(0.2, 1, start, 0.005, 15000, eFieldLine02);
}

void temperature(Doub Ti_start, Doub dT, int iter, Doub R, Orbit& orbit)
{
	ofstream temp;
	temp.open("./output/phivsT.out");
	temp << std::setprecision(10);

	int i(0);
	Doub TiNow = Ti_start;
	while (i <= iter){
		Doub phi   = orbit.pastukhov(TiNow, 1, R);

		// tList[i]   = TiNow;
		// phiList[i] = phi;

		temp << TiNow << ',' << phi << std::endl;

		i++;
		TiNow += dT;
	}
	temp.close();

	return;
}

void particlePush(Orbit& orbit, Doub dr, Doub energy, Doub er, Doub ephi, Doub ez)
{

	// TODO: Dynamically create file names.
	std::string prefix  = "./output/";
	std:: string suffix = "_E_" + std::to_string(energy) + "_dr_" + std::to_string(dr) + ".out";

	std::string coordRZ = prefix + "coordRZ" + suffix;
	std::string coordXYZ = prefix + "coordXYZ" + suffix;
	std::string totenergy   = prefix + "totalEnergy" + suffix;
	std::string mag      = prefix + "Mu" + suffix;


	ofstream coordinatesRZ;
	// coordinatesRZ.open("./output/coordinatesRZ.out");	
	coordinatesRZ.open(coordRZ);
	coordinatesRZ << std::setprecision(10);

	ofstream coordinatesXYZ;
	// coordinatesXYZ.open("./output/coordinatesXYZ.out");	
	coordinatesXYZ.open(coordXYZ);
	coordinatesXYZ << std::setprecision(10);

	ofstream totalEnergy;
	// totalEnergy.open("./output/totalEnergy.out");
	totalEnergy.open(totenergy);
	totalEnergy << std::setprecision(16);

	ofstream magMoment;
	// magMoment.open("./output/magMoment.out");
	magMoment.open(mag);
	magMoment << std::setprecision(10);

	// initialize particle position
	Vector posi(orbit.rrlmtr_ - dr, 0, orbit.zmid_);

	//initialize an hydrogen ion with energy and direction input by user;

	Doub vr = sqrt(energy * er * EVTOJOULE / MI); // thermal velocity
	Doub vphi = sqrt(energy * ephi * EVTOJOULE / MI);
	Doub vz = sqrt(energy * ez * EVTOJOULE / MI);

	// std::cerr << vv << std::endl;
    Vector veli(vr, vphi, vz);


    Particle part(posi, veli, 0);

    // ofstream init;
    // init.open(prefix + "init" + suffix);
    // init << veli << std::endl << posi << std::endl;
    // init.close();

    // A default electric field of 0;
    Vector EField(0, 0, 0);

    Doub dt = 10E-10; // keep this number for ions.

    int step = 0;
    for (step; step < 1000000; ++step) // basically run it till it's lost
    {
    	Vector posNow = part.pos();

    	// if (part.pos().x() < 10E-6){
    	// 	break;
    	// }
    	Vector BNow = orbit.getB(posNow);

    	// std::cerr << "BField:" << BNow << std::endl;
    	part.moveCyl(EField, BNow, dt);
    	// part.move(EField, BNow, dt);

    	// Doub rNow = part.pos().x();
    	// Doub phiNow = part.pos().y();
    	// Doub zNow = part.pos().z();

    	Doub xNow = part.pos().x();
    	Doub yNow = part.pos().y();
    	Doub zNow = part.pos().z();

    	Doub rNow = sqrt( xNow * xNow + yNow * yNow );
    	Doub phiNow = atan(yNow / xNow);

    	if (orbit.isLimiter(rNow, zNow)){
    		std::cerr << "particle lost to limiter after" << step << "iterations." << std::endl;
    		break;
    	}

    	if (step % 500 == 0){ // output every 500 steps
	    	// Vector vperp = part.vel()
	    	coordinatesRZ  << rNow << "," << phiNow << "," << zNow << std::endl;
	    	coordinatesXYZ << xNow << "," << yNow << "," << zNow << std::endl;

	    	Doub energy = MI * part.vel().dot(part.vel()) / (EVTOJOULE);
    		Doub mu = part.mu(BNow);
	    	totalEnergy    << energy << std::endl;
	    	magMoment      << mu << std::endl;
	    }
    }

    std::cerr << "program terminated after" << step << "iterations." << std::endl;

    coordinatesRZ.close();
    coordinatesXYZ.close();
    totalEnergy.close();
    magMoment.close();

    return;
}

void particleStats(Orbit& orbit, Doub dr, Doub energy, int nparts, int maxiter)
{
	return;
}


//-----------------------------------------------------------------
//------------------------- Housekeeping --------------------------
//-----------------------------------------------------------------

/// Allows for exiting without memory leaks
struct ExitException {
    int code_; ///< Exit code
    ExitException(int code): code_(code) {} ///< Trivial constructor
};

/**
 * \brief Prints how to use the program from command line then exits.
 * \throws ExitException in order to exit the program.
 */
void usage()
{
    std::cerr << "Usage: ./Orbit [-t] [-e] [-g] [-h] [-dat] [-part] [-stat]" << std::endl
                  << "See documentation or use [-h] option for details :)" <<  std::endl << "Exiting now." << std::endl;
    throw ExitException(2);
}


/**
 * \brief Prints detailed help message on the usage of the program.
 *
 */
void help()
{
    //write help message here.
    std::cout << "Command line options:" << std::endl << std::endl;
    std::cout << "[-t]   Run unit testing" << std::endl;
    std::cout << "[-e]   Test particle orbit in a purely toroidal field" << std::endl;
    std::cout << "[-g]   Load field from g-eqdsk file (not in use anymore)" << std::endl;
    std::cout << "[-dat] Load field from .dat file, output data from pastukhov potential calculations" << std::endl;
	std::cout << "[-h]   Print help message." << std::endl << std::endl;

	std::cout << "[-part] [dr] [energy] [E_r] [E_phi] [E_z] Push particles in LTX geometry: " << std::endl;
	std::cout << "        [dr]      < m >  radial distance from limiter at midplane" << std::endl;
	std::cout << "        [energy]  < eV > energy of the particle" << std::endl;
	std::cout << "        [E_r]     < - >  fraction of energy in r direction" << std::endl;
	std::cout << "        [E_phi]   < - >  fraction of energy in phi direction" << std::endl;
	std::cout << "        [E_z]     < - >  fraction of energy in z direction" << std::endl << std::endl;

	std::cout << "[-stat] [dr] [energy] [n] Push a statistical sample of particles" << std::endl;
	std::cout << "        [dr]      < m >  radial distance from limiter at midplane" << std::endl;
	std::cout << "        [energy]  < eV > energy of the particle" << std::endl;
	std::cout << "        [n]       < - >  number of particles for the statistical sample" << std::endl << std::endl;

	std::cout << "Output Files:" << std::endl << std::endl;
	std::cout << "'coordRZ_E_[energy]_dr_[dr].out'     > Coordinate of the test particle in R_phi_Z" << std::endl;
	std::cout << "'coordXYZ_E_[energy]_dr_[dr].out'    > Coordinate of the test particle in X_Y_Z" << std::endl;
	std::cout << "'totalEnergy_E_[energy]_dr_[dr].out' > Total kinetic energy of the particle" << std::endl;
	std::cout << "'Mu_E_[energy]_dr_[dr].out'          > Magnetic moment of the particle" << std::endl;

	// TODO: actually implement these options inputs
	return;
}


int main(int argc, const char** argv)
{
	std::cerr << "Hello World" << std::endl << std::endl;;
	std::cout << std::setprecision(10);

	std::cerr << "Magnetic geometry initialization: " << std::endl;
	// read from uniform grid data file.
	// Hard code in path for now
	std::string field = "./input/LTX_Apr29_474-fields.dat";
	std::cerr << "magnetic field path read < " << field << std::endl;

	std::string limiter = "./input/pos_coord_invert.csv";
	std::cerr << "limiter path read < " << limiter << std::endl;

	Orbit orbit(field, limiter);

	std::cerr << "Initialization successful." << std::endl << std::endl;

	// // Process command-line options
    // --argc; // Skip past 0th arg (program name)
    // ++argv;

    std::list<std::string> options(argv + 1, argv + argc);
    std::list<std::string> args;

    if (options.empty()){
        usage();
    }

    std::string controller;
    // Process the options and input arguments
    while (!options.empty()){

    	std::string option = options.front();
        options.pop_front();

        if (option[0] == '-'){
        	controller = option;
        }
        else //if (option.empty() || option[0] != '-') 
        {
            args.push_back(option); // all input arguments are stored in this list as std::string
        }
    } 

    // execute the options

    if (controller == std::string("-t")) {
    	// run unit testing
    	// Orbit orbit;
        orbit.test();
    } 
    else if (controller == std::string("-e")) {
    	Orbit emptyOrbit;
    	emptyOrbit.emptytest();
    }
    else if (controller == std::string("-g"))
    {
    	// read from eqdsk g file. Hard code in path for now
    	std::string path = "./input/LTX_1504291255_47211.eqdsk";
		std::cerr << "path read" << path << std::endl;

		std::cerr << "eqdsk no longer supported. Include gConstructor.cc for eqdsk file read." << std::endl;
		// Orbit orbit(path, 1);
    }
    else if (controller == std::string("-dat"))
    {
    	// TODO implement command line input to Ti Te, and starting point on limiter
    	printData(orbit);
    	temperature(0.1, 0.01, 200, 3,orbit);

    }
    else if (controller == std::string("-part"))
    {
    	// TODO implement this
    	std::cerr << "Let's push some particles!" << std::endl;
    	double dr, energy, er(0), ephi(0), ez(0);

    	if (args.size() != 5){
    		std::cerr << "Input argument list incomplete. Let's try again:" << std::endl;
    		std::cerr << "Give a starting radial position, dr from limiter:" << std::endl;
    		std::cin >> dr;
    		std::cerr << "Energy of the particle:" << std::endl;
    		std::cin >> energy;
    		

    		while ( abs(er + ephi + ez - 1) > 1E-16){
    			std::cerr << "Energy fractions need to add up to 1! " << std::endl;
    			std::cerr << "Fraction of energy in r direction:" << std::endl;
	    		std::cin >> er;
	    		std::cerr << "Fraction of energy in phi direction:" << std::endl;
	    		std::cin >> ephi;
	    		std::cerr << "Fraction of energy in z direction:" << std::endl;
	    		std::cin >> ez;
    		}
    		std::cerr << "All set." << std::endl;

    	} else {
    		dr = stod(args.front());
    		args.pop_front();

    		energy = stod(args.front());
    		args.pop_front();

    		er = stod(args.front());
    		args.pop_front();

    		ephi = stod(args.front());
    		args.pop_front();

    		ez = stod(args.front());
    		args.pop_front();
    	}
    	std::cerr << "Particle successfully initialized. Pushing now." << std::endl;

    	//dr = 0.06 for a big SOL banana!
    	particlePush(orbit, dr, energy, er, ephi, ez);

    	std::cerr << "Done." << std::endl;
    }
    else if (controller == std::string("-stat"))
    {
    	double dr, energy;
    	int nparts, maxiter;
    	// TODO implement this
    	if (args.size() == 4){
    		dr = stod(args.front());
    		args.pop_front();
    		energy = stod(args.front());
    		args.pop_front();
    		nparts = stod(args.front());
    		args.pop_front();
    		maxiter = stod(args.front());
    		args.pop_front();
    	} else {
    		std::cerr << "Input argument list incomplete. Let's try again:" << std::endl;
    		std::cerr << "Give a starting radial position, dr from limiter:" << std::endl;
    		std::cin >> dr;
    		std::cerr << "Energy of the particle:" << std::endl;
    		std::cin >> energy;
    		std::cerr << "Number of particles to push: " << std::endl;
    		std::cin >> nparts;
    		std::cerr << "Maximum iteration for each particle: " << std::endl;
    		std::cin >> maxiter;
    		std::cerr << "All set." << std::endl;
    	}
    	std::cerr << "Particle successfully initialized. Pushing now." << std::endl;

    	particleStats(orbit, dr, energy, nparts, maxiter);

    }
    else if (controller == std::string("-h")){
    	help();
    }
    else 
    {
    	std::cerr << "Unrecognized option" << std::endl;
    	usage();
    }
    
	return 0;
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
			// std::cout << xx[i] << ',' << yy[i] << std::endl;
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
		Vector vv(4, 3, 0);
		Vector B(10, 0, 0);

		Vector vpara(4, 0, 0);
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

    // Vector veli(20000, 0, 0);
    // // Vector veli(2E5, 0, 0);

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

