// Programmer class includes
#include "Orbit.h"
#include <omp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>
#include <numeric> //std::accumulate
#include <algorithm>
#include <iterator>

#define BMAGAXIS  2000                    // mod(B) = 0.2 T = 2000 Gauss at magnetic axis, characteristic field strength
#define NPERORBIT 20                      // steps per Lamor orbit




//-----------------------------------------------------------------
//------------------------- Controllers  --------------------------
//-----------------------------------------------------------------

void printData(Orbit& orbit)
{
	/*
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

	/*
	MatDoub pass(orbit.nw_, orbit.nh_);

	*/
	ofstream passout1; // LOLOLOLOL
	passout1.open("./output/passTest.out");
	// input Ti, Te here. keep Te 1.
	orbit.writePastukhov(0.2, 1, passout1); // pass modified
	passout1.close();

	/*
	ofstream midplane1;
	midplane1.open("./output/midplane02.out");
	orbit.midPlanePhi(0.2, 1, pass, midplane1);
	midplane1.close();

	// MatDoub pass(orbit.nw_, orbit.nh_);

	ofstream passout2; // LOLOLOLOL
	passout2.open("./output/pass1.out");
	// input Ti, Te here. keep Te 1.
	orbit.writePastukhov(1, 1, pass, passout2); // pass modified
	passout2.close();

	ofstream midplane2;
	midplane2.open("./output/midplane1.out");
	orbit.midPlanePhi(1, 1, pass, midplane2);
	midplane2.close();

	// MatDoub pass(orbit.nw_, orbit.nh_);
	ofstream passout3; // LOLOLOLOL
	passout3.open("./output/pass2.out");
	// input Ti, Te here. keep Te 1.
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
	*/
}

/** TODO: make this a member function */
void temperature(Doub Ti_start, Doub dT, int iter, Doub R, Orbit& orbit)
{
	ofstream temp;
	temp.open("./output/phivsT.out");
	temp << std::setprecision(10);

	int i(0);
	Doub TiNow = Ti_start;
	while (i <= iter){
		Doub phi   = orbit.pastukhov(TiNow, 1, R);
		temp << TiNow << ',' << phi << std::endl;
		i++;
		TiNow += dT;
	}
	temp.close();

	return;
}

/** TODO: make this a member function */
void particlePush(Orbit& orbit, Doub dr, Doub energy, Doub er, Doub ephi, Doub ez)
{
	std::string prefix  = "./output/";
	std::string suffix = "_E_" + std::to_string(energy) + "_dr_" + std::to_string(dr) + ".out";

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

	// initialize particle position
	Vector posi(orbit.rrlmtr_ - dr, 0, orbit.zmid_);

	//initialize an hydrogen ion with energy and direction input by user;
	Doub vr = sqrt(energy * er * EVTOJOULE / MI); // thermal velocity
	Doub vphi = sqrt(energy * ephi * EVTOJOULE / MI);
	Doub vz = sqrt(energy * ez * EVTOJOULE / MI);

    Vector veli(vr, vphi, vz);
    Particle part(posi, veli, 0);

    // A default electric field of 0;
    Vector EField(0, 0, 0);

    Doub dt = 10E-10; // keep this number for ions.

    int step = 0;
    for (step; step < 1000000; ++step) // basically run it till it's lost
    {
    	Vector posNow = part.pos();
    	Vector BNow = orbit.getB(posNow);

    	part.moveCyl(EField, BNow, dt);

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

/** TODO: make this a member function. No reason why this isn't a member function */
/**
 * \brief         Push particles in paralell at given midplane position, with Gaussian distributed 
 *                initial velocities
 * \return        Sum of all velocities for all particles lost to the limiter
 * \param orbit   Input Orbit object, carries magnetic geometry and more
 * \param dr      Radial location for particle loading
 * \param energy  Temperature of test particles (in eV), to initialize Maxwellian dist.
 * \param spec    Species of the particle. 0 for H, 1 for e.
 * \param nparts  Number of particles
 * \param maxiter Maximum number of iterations for each particle.
 * \param write   Whether to write list of initial and lost velocities to file, default to false.
 */
Doub particleStats(Orbit& orbit, Doub dr, Doub energy, bool spec, int nparts, int maxiter, bool write = false)
{
	std::list<Vector> initVel;
	std::list<Vector> finlVel;

	std::list<Doub> paraVel;

	std::cerr<< "mirror ratio: " << orbit.getMirrorRatio(orbit.rrlmtr_ - dr, orbit.zmid_) << std::endl;

	Doub mass = MI * (1 - spec) + ME * spec; // logical statement, choosing between ion and electron mass.
	Doub vbar = sqrt(energy  *  EVTOJOULE / mass); // thermal velocity, <v^2> in distribution, sigma.

	Doub fLamor = ( 1520 * (1 - spec) + 2.8E6 * spec ) * BMAGAXIS; // another logical, constants from NRL p28
	Doub TLamor = 1/fLamor;
	Doub dt = TLamor / NPERORBIT;

	std::default_random_engine generator(int(time(NULL)));

    std::normal_distribution<double> distribution(0.0, vbar); // generate a Gaussian distributed velocity

	#pragma omp parallel
	{
		// srand( int(time(NULL)) ^ omp_get_thread_num() ); // seed each thread independently
		generator.seed( int(time(NULL)) ^ omp_get_thread_num() ); // seed the distribution generator
		// std::cerr<< "num of thread:" << omp_get_thread_num() << std::endl;

		int er, ephi, ez, diff;
		Doub vr, vphi, vz, xNow, yNow, zNow, rNow, phiNow;
		Vector veli, posi, posNow, BNow, vPara, vPerp, vGC;
		Particle part;

		part.setSpec(spec);

		// A default electric field of 0;
	    Vector EField(0, 0, 0);

	    std::list<Vector> initVel_private; // A list of (vpara, vperp)
	    std::list<Vector> finlVel_private;
	    std::list<Doub>   paraVel_private; // A list of parallel velocity at exit

	    // std::cerr << "number of threads: " << omp_get_num_threads() << std::endl;

		#pragma omp for private(part, er, ephi, ez, vr, vphi, vz, veli, posi, xNow, yNow, zNow, rNow, phiNow, diff, posNow, BNow, vGC)
		
			for (int i=0; i < nparts; ++i ){

				// initialize particle position
				posi = Vector(orbit.rrlmtr_ - dr, 0, orbit.zmid_);

				//initialize an hydrogen ion with energy input by user, in random directions;
				vr = distribution(generator); // generate 3 normal distributed velocities.
				vphi =  distribution(generator);
				vz = distribution(generator);

			    veli = Vector(vr, vphi, vz);

			    BNow = orbit.getB(posi);
			    vPara = veli.parallel(BNow);
			    vPerp = veli.perp(BNow);
			    vGC = Vector(vPara.mod(), vPerp.mod(), 0); // Guiding Center velocity in (vpara, vperp)

			    initVel_private.push_back(vGC); 

			    part.setPos(posi);
			    part.setVel(veli);

			    for (int step = 0; step < maxiter; ++step){ 
					posNow = part.pos();
			    	BNow = orbit.getB(posNow);

			    	// std::cerr << "BField:" << BNow << std::endl;
			    	part.moveCyl(EField, BNow, dt);

			    	xNow = part.pos().x();
			    	yNow = part.pos().y();
			    	zNow = part.pos().z();

			    	rNow = sqrt( xNow * xNow + yNow * yNow );
			    	phiNow = atan(yNow / xNow);

			    	if (orbit.isLimiter(rNow, zNow)){
			    		part.lost();
					    finlVel_private.push_back(vGC); // a list of initial velocities that are lost
			    		// std::cerr << "particle lost to limiter after " << step << " iterations." << std::endl;
			    		
			    		// collect the final paralell velocity here
					    vPara = part.vel().parallel(BNow);
					    paraVel_private.push_back(vPara.mod());
			    		break;
			    	}
			    }
			    /** Apparently as soon as you get out of the for loop things are scrambled. This push_back has
			        to be inside the serial for loop */
			 //    if (part.isLost()){
				//     finlVel_private.push_back(veli); // a list of initial velocities that are not lost
				// }


			}
		#pragma omp critical
			// collect everything from all threads back into the main structure
			initVel.insert(initVel.end(), initVel_private.begin(), initVel_private.end());
			finlVel.insert(finlVel.end(), finlVel_private.begin(), finlVel_private.end());
			paraVel.insert(paraVel.end(), paraVel_private.begin(), paraVel_private.end());

	}

	std::cerr << "initial velocities: " << initVel.size() << std::endl;	
	std::cerr << std::endl << "ones that were lost: " << finlVel.size() << std::endl;


	std::string suffix = "_E_" + std::to_string(energy).substr(0, 5) + "_dr_" + std::to_string(dr).substr(0, 4) + ".out";

	// Doub sumVel = 
	if (write){
		ofstream initial;
		initial.open("output/initial" + suffix);
		ofstream final;
		final.open("output/final" + suffix);
		while (!initVel.empty()){
			initial << initVel.front() << std::endl;
			initVel.pop_front();
		}	

		while (!finlVel.empty()){
			final << finlVel.front() << std::endl;
			finlVel.pop_front();
		}
	} 

	Doub gammaOut = 0;
	while(!paraVel.empty()){
		gammaOut += paraVel.front();
		paraVel.pop_front();
	}
	
	return gammaOut;
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

//-----------------------------------------------------------------
//------------ Initialization and User Interface  -----------------
//-----------------------------------------------------------------

int main(int argc, const char** argv)
{
	std::cerr << "Hello World" << std::endl << std::endl;
	// std::cout << std::setprecision(10);
	std::cout << std::scientific;

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
    	// temperature(0.1, 0.01, 200, 3,orbit);

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
    	bool spec;
    	int nparts, maxiter;
    	bool write;
    	// TODO implement this
    	if (args.size() == 6){
    		dr = stod(args.front());
    		args.pop_front();
    		energy = stod(args.front());
    		args.pop_front();
    		spec = stoi(args.front());
    		args.pop_front();
    		nparts = stod(args.front());
    		args.pop_front();
    		maxiter = stod(args.front());
    		args.pop_front();
    		write = stoi(args.front());
    		args.pop_front();
    	} else {
    		std::cerr << "Input argument list incomplete. Let's try again:" << std::endl;
    		std::cerr << "Give a starting radial position, dr from limiter:" << std::endl;
    		std::cin >> dr;
    		std::cerr << "Energy of the particle:" << std::endl;
    		std::cin >> energy;
    		std::cerr << "Species of the particle, 1 for electron, 0 for hydrogen:" << std::endl;
    		std::cin >> spec;
    		std::cerr << "Number of particles to push: " << std::endl;
    		std::cin >> nparts;
    		std::cerr << "Maximum iteration for each particle: " << std::endl;
    		std::cin >> maxiter;
    		std::cerr << "Write loss cone to file? 1 for yes, 0 for no." << std::endl;
    		std::cin >> write;
    		std::cerr << "All set." << std::endl;
    	}
    	std::cerr << "Particle successfully initialized. Pushing now." << std::endl;

    	Doub gammaOut = particleStats(orbit, dr, energy, spec, nparts, maxiter, write);

    	std::cout << '(' << nparts << ',' << gammaOut << ")," << std::endl;

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



