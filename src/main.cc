// Programmer class includes
#include "Orbit.h"
#include "Mirror.h"
#include "Pusher.h"

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
	// std::cerr << "Hello World" << std::endl << std::endl;
	// std::cout << std::setprecision(10);
	std::cout << std::scientific;

	// std::cerr << "Magnetic geometry initialization: " << std::endl;

	// read from uniform grid data file.
	// Hard code in path for now
	std::string field = "./input/LTX_Apr29_474-fields.dat";
	// std::cerr << "magnetic field path read < " << field << std::endl;

	std::string limiter = "./input/pos_coord_invert.csv";
	// std::cerr << "limiter path read < " << limiter << std::endl;

	Orbit orbit(field, limiter);

	// std::cerr << "Initialization successful." << std::endl << std::endl;

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
    	// orbit.printData();
    	// temperature(0.1, 0.01, 200, 3,orbit);
        orbit.writePassing(1, 0.2, 1.1, 0.2, 5, 19);
    }
    else if (controller == std::string("-part"))
    {
    	// TODO implement this
    	std::cerr << "Let's push some particles!" << std::endl;
    	double dr, energy, mult, er(0), ephi(0), ez(0);
    	bool spec;

    	if (args.size() != 7){
    		std::cerr << "Input argument list incomplete. Let's try again:" << std::endl;
    		std::cerr << "Give a starting radial position, dr from limiter:" << std::endl;
    		std::cin >> dr;
    		std::cerr << "Energy of the particle:" << std::endl;
    		std::cin >> energy;
    		std::cerr << "species of the particle, 1 for electron, 0 for hydrogen:" << std::endl;
    		std::cin >> spec;
    		

    		while ( abs(er + ephi + ez - 1) > 1E-16){
    			std::cerr << "Energy fractions need to add up to 1! " << std::endl;
    			std::cerr << "Fraction of energy in r direction:" << std::endl;
	    		std::cin >> er;
	    		std::cerr << "Fraction of energy in phi direction:" << std::endl;
	    		std::cin >> ephi;
	    		std::cerr << "Fraction of energy in z direction:" << std::endl;
	    		std::cin >> ez;
    		}
    		std::cerr << "Pastukhov multiplier:" << std::endl;
    		std::cin >> mult;

    		std::cerr << "All set." << std::endl;

    	} else {
    		dr = stod(args.front());
    		args.pop_front();

    		energy = stod(args.front());
    		args.pop_front();

    		spec = stoi(args.front());
    		args.pop_front();

    		er = stod(args.front());
    		args.pop_front();

    		ephi = stod(args.front());
    		args.pop_front();

    		ez = stod(args.front());
    		args.pop_front();

    		mult = stod(args.front());
    		args.pop_front();
    	}
    	std::cerr << "Particle successfully initialized. Pushing now." << std::endl;

    	//dr = 0.06 for a big SOL banana!
    	orbit.particlePush(dr, energy, spec, er, ephi, ez, mult);

    	std::cerr << "Done." << std::endl;
    }
    else if (controller == std::string("-stat"))
    {
    	double dr, energy, mult, Ti, Te;
    	bool spec;
    	int nparts, maxiter;
    	bool write;
    	// TODO implement this
    	if (args.size() == 9){
    		dr = stod(args.front());
    		args.pop_front();

    		energy = stod(args.front());
    		args.pop_front();

    		spec = stoi(args.front());
    		args.pop_front();

    		nparts = stod(args.front());
    		args.pop_front();

            Ti = stod(args.front());
            args.pop_front();

            Te = stod(args.front());
            args.pop_front();

    		mult = stod(args.front());
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

            std::cerr << "Ti: " << std::endl;
            std::cin >> Ti;
            std::cerr << "Te: " << std::endl;
            std::cin >> Te;
    		std::cerr << "Multiplier on potential: " << std::endl;
    		std::cin >> mult;

    		std::cerr << "Maximum iteration for each particle: " << std::endl;
    		std::cin >> maxiter;
    		std::cerr << "Write loss cone to file? 1 for yes, 0 for no." << std::endl;
    		std::cin >> write;
    		std::cerr << "All set." << std::endl;
    	}
    	// std::cerr << "Particle successfully initialized. Pushing now." << std::endl;

    	Doub gammaOut = orbit.particleStats(dr, energy, spec, nparts, Ti, Te, mult, maxiter, write);

    	// std::cout << '(' << spec << ',' << mult << ',' << nparts << ',' << gammaOut << ")," << std::endl;
        std::cout << spec << ',' << mult << Ti << ',' << Te << ',' << ',' << nparts << ',' << gammaOut << std::endl;

    }
    else if (controller == std::string("-straight")){
        Mirror mirror(0.8, 1, 0.7, 4, 2, 101); // setting up a straight box
        
        mirror.printData(getB, std::cerr);
        mirror.printData(getPhi, std::cerr);

        Pusher<Mirror> pusher(mirror); // construct a Pusher object

        //pusher.gridBurst(1.8, 0.116, 5000, 1);
        //pusher.conicBurst(1.8, 0.116, 0, 5000, 8, 1);
    }
    else if (controller == std::string("-passing")) {
        std::cerr << "calculating passing particle potential" << std::endl;
        orbit.setPassing(1, 1, 1); // set passing for Ti/Te = 0.2 TODO take command line
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



