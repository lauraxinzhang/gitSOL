/**
 * \file    Orbit.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Declares the Orbit class, toy model for particle orbits in LTX
 *          SOL w/ Pastukhov potential solved analytically.
 * 
 */

#ifndef ORBIT_H_INCLUDED
#define ORBIT_H_INCLUDED 1

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath> 
#include <random>
#include <stdexcept>
#include <list>
#include <sstream>
#include <string>
#include <stdio.h>
#include <cctype>
#include <cassert>
#include <math.h>
#include <omp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

// Programmer class includes
#include "Vector.h"
#include "Particle.h"

// Numerical Recipe routines includes
#include "nr3.h"
#include "interp_1d.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "dfridr.h"
#include "roots.h"


// Physical contants. Don't change unless you find the value to be wrong.
#define PI 3.14159265358979323846		 // pi, no explanation needed...
#define CC 299792458.0                   // speed of light
// #define QE 4.8E-10                    // electron charge, in Fr

#define ME 9.1E-31                       // electron mass, in kilogram
#define MI 1.6726219E-27                 // ion mass, in kilogram
#define QE 1.60217662E-19                // elementary charge, in coulomb
#define EVTOJOULE 1.60217662E-19         // energy conversion factor

// Numbers specific to LTX geometry
#define BMAGAXIS  2000                   // mod(B) = 0.2 T = 2000 Gauss at magnetic axis, characteristic field strength
#define NPERORBIT 20                     // steps per Lamor orbit, from Boris convergence test.


/** Choose which interpolater to use throughout the class */
// Use these two lines to compile with linear interpolations

typedef Linear_interp INTERP1D;
typedef Bilin_interp  INTERP2D;


// Use these two lines to compile with Spline interpolations 
// (very slow. use only in production mode when very pretty pictures are needed)
// typedef Spline_interp    INTERP1D;
// typedef Spline2D_interp  INTERP2D;

/**
 * \brief  
 *
 * \NOTE   
 */
class Orbit
{
	public:
		/**
		 * \brief A default constructor to produce a simple field for testing
		 */
		Orbit();

		/**
		 * \brief Another constructor to read uniform grid magnetic field
		 * \param field Path of input field profile.
		 * \param limiter Path of limiter coordinates.
		 */
		Orbit(const std::string& field, const std::string& limiter);

		/**
		 * \brief Read limiter file and transfer to data member.
		 * \note  Called by the constructor.
		 */
		void readLimiter(std::ifstream &input);

		/**
		 * \brief Destructor. Manually delete dynamically allocated memory.
		 */
		~Orbit();

		/**
		 * \brief Configure the mirror ratio matrix. Fill in data to RRatio_ private member
		 * \note  Need to be called when mirror ratio is needed.
		 */
		void configMirror();

		/**
		 * \brief Calculate the mirror ratio on grid points.
		 * \param ratio An empty MatDoub reference to write results in.
		 *
		 * \details For each grid point, use root finding to find the corresponding (r, z)
		 *          on the limiter with the same flux value. Then calculate mirror ratio with
		 *          respect to that point.
		 *
		 * \note    TODO: deal with the up-down asymmetry. Define mirror ratio respect to the 
		 *          lower field on the limiter.
		 */
		void mirrorRatio(MatDoub& ratio);

		/**
		 * \brief Get the magnetic field at a position
		 * \param pos Current position as a Vector class member in cartesian coordinate.
		 * \return Magnetic field as a 3D Vector in cylindrical coordinate.
		 */
		Vector getB(const Vector& pos);

		/**
		 * \brief Get mod(B) for a position given in the RZ plane
		 * \param rr R position.
		 * \param zz Z position.
		 */
		Doub getModB(Doub rr, Doub zz);

		Doub getMirrorRatio(Doub rr, Doub zz);

		/**
		 * \brief Check whether input rr, zz is on or beyond the limiter
		 * \note  TODO maybe should only retun *beyond* to allow for the points exactly on limiter?
		 */
		bool isLimiter(Doub rr, Doub zz);

		/**
		 * \brief Returns the vertical coordinate of limiter surface for any r
		 */
		Doub getzLimiter(Doub rr);


		//--------------------------------------------------------------------------------------------
		//--------------------------- End of Basic Class Implementations -----------------------------
		//------------------------------ Start potential calculations --------------------------------
		//--------------------------------------------------------------------------------------------

		/**
		 * \brief Calculate the pastukhov potential for given input parameters
		 * \details Uses nr root finding to solve for normalized potential. RHS of root finding
		 *          is implemented in Struct pastukhovHelp.
		 */
		Doub pastukhov(Doub Ti, Doub Te, Doub R);

		/**
		 * \brief Calculates pastukhov and writes to private member Phi_
		 * \note configMirror() is prerequired. TODO check Rmirror_ == nullptr.
		 */
		void setPastukhov(Doub Ti, Doub Te, Doub multiplier = 1);


		//--------------------------------------------------------------------------------------------
		//------------------------------ End of potential calculations -------------------------------
		//---------------------------------- Start data outputting -----------------------------------
		//--------------------------------------------------------------------------------------------
		
		/** 
		 * \brief Starting from r0, phi0, z0, trace out a magnetic field line
		 *        Advance position in the direction of B, length dl at a time.
		 */
		void fieldLines(std::ofstream &rBlist, std::ofstream &zBlist, std::ofstream &phiBlist, const Vector& init, Doub dl, int iter);

		/**
		 * \brief Calculate the Pastukhov potential along a single field line, 
		 *        outputs arc length along field line l, ephi/Te, and electric field
		 */
		void eField(Doub Ti, Doub Te, Vector init, Doub dl, int iter, std::ofstream &output);

		/**
		 * \brief Writes {R, Z, |B|} to output
		 */
		void writeModB(std::ofstream &output);


		/**
		 * \brief Writes {R, Z, psi} to output
		 */
		void writeFlux(std::ofstream &output);

		/**
		 * \brief Writes {R, Z, Rmirror} to output
		 */
		void writeMirrorRatio(std::ofstream &output);

		/**
		 * \brief Writes {R, Z, ephi/Te} to output
		 */
		void writePastukhov(Doub Ti, Doub Te, std::ofstream &output);

		/**
		 * \brief Writes potential on the mid plane
		 */
		void midPlanePhi(Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output);

		/**
		 * \brief Outputs potential as a function of ion temperature.
		 */
		void temperature(Doub Ti_start, Doub dT, int iter, Doub R);

		/**
		 * \brief Print field data to file.
		 */
		void printData();

		/**
		 * \brief Pushes a particle in background field
		 * \output Writes coordinates in RZ and XYZ, energy, and mu to file
		 * \details outputs a number every 500 steps
		 * TODO: fix this magic number 500
		 */
		void particlePush(Doub dr, Doub energy, Doub er, Doub ephi, Doub ez);

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
		Doub particleStats(Doub dr, Doub energy, bool spec, int nparts, int maxiter, bool write = false);


		/**
		 * \brief A unit testing suite
		 */
		void test();

		/**
		 * \brief Constructs a test orbit with straight field lines.
		 */
		void emptytest();



	private:
	    
	    // declare member as pointers, objects are created on the heap and deleted
	    // in the destruction operator.
	    MatDoub * Br_;
	    MatDoub * Bz_;
	    MatDoub * Btor_; 
	    MatDoub * Bmod_;
	    MatDoub * Rratio_; // initialized as nullptr until Orbit::configMirror is called.
	    MatDoub * Phi_;  // initialized as nullptr until Orbit::setPastukhov is called.

	    MatDoub * psiRZ_;
	    MatDoub * pRZ_;

	    VecDoub * rGrid_; // r grid points
	    VecDoub * zGrid_; // z grid points

	    VecDoub * rLimit_; // r coordinates of limiter
	    VecDoub * zLimit_; // z coordinates of limiter

		int nw_, nh_;

	    Doub rdim_, zdim_, rleft_, zmid_;
	    Doub rllmtr_, rrlmtr_;
	    Doub zllmtr_, zrlmtr_;


};

// Helper Functors for numerical routines

struct Psir
{
	INTERP2D function_;
	const Doub   z_;

	Psir(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub z)
		:function_(rG, zG, f), z_(z)
	{
		/** interpolation objects need to be initialized under the colon 
		 initializer. Otherwise the data members will be messed up. */
		// function_ = function;
	}

	Doub operator() (Doub r)
	{
		return function_.interp(r, z_);
	}
};

struct Psiz
{
	INTERP2D function_;
	const Doub   r_;

	Psiz(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub r)
		:function_(rG, zG, f), r_(r)
	{
		// function_ = function;
	}

	Doub operator() (Doub z)
	{
		return function_.interp(r_, z);
	}
};

struct psiLimiter
{
	INTERP2D psifull_;
	INTERP1D zOfR_;
	Doub RHS_; // RHS is user supplied psiZero at point(r, z) of the full grid

	psiLimiter(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, \
		const VecDoub& rL, const VecDoub& zL)
		: psifull_(rG, zG, f), zOfR_(rL, zL), RHS_(0)
	{
		// nothing to do here?
	}

	void setRHS(Doub r, Doub z)
	{
		RHS_ = psifull_.interp(r, z);
		return;
	}

	// For a given r on the limiter, return the difference in flux 
	// between limiter and given RHS
	Doub operator() (const Doub& r)
	{
		Doub zNow( zOfR_.interp(r) );

		Doub psiLimiter( psifull_.interp(r, zNow) );

		return psiLimiter - RHS_;
	}

	Doub getZ(Doub& r)
	{
		return zOfR_.interp(r);
	}
};

struct pastukhovHelp
{
	Doub Ti_;
	Doub Te_;

	Doub LHS_;

	pastukhovHelp(Doub Ti, Doub Te, Doub R)
		:Ti_(Ti), Te_(Te)
	{
		Doub constant   = sqrt(2.0/PI)*sqrt(MI/ME) * pow( (Ti/Te), (3.0/2.0) );
		// MI and ME defined in particle.cc
		Doub gOfR       = (2.0*R+1.0)/(2.0*R) * log(4.0*R+2.0);
		Doub LHS        =  constant * log10(R)/gOfR;

		LHS_ = LHS;
	}

	Doub operator()(const Doub& x)
	{
		Doub I = 1 + ( sqrt(PI) * exp(x) * erfc( sqrt(x) ) / (2 * sqrt(x)) );
		Doub fOfX = x * exp(x)/I;

		return fOfX - LHS_;
	}
};

/**
 * \brief constructs a linearly interpolated potential along an given l axis
 * \return phi(l). Callable for use in differentiation.
 */
struct eFieldHelp
{
	INTERP1D helper_; // shouldn't have any problem if initialized under ':'
	// VecDoub lList_, potList_;

	eFieldHelp(VecDoub lList, VecDoub potList)
		:helper_(lList, potList)
	{
		// lList_ = lList;
		// potList_ = potList;

		//nothing else to do
	}

	Doub operator() (Doub l)
	{
		// INTERP1D helper(lList_, potList_);
		return helper_.interp(l);
	}
};

#endif // ORBIT_H_INCLUDED