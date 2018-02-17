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
#include <stdlib.h>
#include <cctype>
#include <cassert>
#include <math.h>

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
#define PI 3.14159265358979323846
#define CC 299792458.0                       // speed of light
// #define QE 4.8E-10                        // electron charge, in Fr

#define ME 9.1E-31                        // electron mass, in kilogram
#define MI 1.6726219E-27                  // ion mass, in kilogram
#define QE 1.60217662E-19                 // elementary charge, in coulomb
#define EVTOJOULE 1.60217662E-19          // conversion factor

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

		// *
		//  * \brief A constructor to read Eqdsk magnetic field
		//  * \note  Uses a dummy int to choose constructor
		 
		// Orbit(const std::string& path, int flag);


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
//---------------------- Start potential calculations and outputting -------------------------
//--------------------------------------------------------------------------------------------

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
		 * \brief Calculate the pastukhov potential for given input parameters
		 * \details Uses nr root finding to solve for normalized potential. RHS of root finding
		 *          is implemented in Struct pastukhovHelp.
		 */
		Doub pastukhov(Doub Ti, Doub Te, Doub R);


		// some outputting helpers

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
		void writePastukhov(Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output);


		/**
		 * \brief Writes potential on the mid plane
		 */
		void midPlanePhi(Doub Ti, Doub Te, MatDoub& pas, std::ofstream &output);


		/**
		 * \brief A unit testing suite
		 */
		void test();

		/**
		 *
		 */
		void emptytest();

		// public data members
		int nw_, nh_;

	    Doub rdim_, zdim_, rleft_, zmid_;
	    Doub rllmtr_, rrlmtr_;
	    Doub zllmtr_, zrlmtr_;

	private:
	    
	    // declare member as pointers, objects are created on the heap and deleted
	    // in the destruction operator.
	    MatDoub * Br_;
	    MatDoub * Bz_;
	    MatDoub * Btor_; 
	    MatDoub * Bmod_;

	    MatDoub * psiRZ_;
	    MatDoub * pRZ_;

	    VecDoub * rGrid_; // r grid points
	    VecDoub * zGrid_; // z grid points

	    VecDoub * rLimit_; // r coordinates of limiter
	    VecDoub * zLimit_; // z coordinates of limiter
};

// Helper Functors for numerical routines

struct Psir
{
	Bilin_interp function_;
	const Doub   z_;

	Psir(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub z)
		:function_(rG, zG, f), z_(z)
	{
		// function_ = function;
	}

	Doub operator() (Doub r)
	{
		return function_.interp(r, z_);
	}
};

struct Psiz
{
	Bilin_interp function_;
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
	Bilin_interp psifull_;
	Linear_interp zOfR_;
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

struct eFieldHelp
{
	// Linear_interp helper_;
	VecDoub lList_, potList_;

	eFieldHelp(VecDoub lList, VecDoub potList)
	{
		lList_ = lList;
		potList_ = potList;
		//nothing else to do
	}

	Doub operator() (Doub l)
	{
		Linear_interp helper(lList_, potList_);
		return helper.interp(l);
	}
};

#endif // ORBIT_H_INCLUDED