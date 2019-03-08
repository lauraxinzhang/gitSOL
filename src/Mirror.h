/**
 * \file    Mirror.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    Feb 2019
 *
 * \brief   Declares the Mirror class, a straight mirror field background for iteratively
 *          solving the electrostatic potential.
 * 
 * \note    Uses the same Particle and Vector classes to perform basic operations
 */

#ifndef MIRROR_H_INCLUDED
#define MIRROR_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <stdexcept>

// Programmer class includes
#include "Vector.h"
#include "Particle.h"

// Physical contants. Don't change unless you find the value to be wrong.
// TODO: if both Orbit and Mirror are compiled, are these repeated defs going to give problems?

#define PI 3.14159265358979323846		 // pi, no explanation needed...
#define CC 299792458.0                   // speed of light
// #define QE 4.8E-10                    // electron charge, in Fr

#define ME 9.1E-31                       // electron mass, in kilogram
#define MI 1.6726219E-27                 // ion mass, in kilogram
#define QE 1.60217662E-19                // elementary charge, in coulomb
#define EVTOJOULE 1.60217662E-19         // energy conversion factor

class Mirror
{
	public:

		/**
		 * \brief A constructor for uniform magnetic field
		 *
		 */
		Mirror(double xlim, double ylim, double zlim, double dx, double Buniform);

		~Mirror();

		/**
		 * \brief Initializes the electric potential in the grid
		 * \note  Initializes to 0 potential everywhere if no argument is given
		 * 
		 */
		void setPotential();

		/**
		 * \brief Initializes potential grid given midplane potential
		 */
		void setPotential(double Rratio);

		/**
		 * \brief Find ambipolar potential at midplane given mirror ratio
		 * \param Rratio Mirror ratio at midplane
		 */
		double findPhiMid(double Rratio);

		/**
		 * \brief Calculate electric field everywhere
		 * \note  Potential grid needs to be already setup, if not, calls
		 *        setPotential() to setup 0 potential everywhere.
		 */
		void setEField();

		/** 
		 * \brief Calculate electric field given current position
		 * \param pos Current position of particle
		 *
		 * \note Electric field assumes linearly interpolated potential.
		 */
		Vector getE(const Vector& pos);

		/**
		 * \brief Calculate magnetic field given current position
		 */
		Vector getB(const Vector& pos);

		/**
		 * \brief Return uniform magnetic field if no argument is given
		 */
		Vector getB();

		/**
		 * \brief A dummy function to be called by main.cc to run the simulation.
		 */
		void run();

	private:
		Mirror(); // disable default constructor

		double Ti_;
		double Te_;

		double xlim_; // size of computation box
		double ylim_;
		double zlim_;

		int nx_; // size of grid in x direction
		VecDoub * xGrid_;

		VecDoub * potential_; // potential is only a function of x

		double Buniform_;
};


#endif // MIRROR_H_INCLUDED