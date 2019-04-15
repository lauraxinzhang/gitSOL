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
#include <cassert>
#include <fstream>
#include <string>

// Programmer class includes
#include "Vector.h"
#include "Particle.h"
#include "Constants.h"

class Mirror
{
	public:

		/**
		 * \brief A constructor for uniform magnetic field
		 *
		 */
		Mirror(double Ti, double Te, double xlim, double ylim, double zlim, int nx, double Buniform);

		~Mirror();

		void setGridShift();

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
		 * \brief Check whether pos is beyond the computation box.
		 */
		bool isLimiter(const Vector& pos);
		

	private:
		Mirror(); // disable default constructor

		// physical parameters first
		double Ti_;  
		double Te_;
		double Buniform_;
		double R_; 
		double L_; 

		// simulation parameters
		double xlim_; // size of computation box (0, xlim_)
		double ylim_; // (-ylim_, ylim)
		double zlim_; // (-zlim_, zlim)

		int nx_; // size of grid in x direction
		VecDoub * xGrid_;
		VecDoub * xShift_;

		VecDoub * potential_; // potential is only a function of x




};


#endif // MIRROR_H_INCLUDED
