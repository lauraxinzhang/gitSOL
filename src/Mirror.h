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
#include "Struct.h"
#include "Constants.h"

class Mirror
{
    public:

        /**
         * \brief A constructor for uniform magnetic field
         *
         */
        Mirror(double Ti, double Te, double Buniform, double R, double mult, double L, int nx);

        ~Mirror();

        void setGridShift();

        /**
         * \brief Initializes the electric potential in the grid
         * \note  Initializes to 0 potential everywhere if no argument is given
         * 
         */
        void setPotential();

        /**
         * \brief    Initializes potential grid given midplane potential
         * \details  Assume a cos shape for convenience: 
         *             phi (x) = phiMid * cos(x * pi / 2 * xlim)
         *           
         */
        void setPotential(double Rratio, double mult);

        /**
         * \brief Find ambipolar potential at midplane given mirror ratio
         * \param Rratio Mirror ratio at midplane
         */
        double findPhiMid(double Rratio);

        /**
         * \brief Get local potential for arbitrary posion vector pos
         *
         */
        double getPhi(const Vector& pos);

        /**
         * \brief    Calculate electric field everywhere w/ linear finite difference
         * \details  E[i] = - (phi[i+1] - phi[i]) / dx 
         * \note     Potential grid needs to be already setup, if not, calls
         *           setPotential() to setup 0 potential everywhere.
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
         * \brief Gets mod(B) at given location
         */
        double getModB(const Vector& pos);

        /**
         * \brief  Returns grid spacing dx
         */
        double gridSize();

        /**
         * \brief Check whether pos is beyond the computation box.
         */
        bool isLimiter(const Vector& pos);
        
        /**
         * \brief  Calculates which grid the particle is in and increment the 
         *         corresponding value in density_
         */
        void addToBin(Vector& pos);

        /** 
         * \brief Print Vector data along the x axis
         * \param func  A member function with void return type that returns funx(Vector) as a double
         * \param out   Output stream to direct the printed data
         * \note  options include: Bx, phi, modB, density
         */
        void printData(std::string& option, std::ostream &os);



    private:
        Mirror(); // disable default constructor

        // physical parameters first
        double Ti_;  
        double Te_;
        double Buniform_;
        double R_; 
        double L_; 

        // simulation parameters
        double xlim_; // size of computation box (-xlim_, xlim_)
        double ylim_; // (-ylim_, ylim)
        double zlim_; // (-zlim_, zlim)

        int nx_; // size of grid in x direction
        VecDoub * xGrid_;
        VecDoub * xShift_;

        VecDoub * potential_; // potential is only a function of x
        VecDoub * EField_;

        VecDoub * density_; // defined on xShift_

};


#endif // MIRROR_H_INCLUDED
