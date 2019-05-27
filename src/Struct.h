/**
 * \file    Struct.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    March, 2019
 *
 * \brief   Declares the structs for root finding.
 *
 * \NOTE    See header file for detailed explainations for each 
 *          function \n
 *
 *          Uses Numerical Recipe nr3.h data containers for Vector and Matrix
 *          functionalities; uses Vector.h for 3D vectors with arithmetic operations.
 *      
 *          Uses nr3 routines for interpolations.     
 *      
 *          Uses double instead of float to avoid floating number errors.
 */
#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED 1

#include <cassert>
#include "Constants.h"

struct PassingHelp
{
	Doub Ti_;
	Doub Te_;
	Doub R_;

	VecDoub * xGrid_;
	VecDoub * RGrid_;
	MatDoub * integralTable_;

//	friend class Orbit;

	/**
	 * Constructor of the helper struct for finding passing potential
	 */
	PassingHelp(Doub Ti, Doub Te, Doub R);

	~PassingHelp();

	/**
	 * \brief Reads the table of I(x, R) integrals into member integralTable_.
	 * \param input String to the location of integral table.
	 */
	void readTable(std::string input, int nx, int nR);

	/**
	 * \brief Returns Je - Ji for a given phi.
	 * \note  Setup interpolated I(x, R) function from integralTable_.
	 * \return Current difference for current phi.
	 */
	Doub operator() (Doub phi);

};

struct Psir
{
	INTERP2D function_;
	const Doub   z_;

//	friend class Orbit;

	Psir(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub z);

	Doub operator() (Doub r);
};

struct Psiz
{
	INTERP2D function_;
	const Doub   r_;

//	friend class Orbit;

	Psiz(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, Doub r);

	Doub operator() (Doub z);
};

struct psiLimiter
{
	INTERP2D psifull_;
	INTERP1D zOfR_;
	Doub RHS_; // RHS is user supplied psiZero at point(r, z) of the full grid

//	friend class Orbit;

	psiLimiter(const VecDoub& rG, const VecDoub& zG, const MatDoub& f, \
		const VecDoub& rL, const VecDoub& zL);

	void setRHS(Doub r, Doub z);

	// For a given r on the limiter, return the difference in flux 
	// between limiter and given RHS
	Doub operator() (const Doub& r);

	Doub getZ(Doub& r);
};

struct pastukhovHelp
{
	Doub Ti_;
	Doub Te_;

	Doub LHS_;

//	friend class Orbit;

	pastukhovHelp(Doub Ti, Doub Te, Doub R);

	Doub operator()(const Doub& x);
};

/**
 * \brief constructs a linearly interpolated potential along an given l axis
 * \return phi(l). Callable for use in differentiation.
 */
struct eFieldHelp
{
	INTERP1D helper_; // shouldn't have any problem if initialized under ':'
	// VecDoub lList_, potList_;

//	friend class Orbit;

	eFieldHelp(VecDoub lList, VecDoub potList);

	Doub operator() (Doub l);
};

struct Histogram
{
	Doub min_, max_, gridsize_;
	VecDoub * bins_;
	
	Histogram(Doub min, Doub max, int numBin);

	void addToBin(Doub val);

	friend class Mirror;
};

#endif

