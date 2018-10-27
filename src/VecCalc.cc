/**
 * \file    VecCalc.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2018
 *
 * \brief   Implements the VecCalc class, uitility class for vector calculus
 * 
 */

#include "VacCalc.h"

template < class MAT, class VEC >
VecCalc::VecCalc( const MAT &data )
	: data_( data )
{
	//type MAT is assumed to have a valid copy constuctor.
	// nothing else to do here.
}

template < class MAT, class VEC >
VecCalc::gradFinDiff( MAT &out, const VEC &x1, const VEC &x2 )
{
	VEC x1Shift( x1.begin(), --( x2.end()) );
	VEC x1Shift( x1.begin(), --( x2.end()) );
}