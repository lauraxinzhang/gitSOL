/**
 * \file    Matrix.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2019
 *
 * \brief   Declares the Matrix class
 * \details A 3D matrix is stored as 3 vectors
 */

#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <list>
#include <sstream>

class Matrix
{
	public:
		Matrix();

		Matrix(Vector rowOne, Vector rowTwo, Vector rowThree);

		Matrix(Matrix& right);

		Vector getRow(int index) const;

		Vector dot(Vector& right);

		/**
		 * \brief Vector - scalar arithmetics
		 */
		Matrix operator+(const Matrix& right);
		Matrix operator-(const Matrix& right);
		Matrix operator/(const double denom);
		Matrix operator*(const double mult);

	private:
		Vector rowOne_;
		Vector rowTwo_;
		Vector rowThree_;


}