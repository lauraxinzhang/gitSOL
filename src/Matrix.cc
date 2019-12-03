/**
 * \file    Matrix.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2019
 *
 * \brief   Implements the mathematical Matrix class. Provides vector tensor
 *          products and other arithmetic operations.
 * 
 */

#include "Matrix.h"

Matrix::Matrix()
		: rowOne_(), rowTwo_(), rowThree_()
{
	//nothing to do here.
}

Matrix::Matrix(Vector rowOne, Vector rowTwo, Vector rowThree)
		: rowOne_(rowOne), rowTwo_(rowTwo), rowThree_(rowThree)
{
	// nothing to do here
}

Matrix::Matrix(Matrix& right)
		: rowOne_(right.getRow(0)), rowTwo_(right.getRow(1)), 
		rowThree(right.getRow(2))
{
	// nothing to do here
}

Vector Matrix::getRow(int index) const
{
	switch(index) {
		case 0: return rowOne_;
		case 1: return rowTwo_;
		case 2: return rowThree_;
	}
}

Vector Matrix::dot(Vector& right)
{
	double x = rowOne_.dot(right);
	double y = rowTwo_.dot(right);
	double z = rowThree_.dot(right);

	Vector result(x, y, z);
	return result;
}


Matrix Matrix::operator+(const Matrix& right)
{
    Vector rowOne = getRow(0) + right.getRow(0);
    Vector rowTwo = getRow(1) + right.getRow(1);
    Vector rowThree = getRow(2) + right.getRow(2);
    Matrix result(rowOne, rowTwo, rowThree)
    return result;
}

Matrix Matrix::operator-(const Matrix& right)
{
    Vector rowOne = getRow(0) - right.getRow(0);
    Vector rowTwo = getRow(1) - right.getRow(1);
    Vector rowThree = getRow(2) - right.getRow(2);
    Matrix result(rowOne, rowTwo, rowThree)
    return result;
}

Matrix Matrix::operator/(const double denom)
{
	Vector rowOne = getRow(0) / denom;
    Vector rowTwo = getRow(1) / denom;
    Vector rowThree = getRow(2) / denom;
    Matrix result(rowOne, rowTwo, rowThree)

	return result;
}

Matrix Matrix::operator*(const double mult)
{
	Vector rowOne = getRow(0) * mult;
    Vector rowTwo = getRow(1) * mult;
    Vector rowThree = getRow(2) * mult;
    Matrix result(rowOne, rowTwo, rowThree)

	return result;
}
