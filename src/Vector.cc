/**
 * \file    Vector.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Implements the mathematical Vector class. Provides dot
 *          and cross products and other arithmetic operations.
 * 
 */

#include "Vector.h"

Vector::Vector()
		: x_(0), y_(0), z_(0)
		{
			//nothing to do here.
		}

Vector::Vector(double x, double y, double z)
		: x_(x), y_(y), z_(z)
		{
			// nothing to do here
		}

double Vector::x() const
{
	return x_;
}

double Vector::y() const
{
	return y_;
}

double Vector::z() const
{
	return z_;
}

void Vector::setX(double x)
{
	x_ = x;
}

void Vector::setY(double y)
{
	y_ = y;
}

void Vector::setZ(double z)
{
	z_ = z;
}

Vector Vector::operator+(const Vector& right)
{
    Vector result;
    result.setX(x() + right.x());
    result.setY(y() + right.y());
    result.setZ(z() + right.z());
    return result;
}

Vector Vector::operator-(const Vector& right)
{
    Vector result;
    result.setX(x() - right.x());
    result.setY(y() - right.y());
    result.setZ(z() - right.z());
    return result;
}

Vector Vector::operator/(const double denom)
{
	Vector result;
	result.setX(x() / denom);
	result.setY(y() / denom);
	result.setZ(z() / denom);

	return result;
}

Vector Vector::operator*(const double mult)
{
	Vector result;
	result.setX(x() * mult);
	result.setY(y() * mult);
	result.setZ(z() * mult);

	return result;
}

bool Vector::operator==(const Vector& right)
{
    return x() == right.x() && y() == right.y() && z() == right.z();
}

double Vector::dot(const Vector& right)
{
	return x() * right.x() + y() * right.y() + z() * right.z();
}

Vector Vector::cross(const Vector& right)
{
	Vector result;
	double u1(x());
	double u2(y());
	double u3(z());
	double v1(right.x());
	double v2(right.y());
	double v3(right.z());

	double xx;
	double yy;
	double zz;

	xx = u2 * v3 - u3 * v2;
	yy = u3 * v1 - u1 * v3;
	zz = u1 * v2 - u2 * v1;
	result.setX(xx);
	result.setY(yy);
	result.setZ(zz);

	return result;
}

double Vector::mod()
{
	return sqrt(pow(x(),2) + pow(y(),2) + pow(z(),2));
}

Vector Vector::normalize()
{
	double norm = mod();
	Vector result(x_/norm, y_/norm, z_/norm);
	return result;

}

Vector Vector::parallel(Vector& right)
{
	Vector result;

	double modRight = right.mod();
	double norm = mod();
	double costheta = dot(right) / (modRight * norm);
	double magnitude = norm * costheta;

	result =  right.normalize() * magnitude;
	return result;
}

Vector Vector::perp(Vector& right)
{
	Vector result;
	Vector para = parallel(right);
	result = (*this) - para;
	return result;
}

void Vector::cyl2Cart(const Vector& pos, Vector& vCart)
{
	double xynorm = sqrt( pow(pos.x(),2) + pow(pos.y(),2) );
	if ( xynorm < 1e-10 ) {
		vCart.setX(0);
		vCart.setY(0);
	} else {
		double newX = pos.x()/xynorm * x() - pos.y()/xynorm * y();
		double newY = pos.y()/xynorm * x() + pos.x()/xynorm * y();
		vCart.setX(newX);
		vCart.setY(newY);
	}
	vCart.setZ(z());

	return;
}

Vector Vector::turn(Vector& axis, bool sign)
{
	// assume axis is an already normalized vector; turning (sign * 90) 
	// degrees around axis
	Vector dotted = axis * dot(axis);
	Vector crossed = cross(axis) * (-1); // k x v = - v x k
	if (sign){
		return crossed + dotted;
	} else {
		return dotted - crossed;
	}
}

std::ostream& operator<<(std::ostream &os, const Vector& v) 
{ 
	os << v.x() << ", " << v.y() << ", " << v.z();
	return os;
}

