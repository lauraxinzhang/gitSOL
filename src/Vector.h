/**
 * \file    Vector.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    August, 2017
 *
 * \brief   Declares the Vector class
 * 
 */

#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <list>
#include <sstream>

/**
  *
  *
  */

class Vector
{
	public:
		/**
		 * \brief Overloading default constructor. 
		 * \details set all components to 0
		 */
		Vector(); // default constructor

		/**
    	 * \brief Constructor of Vector class
    	 * \param x  x component of the vector
    	 * \param y  y component of the vector
    	 * \param z  z component of the vector
		 *
		 * \details Transfers input to private data members.
    	 */
		Vector(double x, double y, double z);

		double x() const;
		double y() const;
		double z() const;

		void setX(double x);
		void setY(double y);
		void setZ(double z);

		Vector operator+(const Vector& right);
		Vector operator-(const Vector& right);
		Vector operator/(const double denom);
		Vector operator*(const double mult);
		bool operator==(const Vector& right);

		double dot(const Vector& right);
		Vector cross(const Vector& right);
		double mod();
		Vector normalize();

		Vector parallel(Vector& right);
		Vector perp(Vector& right);

		void cyl2Cart(const Vector&pos, Vector& vCart);
		
		/**
		 * \brief Rotate the pitch angle of self with respect to axis by 90 degrees
		 * \param axis The axis of projection; This is usually magnetic field B for pitch 
		 *             angle scattering of particle velocities
		 * \param sign The (+-) sign of the 90 degree angle.
		 */
		Vector turn(const Vector& axis, bool sign)

		friend std::ostream& operator<<(std::ostream &os, const Vector& v);

	private:

		double x_;
		double y_;
		double z_;


};




#endif // VECTOR_H_INCLUDED
