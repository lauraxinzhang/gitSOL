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

class Matrix;
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

		/**
		 * \brief Copy constructor
		 */
//		Vector(Vector& right);
		/**
		 * \brief Data getter for x component
		 */
		double x() const;

		/**
		 * \brief Data getter for y component
		 */
		double y() const;

		/**
		 * \brief Data getter for z component
		 */
		double z() const;

		/**
		 * \brief Data setter for x component
		 */
		void setX(double x);

		/**
		 * \brief Data setter for y component
		 */
		void setY(double y);

		/**
		 * \brief Data setter for z component
		 */
		void setZ(double z);


		/**
		 * \brief Vector - scalar arithmetics
		 */
		Vector operator+(const Vector& right);
		Vector operator-(const Vector& right);
		Vector operator/(const double denom);
		Vector operator*(const double mult);

		/**
		 * \brief Comparison operator
		 */
		bool operator==(const Vector& right);

		/**
		 * \brief Vector - Vector arithmetics
		 */
		double dot(const Vector& right);
		Vector cross(const Vector& right);

		/**
		 * \brief Vector tensor product
		 */
		Matrix tensor(Vector& right);

		double mod();
		Vector normalize();

		Vector parallel(Vector& right);
		Vector perp(Vector& right);

		/**
		 * \brief Converts a (R, phi, Z) vector to (x, y, z)
		 * \note  Used to convert field vectors defined on (R, Z) to cartesian 
		 *        coordinated used by particle push.
		 */
		void cyl2Cart(const Vector&pos, Vector& vCart);


		/**
		 * \brief Converts a (x, y, Z) vector to (R, phi, z)
		 * TODO: implement this
		 */
		void cart2Cyl(const Vector&pos, Vector& vCyl);
		
		/**
		 * \brief Rotate the pitch angle of self with respect to axis by 90 degrees
		 * \param axis The axis of projection; This is usually magnetic field B for pitch 
		 *             angle scattering of particle velocities
		 * \param sign The (+-) sign of the 90 degree angle.
		 * \details Rotation is done with Rodrigues formula at 90 degrees.
		 */
		Vector turn(Vector& axis, bool sign);

		/**
		 * \brief Damp (reduce) the speed of the particle according to slowing down
		 *        frequency nu_s
		 * \param nu_s Slowing down frequency calculated from thermal plasma profile
		 * \param dt   Time step of particle push.
		 * \details    v_(t+1) = (1 - dt * nu_s) * v_t
		 */
		Vector damp(double nu_s, double dt);

		friend std::ostream& operator<<(std::ostream &os, const Vector& v);

	private:

		double x_;
		double y_;
		double z_;


};




#endif // VECTOR_H_INCLUDED
