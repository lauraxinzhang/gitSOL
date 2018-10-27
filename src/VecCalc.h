/**
 * \file    VecCalc.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2018
 *
 * \brief   Declares the VecCalc class, uitility class for vector calculus
 * 
 */

#ifndef VECCALC_H_INCLUDED
#define VECCALC_H_INCLUDED 1

/**
  *
  *
  */
template < class MAT, class VEC>
class VecCalc
{
	public:
		/**
		 * 
		 */
		VecCalc(); // default constructor

		/**
    	 * \brief Constructor of VecCalc class
    	 * \param lr  length of row
		 * \paran lc  length of column
		 * \details Transfers input to private data members.
    	 */
		VecCalc( const MAT &data , int lr, int lc);

		/**
		 * \brief Calculate the 2-D gradient of data_ (assuming Cartesian), 
		 *        using finite difference in each direction
		 * \details axes corresponds to data[x1][x2].
		 */
		gradFinDiff( MAT &out, const VEC &x1, const VEC &x2 );

		friend std::ostream& operator<<(std::ostream &os, const Vector& v);

	private:
		MAT data_;

};




#endif // VECTOR_H_INCLUDED