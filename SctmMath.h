/**
* @file SctmMath.h 
* @brief This file contains the general math problems used in the simulation.
*
*
*
* @author
* @version 
* @date 2013-7-1   17:13
* @note
* @todo
*/
#ifndef _SCTMMATH_H_
#define _SCTMMATH_H_

#define _USE_MATH_DEFINES
//#ifndef _GENERALMATH_H_
//#define _GENERALMATH_H_

#include <math.h>
#include <cmath>
#include <iostream>


/// @brief SctmMath is a namespace in which the math methods are defined for this simulation.
///
/// 
/// This class is designed to provide probable optimization of the math function used in the simulation. 
/// Currently, some of the implementations are simply calling the method from cmath library.
namespace SctmMath
{
	const double PI = M_PI;
	// =====================================================================================
	// DECLARATIONS
	// =====================================================================================
	inline double abs(double val);
	inline double square(double val);
	inline double sqrt(double val);
	inline double ln(double val);
	inline double exp(double val);

	// =====================================================================================
	// IMPLEMENTATIONS
	// =====================================================================================
	inline double abs(double val) {return val > 0 ? val : -val;}
	inline double square(double val) {return val * val;}
	inline double sqrt(double val)  { return std::sqrt(val); }
	inline double ln(double val) { return std::log(val); }
	inline double exp(double val) { return std::exp(val); }

	class DirectionVector
	{
	public:
		DirectionVector(){ this->vX = 0; this->vY = 0; }
		DirectionVector(double _vx, double _vy)
		{
			if ( (_vx == 0) && (_vy == 0) )
			{
				this->vX = 0;
				this->vY = 0;
			}
			else
			{
				this->vX = _vx / sqrt( square(_vx) + square(_vy) );
				this->vY = _vy / sqrt( square(_vx) + square(_vy) );
			}
		}
		DirectionVector &operator=(const DirectionVector &_dv)
		{
			this->vX = _dv.X();
			this->vY = _dv.Y();
			return *this;
		}
		friend std::ostream &operator<<(std::ostream &os, const DirectionVector &_dv)
		{
			os << "(" << _dv.X() << "," << _dv.Y() << ")";
			return os;
		}
		double X() const { return this->vX; }
		double Y() const { return this->vY; }
		double Valid() const { return ( (this->vX != 0)||(this->vY != 0) ); }
	protected:
		double vX;
		double vY;
	};
}

#endif