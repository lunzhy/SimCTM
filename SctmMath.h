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
#include <limits>


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
	inline double min_val(double val1, double val2);
	double Bernoulli_Potential(double potVal);
	double SolveQuadEquation(double coeffA, double coeffB, double coeffC, int rootChoice = 0);

	// =====================================================================================
	// IMPLEMENTATIONS
	// =====================================================================================
	inline double abs(double val) {return val > 0 ? val : -val;}
	inline double square(double val) {return val * val;}
	inline double sqrt(double val)  { return std::sqrt(val); }
	inline double ln(double val) { return std::log(val); }
	inline double exp(double val) { return std::exp(val); }
	inline double exp10(double val) { return std::pow(10, val); }
	inline double pow(double val, double arg) { return std::pow(val, arg); }
	inline double asinh(double val) { return std::asinh(val); }
	inline double min_val(double val1, double val2) { return val1 < val2 ? val1 : val2; }
	inline double log10(double val) { return std::log10(val); }
	inline double floor(double val) { return std::floor(val); }
	inline bool logically_equal(double a, double b, double error_factor = 1.0)
	{
		return a == b ||
			abs(a - b) < std::abs(min_val(a, b)) * std::numeric_limits<double>::epsilon() *
			error_factor;
	}

	class VectorValue
	{
	public:
		VectorValue();
		VectorValue(double _vx, double _vy);
		VectorValue &operator=(const VectorValue &_dv);
		friend std::ostream &operator<<(std::ostream &os, const VectorValue &_dv)
		{
			os << "(" << _dv.X() << "," << _dv.Y() << ")";
			return os;
		}
		void Normalize();
		double X() const;
		double Y() const;
		double NormX() const;
		double NormY() const;
		double DirectionValid() const;
	protected:
		double vX;
		double vY;
	};
}

#endif