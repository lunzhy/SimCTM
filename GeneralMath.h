/**
* @file GeneralMath.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-1   17:13
* @note
* @todo
*/
#ifndef _GENERALMATH_H_
#define _GENERALMATH_H_

#include <cmath>

namespace GeneralMath
{
	const double nm_in_cm = 1e-7;
	// =====================================================================================
	// DECLARATIONS
	// =====================================================================================
	inline double abs(double val);
	inline double square(double val);
	inline double sqrt(double val);

	// =====================================================================================
	// IMPLEMENTATIONS
	// =====================================================================================
	inline double abs(double val) {return val > 0 ? val : -val;}
	inline double square(double val) {return val * val;}
	inline double sqrt(double val)  { return std::sqrt(val); }
}

#endif // !_GENERALMATH_H_