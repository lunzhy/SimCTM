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
#pragma once
#include <cmath>
#include "DomainDetails.h"

namespace GeneralMath
{
	const double nm_in_cm = 1e-7;
	// =====================================================================================
	// DECLARATIONS
	// =====================================================================================
	inline double abs(double val);
	inline double square(double val);
	inline double sqrt(double val);
	inline double distance(FDVertex *vertex1, FDVertex *vertex2);

	// =====================================================================================
	// IMPLEMENTATIONS
	// =====================================================================================
	inline double abs(double val) {return val > 0 ? val : -val;}
	inline double square(double val) {return val * val;}
	inline double sqrt(double val)  { return std::sqrt(val); }
	inline double distance(FDVertex *vertex1, FDVertex *vertex2)
	{
		return sqrt(square(vertex1->X - vertex2->X) + square(vertex1->Y - vertex2->Y));
	}
}