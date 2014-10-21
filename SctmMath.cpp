/**
* @file SctmMath.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-10-29   16:05
* @note
* @todo
*/

#include "SctmMath.h"
#include "SctmUtils.h"
#include "SctmPhys.h"
#include "Normalization.h"

namespace SctmMath
{
	VectorValue::VectorValue()
	{
		this->vX = 0; this->vY = 0;
	}

	VectorValue::VectorValue(double _vx, double _vy)
	{
		this->vX = _vx;
		this->vY = _vy;
	}

	double VectorValue::DirectionValid() const
	{
		return ( (this->vX != 0)||(this->vY != 0) );
	}

	double VectorValue::NormY() const
	{
		if ( (vX == 0) && (vY == 0) )
		{
			//keep (0, 0) vector for valid inspection
			SCTM_ASSERT(SCTM_ERROR, 10018);
			return 0;
		}
		else
			return vY / sqrt( square(vX) + square(vY) );
	}

	double VectorValue::NormX() const
	{
		if ( (vX == 0) && (vY == 0) )
		{
			//keep (0, 0) vector for valid inspection
			SCTM_ASSERT(SCTM_ERROR, 10018);
			return 0;
		}
		else
			return vX / sqrt( square(vX) + square(vY) );
	}

	void VectorValue::Normalize()
	{
		if ( (vX == 0) && (vY == 0) )
		{
			//keep (0, 0) vector for valid inspection
			SCTM_ASSERT(SCTM_ERROR, 10018);
		}
		else
		{
			double tempX = vX;
			double tempY = vY;
			vX = tempX / sqrt( square(tempX) + square(tempY) );
			vY = tempY / sqrt( square(tempX) + square(tempY) );
		}
	}

	VectorValue & VectorValue::operator=(const VectorValue &_dv)
	{
		this->vX = _dv.X();
		this->vY = _dv.Y();
		return *this;
	}

	double VectorValue::X() const
	{
		return this->vX;
	}

	double VectorValue::Y() const
	{
		return this->vY;
	}

	double Bernoulli_Potential(double potVal)
	{
		using SctmUtils::SctmGlobalControl;
		double ret = 0;

		double temperature = SctmGlobalControl::Get().Temperature;
		double kT_div_q = SctmPhys::k0 * temperature / SctmPhys::q;
		
		using SctmUtils::Normalization;
		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		
		double val = norm.PullPotential(potVal) / kT_div_q;
		double deno = exp(val) - 1; // denominator

		if (deno == 0)
		{
			ret = 1;
		}
		else
		{
			ret = val / deno;
		}
		return ret;
	}

	double SolveQuadEquation(double coeffA, double coeffB, double coeffC, int rootChoice /*= 0*/)
	{
		double sqrtvalue = 0;
		double root = 0;

		sqrtvalue = SctmMath::sqrt(coeffB * coeffB - 4 * coeffA * coeffC);
		if (sqrtvalue < 0)
			SCTM_ASSERT(SCTM_ERROR, 10057);

		if (rootChoice == 0)
		{
			root = (-coeffB - sqrtvalue) / 2 / coeffA;
		}
		else
		{
			root = (-coeffB + sqrtvalue) / 2 / coeffA;
		}

		return root;
	}



}
