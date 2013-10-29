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
			vX = vX / sqrt( square(vX) + square(vY) );
			vY = vY / sqrt( square(vX) + square(vY) );
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

}
