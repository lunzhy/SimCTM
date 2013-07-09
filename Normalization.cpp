/**
* @file Normalization.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-9   9:51
* @note
* @todo
*/

#include "Normalization.h"
#include "GeneralMath.h"
using namespace Utility;

Normalization::Normalization(double temperature)
{
	this->temperature = temperature;
	initFactors();
}

Normalization::Normalization()
{
	this->temperature = ROOM_TEMP;
	initFactors();
}

Normalization::~Normalization(void)
{
}

void Normalization::initFactors()
{
	this->lengthFactor = GeneralMath::sqrt(EPSILON * BOLTZMAN * this->temperature / CHARGE / CHARGE / INTRINSIC_CONC_SI);
	this->potentialFactor = BOLTZMAN * this->temperature / CHARGE;
	this->elecFieldFactor = potentialFactor / lengthFactor;
	this->concFactor = INTRINSIC_CONC_SI;
}

void Normalization::ConverseLengthVector( std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction )
{
	switch (direction)
	{
	case Push:
		for (std::vector<double>::size_type ix = 0; ix != real.size(); ++ix)
		{
			norm[ix] = this->PushLength(real[ix]);
		}
		break;
	case Pull:
		for (std::vector<double>::size_type ix = 0; ix != norm.size(); ++ix)
		{
			real[ix] = this->PullLength(norm[ix]);
		}
		break;
	}
}

void Normalization::ConversePotentialVector( std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction )
{
	switch (direction)
	{
	case Push:
		for (std::vector<double>::size_type ix = 0; ix != real.size(); ++ix)
		{
			norm[ix] = this->PushPotential(real[ix]);
		}
		break;
	case Pull:
		for (std::vector<double>::size_type ix = 0; ix != norm.size(); ++ix)
		{
			real[ix] = this->PullPotential(norm[ix]);
		}
		break;
	}
}

void Normalization::ConverseElecFieldVector( std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction )
{
	switch (direction)
	{
	case Push:
		for (std::vector<double>::size_type ix = 0; ix != real.size(); ++ix)
		{
			norm[ix] = this->PushElecField(real[ix]);
		}
		break;
	case Pull:
		for (std::vector<double>::size_type ix = 0; ix != norm.size(); ++ix)
		{
			real[ix] = this->PullElecField(norm[ix]);
		}
		break;
	}
}

void Normalization::ConveseConcVector( std::vector<double> &real, std::vector<double> &norm, ConverseDirection direction )
{
	switch (direction)
	{
	case Push:
		for (std::vector<double>::size_type ix = 0; ix != real.size(); ++ix)
		{
			norm[ix] = this->PushConcentration(real[ix]);
		}
		break;
	case Pull:
		for (std::vector<double>::size_type ix = 0; ix != norm.size(); ++ix)
		{
			real[ix] = this->PullConcentration(norm[ix]);
		}
		break;
	}
}
