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
#include "SctmMath.h"
using namespace SctmUtils;

Normalization::Normalization(double temperature)
{
	this->temperature = temperature;
	initFactors();
}

Normalization::Normalization()
{
	this->temperature = SctmPhys::T0;
	initFactors();
}

Normalization::~Normalization(void)
{
}

void Normalization::initFactors()
{
	//this is not used.
	/*
	this->lengthFactor = SctmMath::sqrt(EPSILON * BOLTZMAN * this->temperature / CHARGE / CHARGE / INTRINSIC_CONC_SI);
	this->potentialFactor = BOLTZMAN * this->temperature / CHARGE;
	this->elecFieldFactor = potentialFactor / lengthFactor;
	this->concFactor = INTRINSIC_CONC_SI;
	*/

	const double eps0 = SctmPhys::eps / ( 1 / SctmPhys::cm_in_m); // in [F/cm], important
	const double k0 = SctmPhys::k0;
	const double q = SctmPhys::q;
	const double ni = SctmPhys::ni;

	this->lengthFactor = SctmMath::sqrt(eps0 * k0 * this->temperature / q / q / ni);
	this->potentialFactor = k0 * this->temperature / q;
	this->elecFieldFactor = this->potentialFactor / this->lengthFactor;
	this->concFactor = ni;
	this->diffusionFactor = 1; // in [cm^2/s] the diffusion coefficient is normalized using D0 = 1 cm^2/s
	this->mobilityFactor = this->diffusionFactor / this->potentialFactor;
	this->timeFactor = this->lengthFactor * this->lengthFactor / this->diffusionFactor;
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
