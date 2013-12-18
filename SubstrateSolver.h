/**
* @file SubstrateSolver.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-12-16   21:22
* @note
* @todo
*/

#ifndef _SUBSTRATESOLVER_H_
#define _SUBSTRATESOLVER_H_

class FDDomain;

class SubstrateSolver
{
public:
	enum SubsType
	{

	};
	SubstrateSolver(FDDomain *_domain);
	void SolveSurfacePot();
protected:
	FDDomain *domain;
	double temperature;

	double subsDopConc; ///< substrate doping concentration
	double hDensEqui; ///< p0, equilibrium hole density, in [ni]
	double eDensEqui; ///< n0, equilibrium electron density, in [ni]

	double gateVoltage;
	double flatbandVoltage;
	double gateCapacitance;

	double func_SurfPot;
	double funcDeriv_SurfPot;

	void initializeSolver(); 
	void calcFuncAndItsDeriv(double surfpot);
};

#endif