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

class OneDimSubsSolver
{
public:
	enum DopType
	{
		NType, ///< N-type of substrate
		PType, ///< P-type of substrate
	};
	OneDimSubsSolver(FDDomain *_domain);
	void SolveSurfacePot();
	void ReturnResult(double &fermiAbove, double &channelPot);
protected:
	FDDomain *domain;
	double temperature;

	double subsDopConc; ///< substrate doping concentration, in [ni]
	double hDensEqui; ///< p0, equilibrium hole density, in [ni]
	double eDensEqui; ///< n0, equilibrium electron density, in [ni]

	DopType subsType;
	double gateVoltage;
	double flatbandVoltage;
	double gateCapacitance;

	double func_SurfPot;
	double funcDeriv_SurfPot;
	double surfacePotBend;
	double fermiAbove;
	double channelPot;

	void initializeSolver(); 
	void calcFuncAndItsDeriv(double surfpot);
	void calcFlatbandVoltage();
	double solve_NewtonMethod();

	void calcFermiAboveCB();
	void calcChannelPotential();
};

#endif