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
#include <map>

typedef std::map<int, double> VertexMapDouble; // <vertID, value>, map of vertex physical value

class FDDomain;
class FDVertex;

namespace SctmUtils
{
	class SctmData;
}

class OneDimSubsSolver
{
	friend class SctmUtils::SctmData;
public:
	enum DopType
	{
		NType, ///< N-type of substrate
		PType, ///< P-type of substrate
	};
	OneDimSubsSolver(FDDomain *_domain);
	void SolveSurfacePot();
	void ReturnResult(VertexMapDouble &_fermiAboveMap, VertexMapDouble &_channelPotMap);
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

	void initializeSolver(); 
	void calcFuncAndItsDeriv(double surfpot);
	double calcFlatbandVoltage(FDVertex *channelVert);
	double calcGateCapacitance(FDVertex *channelVert);
	double solve_NewtonMethod();

	void setFermiAboveCB(FDVertex *channelVert);
	void setChannelPotential(FDVertex *channelVert);

	VertexMapDouble fermiAboveMap;
	VertexMapDouble channelPotMap;
};

#endif