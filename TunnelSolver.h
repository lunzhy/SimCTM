/**
* @file TunnelSolver.h
* @brief This file contains the classes that account for solving tunneling problems
*
*
*
* @author
* @version 
* @date 2013-7-16   10:39
* @note
* @todo
*/
#ifndef _TUNNELSOLVER_H_
#define _TUNNELSOLVER_H_

#include <vector>
class FDVertex;

using std::vector;
class TunnelSolver
{
	//friend class SctmUtils::SctmDebug;
public:
	TunnelSolver();
	virtual void PrepareProblem(FDVertex *startVertex) = 0;
	void SolveTunneling();
	double GetCurrentDensity();

protected:
	virtual double getSupplyFunction(double energy);
	double getTransmissionCoefficient(double energy);
	void initializeSolver();
	void calcDTFNtunneling();
	void calcThermalEmission();

	vector<double> cbegde;
	vector<double> emass;
	vector<double> deltaX;
	double cbedgeTunnelFrom; ///< left electrode conduction band edge
	double cbedgeTunnelTo;
	double fermiEnergyTunnelFrom;
	double fermiEnergyTunnelTo;
	double effTunnelMass; ///< effective mass
	double temperature;
	double areaFactor; ///< the area factor(cross-section) of the tunneling problem

	double currentDensity; ///< the tunneling current density, in [A/cm^2]
};

class SubsToGateEletronTunnel : public TunnelSolver
{
public:
	void PrepareProblem(FDVertex *startVertex);
};

class TunnelTest : public TunnelSolver
{
public:
	TunnelTest();
	void PrepareProblem(FDVertex *startVertex);
	void SolveParamterSet();
	void SolveCalibrate();
private:
	double oxideEmass;
	double siliconBandEdge;
	double oxideThickness;
	double gateVoltage;
	double elecField;
};

#endif