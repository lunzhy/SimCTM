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
#pragma once
#include <vector>
#include "DomainDetails.h"
using namespace std;
class TunnelSolver
{
public:
	virtual void PrepareProblem(FDVertex *startVertex) = 0;

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
};

class SubsToGateEletronTunnel : TunnelSolver
{
public:
	void PrepareProblem(FDVertex *startVertex);
};