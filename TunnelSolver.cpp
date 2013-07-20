/**
* @file TunnelSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-19   10:31
* @note
* @todo
*/

#include "TunnelSolver.h"
#include "SctmMath.h"
#include "SctmPhys.h"

double TunnelSolver::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;
	double EfTunnelTo = this->fermiEnergyTunnelTo;

	double kB = SctmPhys::BoltzmanConstant;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 1;

	integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(-(energy - EfTunnelFrom) / kB / T));
	integralTunnelTo = kB * T * SctmMath::ln(1 + SctmMath::exp(-(energy - EfTunnelTo) / kB / T));
	double supply = integralTunnelFrom - integralTunnelTo;

	return supply;
}

double TunnelSolver::getTransmissionCoefficient(double energy)
{
	double tunnelFactor = 1;
	// The transmission coefficient is proportional to the integral used in the calculation.
	// The difference/factor is used to describe the reflection of the tunneling wave.
	// Here we assume that the reflection is tiny in this problem.

	double integral = 0;
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		integral += SctmMath::sqrt(2 * emass.at(ix) * (cbegde.at(ix) - energy)) * deltaX.at(ix);
	}
	double tc = tunnelFactor * SctmMath::exp(-2 / SctmPhys::hbar * integral);

	return tc;
}

void TunnelSolver::calc_DT_FN_Tunneling()
{
	double currentDensity = 0; // in A/m^3
	double Emin = cbedgeTunnelFrom; // in eV
	double Emax = Emin + 10;

	double kB = SctmPhys::BoltzmanConstant;
	double h = SctmPhys::h;
	double T = this->temperature;
	double q = SctmPhys::ElementaryCharge;
	double kT = kB * T; // dE = kT/q
	double dE = kT; 
	double pi = SctmMath::PI;

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		currentDensity += 4 * pi * this->effTunnelMass * q / h / h / h
			* getTransmissionCoefficient(currEnergy)
			* getSupplyFunction(currEnergy)
			* dE;
		currEnergy += dE;
	}
}

void SubsToGateEletronTunnel::PrepareProblem(FDVertex *startVertex)
{
	// the tunneling direction is north
	this->areaFactor = startVertex->EastLength / 2 + startVertex->WestLength / 2;
	FDVertex *currentVertex = startVertex;
	while (currentVertex->NorthVertex != NULL)
	{
		this->deltaX.push_back((currentVertex->NorthLength + currentVertex->SouthLength) / 2);
		this->cbegde.push_back(currentVertex->Phys.ConductionBandEnergy);
		this->emass.push_back(currentVertex->Phys.ElectronMass);
		currentVertex = currentVertex->NorthVertex; // move the vertex north
	}
}
