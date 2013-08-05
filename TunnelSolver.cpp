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

double TunnelSolver::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;
	double EfTunnelTo = this->fermiEnergyTunnelTo;

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 1;

	integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelFrom) / kB / T));
	integralTunnelTo = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelTo) / kB / T));
	double supply = integralTunnelFrom - integralTunnelTo;

	return supply;
}

double TunnelSolver::getTransmissionCoefficient(double energy)
{
	double tunnelFactor = 1;
	// The transmission coefficient is proportional to the integral used in the calculation.
	// The difference/factor is used to describe the reflection of the tunneling wave.
	// Here we assume that the reflection is tiny in this problem.

	double q = SctmPhys::ElementaryCharge;
	double hbar = SctmPhys::hbar;
	double integral = 0;
	//direct tunneling and Fowler-Nordheim tunneling are included.
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		if (cbegde.at(ix) >= energy)
		{
			integral += SctmMath::sqrt(2 * emass.at(ix) * q * (cbegde.at(ix) - energy)) * deltaX.at(ix);
		}
	}
	double tc = tunnelFactor * SctmMath::exp(-2 / hbar * integral);

	return tc;
}

void TunnelSolver::calcDTFNtunneling()
{
	// the unit of calculated current density is [A/m^2]
	double currentDensity = 0; // in [A/m^2]
	double Emin = cbedgeTunnelFrom; // in [eV]
	double Emax = cbegde.front(); // in [eV]

	double T = this->temperature;
	double pi = SctmMath::PI;
	double h = SctmPhys::h;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T; // dE = kT

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

void TunnelSolver::calcThermalEmission()
{
	double currentDensity = 0; // in [A/m^2]
	double Emin = cbegde.front(); // in [eV]
	double Emax = Emin + 5; // only calculate the energy within 5eV large than the smallest energy to surpass the barrier

	double T = this->temperature;
	double pi = SctmMath::PI;
	double h = SctmPhys::h;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T; // dE = kT

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		currentDensity += 4 * pi * this->effTunnelMass * q / h / h / h
			* 1  //for transmission coefficient
			* getSupplyFunction(currEnergy)
			* dE;
		currEnergy += dE;
	}
}

void SubsToGateEletronTunnel::PrepareProblem(FDVertex *startVertex)
{
	//IMPORTANT! the parameters are in normalization values. They must be converted !
	// the tunneling direction is north
	this->areaFactor = startVertex->EastLength / 2 + startVertex->WestLength / 2;
	FDVertex *currentVertex = startVertex;
	while (currentVertex != NULL)
	{
		this->deltaX.push_back((currentVertex->NorthLength + currentVertex->SouthLength) / 2);
		//the method to get physical property is changed
		//this->cbegde.push_back(currentVertex->Phys.conductionBandEnergy);
		//this->emass.push_back(currentVertex->Phys.electronMass);
		currentVertex = currentVertex->NorthVertex; // move the vertex north
	}
}
