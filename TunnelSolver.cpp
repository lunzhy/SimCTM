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

	double m0 = SctmPhys::ElectronMass;
	double q = SctmPhys::ElementaryCharge;
	double hbar = SctmPhys::hbar;
	double cm_in_m = SctmPhys::cm_in_m;
	double integral = 0;
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		if (cbegde.at(ix) >= energy)
		{
			integral += SctmMath::sqrt(2 * emass.at(ix) * m0 * q * (cbegde.at(ix) - energy)) 
						* deltaX.at(ix) * cm_in_m;
		}
	}
	double tc = tunnelFactor * SctmMath::exp(-2 / hbar * integral);

	return tc;
}

void TunnelSolver::calcDTFNtunneling()
{
	// the unit of calculated current density is [A/m^2]
	// Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	double Emin = cbedgeTunnelFrom; // in [eV]
	double Emax = cbegde.front(); // in [eV]

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double h = SctmPhys::h;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q

	double TC = 0;
	double supply = 0;

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		TC = getTransmissionCoefficient(currEnergy);
		supply = getSupplyFunction(currEnergy);
		//mass is normalized using m0, E is normalized using q
		DTFNdensity += 4 * pi * this->effTunnelMass * m0 * q / h / h / h
			* TC
			* supply
			* q * dE;
		currEnergy += dE;
	}
	this->currentDensity += DTFNdensity * SctmPhys::per_sqr_m_in_per_sqr_cm;
}

void TunnelSolver::calcThermalEmission()
{
	double TEdensity = 0; // in [A/m^2]
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
		TEdensity += 4 * pi * this->effTunnelMass * q / h / h / h
			* 1  //for transmission coefficient
			* getSupplyFunction(currEnergy)
			* dE;
		currEnergy += dE;
	}

	this->currentDensity += TEdensity;
}

TunnelSolver::TunnelSolver()
{
	this->currentDensity = 0;
}

double TunnelSolver::GetCurrentDensity()
{
	return this->currentDensity;
}

void TunnelSolver::SolveTunneling()
{
	calcDTFNtunneling();
}

void SubsToGateEletronTunnel::PrepareProblem(FDVertex *startVertex)
{
	// IMPORTANT! the parameters are in normalization values. They must be converted !
	// the tunneling direction is north
	double val = 0;
	this->areaFactor = startVertex->EastLength / 2 + startVertex->WestLength / 2;
	FDVertex *currentVertex = startVertex;
	while (currentVertex != NULL)
	{
		this->deltaX.push_back((currentVertex->NorthLength + currentVertex->SouthLength) / 2);
		
		//the method to get physical property is changed
		//this->cbegde.push_back(currentVertex->Phys.conductionBandEnergy);
		val = currentVertex->Phys.GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		this->cbegde.push_back(val);
		//this->emass.push_back(currentVertex->Phys.electronMass);
		val = currentVertex->Phys.GetPhysPrpty(PhysProperty::eMass);
		this->emass.push_back(val);
		
		currentVertex = currentVertex->NorthVertex; // move the vertex north
	}
}

void TunnelTest::PrepareProblem(FDVertex *startVertex)
{
	double oxideThickness = 3; // in [nm]
	int gridNumber = 100;
	this->temperature = 300;
	this->effTunnelMass = 0.5; // in m0

	double bandedgeDifference = 3.15; // in [eV], Silicon and oxide band edge difference
	double siliconBandegde = -0.166; // in [eV]
	double elecField = 3.27e7; // in [V/cm]
	double gateVoltage = 10; // in [V]

	this->cbedgeTunnelFrom = siliconBandegde;
	this->fermiEnergyTunnelFrom = 0;
	this->fermiEnergyTunnelTo = -gateVoltage;

	double deltax = 0;
	deltax = oxideThickness / gridNumber * SctmPhys::nm_in_cm;
	double currentX = 0;

	double val = 0;
	for (size_t iVert = 0; iVert != gridNumber + 1; ++iVert)
	{
		currentX = iVert * deltax;

		val = (iVert == 0 || iVert == gridNumber) ? deltax / 2 : deltax;
		this->deltaX.push_back(val);

		val = this->cbedgeTunnelFrom + bandedgeDifference - currentX * elecField;
		this->cbegde.push_back(val);

		this->emass.push_back(0.48);
	}
}
