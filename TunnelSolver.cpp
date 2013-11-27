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
#include "SctmUtils.h"
#include "SctmPhys.h"
#include "FDDomain.h"
#include "Normalization.h"
#include "Material.h"
#include <algorithm>

using SctmPhys::PhysProperty;
using SctmUtils::SctmFileStream;
using SctmUtils::Normalization;

double TunnelSolver::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;
	double EfTunnelTo = this->fermiEnergyTunnelTo;

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 1;

	if (energy - EfTunnelFrom > 5 * kB * T / q)
		integralTunnelFrom =  kB * T * SctmMath::exp(- q * (energy - EfTunnelFrom) / kB / T);
	else
		integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelFrom) / kB / T));

	if (energy - EfTunnelTo > 5 * kB * T / q)
		integralTunnelTo = kB * T * SctmMath::exp(- q * (energy - EfTunnelTo) / kB / T);
	else
		integralTunnelTo = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelTo) / kB / T));
	
	double supply = integralTunnelFrom - integralTunnelTo;
	return supply;
}

double TunnelSolver::getTransCoeff(double energy)
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
	double dX = 0;
	double energyDiff = 0; // V(x) - E, potential energy - current energy

	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		energyDiff = cbEdge.at(ix) - energy;
		if (energyDiff >= 0)
		{
			dX = deltaX.at(ix);
			integral += SctmMath::sqrt( 2 * elecMass.at(ix) * m0 * q * energyDiff ) 
						* dX * cm_in_m;
		}
	}
	double tc = tunnelFactor * SctmMath::exp( -2 / hbar * integral );

	return tc;
}

double TunnelSolver::calcDTFNtunneling()
{
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	//Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	//Emin and Emax are real values, in [eV]
	double Emin = cbedgeTunnelFrom > cbedgeTunnelTo ? cbedgeTunnelFrom : cbedgeTunnelTo; // in normalized value
	double Emax = cbEdge.front();

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double hbar = SctmPhys::hbar;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q

	double TC = 0;
	double supply = 0;

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		TC = getTransCoeff(currEnergy);
		supply = getSupplyFunction(currEnergy);
		//mass is normalized using m0, E is normalized using q
		DTFNdensity += 0.5 * pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
			* TC
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/m^2]
	this->eCurrDens = DTFNdensity;
	return eCurrDens;
}

double TunnelSolver::calcThermalEmission()
{
	double TEdensity = 0; // in [A/m^2]
	double Emin = cbEdge.front(); // in [eV]
	double Emax = Emin + 5; // only calculate the energy within 5eV large than the smallest energy to surpass the barrier

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double hbar = SctmPhys::hbar;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double supply = 0;
	double currEnergy = Emin;
	//mass is normalized using m0, E is normalized using q
	while (currEnergy <= Emax)
	{
		supply = getSupplyFunction(currEnergy);
		TEdensity += 0.5 / pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
			* 1  //for transmission coefficient
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/m^2]
	this->eCurrDens = TEdensity;
	return eCurrDens;
}

TunnelSolver::TunnelSolver(FDDomain *_domain): domain(_domain)
{
	this->eCurrDens = 0;
	this->temperature = 300;
}

void TunnelSolver::ReadInput(VertexMapDouble &fermi)
{
	int vertID = 0;
	double val = 0;
	FDVertex *vert = NULL;
	for (size_t iVert = 0; iVert != vertsStart_Tunnel.size(); ++iVert)
	{
		vert = vertsStart_Tunnel.at(iVert);
		vertID = vert->GetID();

		SCTM_ASSERT(fermi.find(vertID) != fermi.end(), 10020);
		val = fermi[vertID];
		fermiAboveMap[vertID] = val;
	}
}

double TunnelSolver::solveCurrDens_Tunnel()
{
	double ret = 0;
	ret = calcDTFNtunneling();
	ret += calcThermalEmission();
	return ret;
}


void SubsToTrapElecTunnel::initialize()
{
	this->vertsStart_Tunnel.clear();
	//this->vertsEnd_Trap.clear();

	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->domain->GetVertices().size(); ++iVert)
	{
		//the sequence in the vertices vector is its corresponding vertex index
		currVert = this->domain->GetVertex(iVert);
		if ( (currVert->IsAtContact()) && (currVert->Contact->ContactName == "Channel") )
		{
			vertsStart_Tunnel.push_back(currVert);
			//vertsEnd_Trap.push_back(currVert);
		}
	}
	//set the effective tunneling mass
	//TODO: this should be obtained from user input
	this->effTunnelMass = 0.5; // in m0
}

void SubsToTrapElecTunnel::setSolver_Tunnel(FDVertex *startVertex)
{
	//reset the vectors used in the calculation
	cbEdge.clear();
	elecMass.clear();
	deltaX.clear();
	//
	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// the tunneling direction is north
	Normalization norm = Normalization();
	double dx = 0;
	double emass = 0;
	double cbedge = 0;
	FDVertex *currVert = startVertex;
	
	while(true)
	{
		//for boundary vertex, the south length equals to 0
		//CAUTION: the following process is case-dependent.
		//                 x--------x---------x
		//            tunneling interface trapping
		//special treatment has been considered at the interface

		//finish the including process when the vertex is at the trapping boundary
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsEnd_Tunnel.push_back(currVert);
			break;
		}

		if (currVert->NorthVertex->IsAtBoundary(FDBoundary::eDensity))
		{
			dx = currVert->SouthLength / 2 + currVert->NorthLength;
		}
		else
		{
			dx = currVert->SouthLength / 2 + currVert->NorthLength / 2;
		}
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		
		dx = norm.PullLength(dx);
		cbedge = norm.PullEnergy(cbedge);

		this->deltaX.push_back(dx);
		this->elecMass.push_back(emass);
		this->cbEdge.push_back(cbedge);

		currVert = currVert->NorthVertex; // to run the iteration
	}
	//CAUTION this is a temporary method
	//set the silicon band edge, because the difference is fixed
	using namespace MaterialDB;
	double barrier = 0;
	barrier = GetMatPrpty(&MaterialDB::Silicon, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(&MaterialDB::SiO2, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelFrom = cbEdge.front() - barrier;

	//set the conduction band edge in the trapping layer where the tunneling ends.
	barrier = GetMatPrpty(&MaterialDB::Si3N4, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(&MaterialDB::SiO2, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelTo = cbEdge.back() - barrier;
	
	//set the fermi energy of the tunneling-in vertex
	double fermiAbove = norm.PullEnergy(fermiAboveMap[startVertex->GetID()]);
	fermiEnergyTunnelFrom = cbedgeTunnelFrom + fermiAbove;
}

SubsToTrapElecTunnel::SubsToTrapElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	initialize();
}

double SubsToTrapElecTunnel::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 0; // assume the energy level in CB of trapping layer is always empty

	double energyDiff = 0; // in [eV], energy - fermi energy tunneling from
	energyDiff = energy - EfTunnelFrom;

	if (energy - EfTunnelFrom > 5 * kB * T / q)
		integralTunnelFrom =  kB * T * SctmMath::exp(- q * energyDiff / kB / T);
	else
		integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * energyDiff / kB / T));

	double supply = integralTunnelFrom - integralTunnelTo;
	return supply;
}

void SubsToTrapElecTunnel::SolveTunnel()
{
	eCurrDens_Tunnel.resize(vertsStart_Tunnel.size());
	vertsEnd_Tunnel.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsStart_Tunnel.size(); ++iVert)
	{
		setSolver_Tunnel(vertsStart_Tunnel.at(iVert));
		currdens = solveCurrDens_Tunnel();
		eCurrDens_Tunnel.at(iVert) = currdens;
	}
}

void SubsToTrapElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double currDens = 0;
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	Normalization norm = Normalization();
	for (size_t iVert = 0; iVert != vertsEnd_Tunnel.size(); ++iVert)
	{
		currVert = vertsEnd_Tunnel.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel saves the current density, which is in A/m2
		currDens = eCurrDens_Tunnel.at(iVert) * per_m2_in_per_cm2;
		ret[vertID] = norm.PushCurrDens(currDens);
	}
}


TrapToGateElecTunnel::TrapToGateElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	initialize();
}

double TrapToGateElecTunnel::getSupplyFunction(double energy)
{
	//in this class, the supply function is the coefficient of the electron density
	//i.e. the real supply function = ret * eDensity
	double ret = 0;

	double T = this->temperature;
	double kB = SctmPhys::k0;
	double h = SctmPhys::h;
	double q = SctmPhys::q;
	double pi = SctmMath::PI;
	double m0 = SctmPhys::ElectronMass;

	// prefactor A in f(E) = A*exp(-E/KT), and A is a dimensionless value
	// here, the prefactor lacks n (density), the dimension is m^3
	static double prefactor = h * h * h / 4 / pi / (SctmMath::sqrt(pi)/2) /
								(kB*T) / SctmMath::sqrt(kB * T) / 
								SctmMath::sqrt(2*effTunnelMass*m0*effTunnelMass*m0*effTunnelMass*m0);

	
	//the supply function has a dimension of m^3 * J
	ret = prefactor * kB * T * SctmMath::exp(- q * (energy - cbedgeTunnelFrom) / kB / T);
	return ret;
}

void TrapToGateElecTunnel::initialize()
{
	vertsEnd_Tunnel.clear();

	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->domain->GetVertices().size(); ++iVert)
	{
		//the sequence in the vertices vector is its corresponding vertex index
		currVert = this->domain->GetVertex(iVert);
		if (currVert->IsAtContact() && currVert->Contact->ContactName == "Gate")
		{
			vertsEnd_Tunnel.push_back(currVert);
			//vertsEnd_Trap.push_back(currVert);
		}
	}
	//set the effective tunneling mass
	//TODO: this should be obtained from user input
	this->effTunnelMass = 0.5;
}

void TrapToGateElecTunnel::setSolver_Tunnel(FDVertex *endVertex)
{
	//reset the vectors used in the calculation
	cbEdge.clear();
	elecMass.clear();
	deltaX.clear();
	//
	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// the tunneling direction is north, but the vertex find direction is south
	Normalization norm = Normalization();
	double dx = 0;
	double emass = 0;
	double cbedge = 0;
	FDVertex *currVert = endVertex;

	while(true)
	{
		//for boundary vertex, the north length equals to 0
		//CAUTION: the following process is case-dependent.
		//                 x----------x----------x
		//              trapping  interface   blocking
		//special treatment has been considered at the interface

		//finish the including process when the vertex is at the trapping boundary
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsStart_Tunnel.push_back(currVert);
			break;
		}

		if (currVert->SouthVertex->IsAtBoundary(FDBoundary::eDensity))
		{
			dx = currVert->SouthLength + currVert->NorthLength / 2;
		}
		else
		{
			dx = currVert->SouthLength / 2 + currVert->NorthLength / 2;
		}
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);

		dx = norm.PullLength(dx);
		cbedge = norm.PullEnergy(cbedge);

		this->deltaX.push_back(dx);
		this->elecMass.push_back(emass);
		this->cbEdge.push_back(cbedge);

		currVert = currVert->SouthVertex; // to run the iteration
	}

	std::reverse(deltaX.begin(), deltaX.end());
	std::reverse(elecMass.begin(), elecMass.end());
	std::reverse(cbEdge.begin(), cbEdge.end());

	//CAUTION this is a temporary method
	using namespace MaterialDB;
	double barrier = 0;
	barrier = GetMatPrpty(&MaterialDB::Si3N4, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(&MaterialDB::SiO2, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelFrom = cbEdge.front() - barrier;

	//TODO: read this from simulation configuration
	//band energy of metal contact, equals to the applied voltage
	double gateAppliedVoltage = -16; // in [V]
	cbedgeTunnelTo = gateAppliedVoltage;
}

void TrapToGateElecTunnel::SolveTunnel()
{
	eCurrDens_Tunnel.resize(vertsEnd_Tunnel.size());
	vertsStart_Tunnel.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsEnd_Tunnel.size(); ++iVert)
	{
		setSolver_Tunnel(vertsEnd_Tunnel.at(iVert));
		currdens = solveCurrDens_Tunnel();
		eCurrDens_Tunnel.at(iVert) = currdens;
	}
}

void TrapToGateElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double tunCoeff = 0;
	double cm_in_m = SctmPhys::cm_in_m;

	Normalization norm = Normalization();
	for (size_t iVert = 0; iVert != vertsEnd_Tunnel.size(); ++iVert)
	{
		currVert = vertsStart_Tunnel.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel stores the coefficient to calculate current density
		//the calculated coefficient has a dimension of A*m, and should be converted to A*cm
		//value[A*cm] * eDensity[cm^-3] = current density[A/cm^2]
		tunCoeff = eCurrDens_Tunnel.at(iVert) / cm_in_m;
		ret[vertID] = norm.PushTunCoeff(tunCoeff);
	}
}
