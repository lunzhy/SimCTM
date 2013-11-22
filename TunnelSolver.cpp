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
	Normalization norm = Normalization();
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		if (cbEdge.at(ix) >= energy)
		{
			dX = norm.PullLength(deltaX.at(ix));
			energyDiff = norm.PullEnergy(cbEdge.at(ix) - energy);
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

	// the unit of calculated current density is [A/m^2]
	// Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	double Emin = cbedgeTunnelFrom > cbedgeTunnelTo ? cbedgeTunnelFrom : cbedgeTunnelTo; // in normalized value
	double Emax = cbEdge.front(); // in normalized value

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double hbar = SctmPhys::hbar;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double TC = 0;
	double supply = 0;

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		TC = getTransCoeff(currEnergy);
		supply = getSupplyFunction(currEnergy);
		//mass is normalized using m0, E is normalized using q
		DTFNdensity += 1/2 / pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
			* TC
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/cm^2]
	this->eCurrDens = DTFNdensity * per_m2_in_per_cm2;
	return eCurrDens;
}

double TunnelSolver::calcThermalEmission()
{
	double TEdensity = 0; // in [A/m^2]
	double Emin = cbEdge.front(); // in [eV]

	Normalization norm = Normalization();
	Emin = norm.PullEnergy(Emin); // in [eV]
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
		TEdensity += 1/2 / pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
			* 1  //for transmission coefficient
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/cm^2]
	this->eCurrDens = TEdensity * per_m2_in_per_cm2;
	return eCurrDens;
}

TunnelSolver::TunnelSolver(FDDomain *_domain): domain(_domain)
{
	this->eCurrDens = 0;
}

void TunnelSolver::SolveTunnel_Interface()
{
	eCurrDens_Interface.clear();
	//this->currentDensity = 0;
	//calcDTFNtunneling();
	//calcThermalEmission();
	//SctmUtils::UtilsDebug.PrintValue(this->currentDensity);
	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsTunnelStart.size(); ++iVert)
	{
		setSolver_Interface(vertsTunnelStart.at(iVert));
		solve_Interface();
	}
}

void TunnelSolver::ReadInput(VertexMapDouble fermi)
{
	int vertID = 0;
	double val = 0;
	FDVertex *vert = NULL;
	for (size_t iVert = 0; iVert != vertsTunnelStart.size(); ++iVert)
	{
		vert = vertsTunnelStart.at(iVert);
		vertID = vert->GetID();

		SCTM_ASSERT(fermi.find(vertID) != fermi.end(), 10020);
		val = fermi[vertID];
		fermiAboveMap[vertID] = val;
	}
}

double TunnelSolver::solve_Interface()
{
	double ret = 0;
	ret = calcDTFNtunneling();
	return ret;
}

void SubsToGateEletronTunnel::initialize()
{
	this->vertsTunnelStart.clear();
	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->domain->GetVertices().size(); ++iVert)
	{
		//the sequence in the vertices vector is its corresponding vertex index
		currVert = this->domain->GetVertex(iVert);
		if ( (currVert->IsAtContact()) && (currVert->Contact->ContactName == "Gate") )
		{
			vertsTunnelStart.push_back(currVert);
		}
	}
	//set the effective tunneling mass
	this->effTunnelMass = 0.5; // in m0
}

void SubsToGateEletronTunnel::setSolver_Interface(FDVertex *startVertex)
{
	//reset the vectors used in the calculation
	cbEdge.clear();
	//
	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// the tunneling direction is north
	using SctmUtils::Normalization;
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
			//in terms of the conduction band edge where the electron tunneling in, use the cbedge of the second
			//vertex instead, as an approximation
			vertsTunnelEnd_Interface.push_back(currVert);
			//set the conduction band edge of the interface between tunneling layer and trapping layer
			this->cbedgeTunnelTo = currVert->NorthVertex->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
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
		this->deltaX.push_back(dx);
		this->elecMass.push_back(emass);
		this->cbEdge.push_back(cbedge);

		currVert = currVert->NorthVertex; // to run the iteration
	}
	//set the silicon band edge, because the difference is fixed
	using namespace MaterialDB;
	double diff = 0;
	diff = GetMatPrpty(&MaterialDB::Silicon, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(&MaterialDB::SiO2, MatProperty::Mat_ElectronAffinity);
	cbedgeTunnelFrom = startVertex->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
	//set the fermi energy of the tunneling-in vertex
	fermiEnergyTunnelFrom = fermiAboveMap[startVertex->GetID()];
}

SubsToGateEletronTunnel::SubsToGateEletronTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	initialize();
}

double SubsToGateEletronTunnel::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 0; // assume the energy level in CB of trapping layer is always empty

	double energyDiff = 0; // in [eV], energy - fermi energy tunneling from
	Normalization norm = Normalization();
	energyDiff = norm.PullEnergy(energy - EfTunnelFrom);

	if (energy - EfTunnelFrom > 5 * kB * T / q)
		integralTunnelFrom =  kB * T * SctmMath::exp(- q * energyDiff / kB / T);
	else
		integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * energyDiff / kB / T));

	double supply = integralTunnelFrom - integralTunnelTo;
	return supply;
}
