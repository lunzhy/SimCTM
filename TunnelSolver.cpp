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
#include "DomainDetails.h"
#include "Normalization.h"
#include "Material.h"
#include "DDSolver.h"
#include <algorithm>

using SctmPhys::PhysProperty;
using SctmUtils::SctmFileStream;
using SctmUtils::Normalization;
using SctmUtils::SctmGlobalControl;

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

double TunnelSolver::getTransCoeff(double energy, vector<double> &deltax, vector<double> &emass, vector<double> &cbedge, int size /* = 0 */, int startindex /* = 0 */)
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

	if (size == 0)
	{
		size = deltax.size();
	}

	for (vector<double>::size_type ix = startindex; ix != size; ++ ix)
	{
		energyDiff = cbedge.at(ix) - energy;
		if (energyDiff >= 0)
		{
			dX = deltax.at(ix);
			integral += SctmMath::sqrt(2 * emass.at(ix) * m0 * q * energyDiff)
						* dX * cm_in_m;
		}
	}
	double tc = tunnelFactor * SctmMath::exp( -2 / hbar * integral );

	return tc;
}

double TunnelSolver::calcDTFNtunneling(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge, double cbedgeMax)
{
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	//Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	//Emin and Emax are real values, in [eV]
	//TODO: the choice of Emin determines the tunneling direction
	double Emin = cbedgeTunnelFrom > cbedgeTunnelTo ? cbedgeTunnelFrom : cbedgeTunnelTo; // in real value
	//double Emax = cbedge.front();
	double Emax = cbedgeMax;

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
		TC = getTransCoeff(currEnergy, deltaX, emass, cbedge);
		supply = getSupplyFunction(currEnergy);
		//mass is normalized using m0, E is normalized using q
		DTFNdensity += 0.5 / pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
			* TC
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/m^2]
	this->eCurrDens = DTFNdensity;
	return eCurrDens;
}

double TunnelSolver::calcThermalEmission(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge, double cbedgeMin)
{
	double TEdensity = 0; // in [A/m^2]
	//double Emin = cbedge.front(); // in [eV]
	double Emin = cbedgeMin;
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
	using SctmUtils::SctmGlobalControl;
	this->temperature = SctmGlobalControl::Get().Temperature;
	this->eCurrDens = 0;
	initialize();
}

void TunnelSolver::ReadInput(VertexMapDouble &fermi)
{
	int vertID = 0;
	double val = 0;
	FDVertex *vert = NULL;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		vert = vertsTunnelOxideStart.at(iVert);
		vertID = vert->GetID();

		SCTM_ASSERT(fermi.find(vertID) != fermi.end(), 10020);
		val = fermi[vertID];
		fermiAboveMap[vertID] = val;
	}
}


void TunnelSolver::loadBandStructure(FDVertex *startVert)
{
	//clear and reset all the vectors
	cbEdge_Tunnel.clear();
	eMass_Tunnel.clear();
	deltaX_Tunnel.clear();

	cbEdge_Trap.clear();
	eMass_Trap.clear();
	deltaX_Trap.clear();
	eEnergyLevel_Trap.clear();
	verts_Trap.clear();

	cbEdge_Block.clear();
	eMass_Block.clear();
	deltaX_Block.clear();

	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// CAUTION! the tunneling direction is north
	Normalization norm = Normalization(this->temperature);
	double dx = 0;
	double emass = 0;
	double cbedge = 0;
	double trapEnergyLevel = 0;
	double trapDepth = 0;
	FDVertex *currVert = startVert;

	using namespace MaterialDB;
	static Mat::Name tunnelMat = domain->GetRegion("Tunnel")->Mat->MatName();
	static Mat::Name trapMat = domain->GetTrapMatName();//trap region name is different in different cases
	static Mat::Name blockMat = domain->GetRegion("Block")->Mat->MatName();

	//for boundary vertex, the south length equals to 0
	//CAUTION: the following process is case-dependent. (south --> north)
	//                 x--------x---------x
	//            tunneling interface trapping
	//special treatment has been considered at the interface

	//set the tunneling oxide
	while (true)
	{
		//finish the including process when the vertex is at the trapping boundary
		//for boundary vertex of trapping layer, it always at the boundary of eDensity

		//load dx
		dx = 0;
		if (currVert->SouthVertex != NULL)
		{
			dx += currVert->SouthLength / 2;
		}
		if (currVert->Trap == NULL)
		{
			dx += currVert->NorthLength / 2;
		}
		dx = norm.PullLength(dx);

		//load cbedge
		if (currVert->Trap == NULL)
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		}
		else
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, tunnelMat);
		}
		cbedge = norm.PullEnergy(cbedge);

		//load emass
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);

		deltaX_Tunnel.push_back(dx);
		cbEdge_Tunnel.push_back(cbedge);
		eMass_Tunnel.push_back(emass);

		//check if current vertex is the ending vertex for tunneling layer
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsTunnelOxideEnd.push_back(currVert);
			break;
		}
		currVert = currVert->NorthVertex;
	}

	//set the trapping layer
	double dy = 0;
	while (currVert->Trap != NULL)
	{
		//load dx
		DriftDiffusionSolver::getDeltaXYAtVertex(currVert, dx, dy);
		dx = norm.PullLength(dy); //dy means the tunneling direction is south to north

		//load cbedge
		if (!currVert->Phys->HasMultiPrpty(PhysProperty::ConductionBandEnergy))
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		}
		else
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, trapMat);
		}
		cbedge = norm.PullEnergy(cbedge);

		//load emass
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);

		//load trap energy level
		trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::EnergyFromCondBand);
		trapEnergyLevel = cbedge - norm.PullEnergy(trapDepth);

		cbEdge_Trap.push_back(cbedge);
		deltaX_Trap.push_back(dx);
		eMass_Trap.push_back(emass);
		eEnergyLevel_Trap.push_back(trapEnergyLevel);
		verts_Trap.push_back(currVert);

		currVert = currVert->NorthVertex;
	}
	currVert = currVert->SouthVertex; // compensate the last statement in while loop

	//set the blocking oxide
	vertsBlockOxideStart.push_back(currVert);
	if (SctmGlobalControl::Get().LateralTunneling && SlopingTunnelTrapToGate::IsSlopingTunnel(currVert))
	{
		SlopingTunnelTrapToGate slopeTunnel = SlopingTunnelTrapToGate(this->domain, currVert);
		slopeTunnel.LoadBandStructureAlongPath(deltaX_Block, cbEdge_Block, eMass_Block);
		this->tunnelTrapToGateEnable = true;
		vertsBlockOxideEnd.push_back(slopeTunnel.GetGateVertex());
	}
	else
	{
		while (true)
		{
			//load dx
			dx = 0;
			if (currVert->Trap == NULL)
			{
				dx += currVert->SouthLength / 2;
			}
			if (currVert->NorthVertex != NULL)
			{
				dx += currVert->NorthLength / 2;
			}
			dx = norm.PullLength(dx);

			//load cbedge
			if (currVert->Trap == NULL)
			{
				cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
			}
			else
			{
				cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
			}
			cbedge = norm.PullEnergy(cbedge);

			//load emass
			emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);

			deltaX_Block.push_back(dx);
			cbEdge_Block.push_back(cbedge);
			eMass_Block.push_back(emass);

			if (currVert->IsAtContact()) //the gate contact
			{
				tunnelTrapToGateEnable = true;
				vertsBlockOxideEnd.push_back(currVert);
				break;
			}
			//the if below should not be reached
			if (currVert->NorthVertex == NULL) //the north boundary, within the isolation region
			{
				tunnelTrapToGateEnable = false;
				//to guarantee that the vertices vectors have same amount of member and the corresponding vertices can be obtained with the same index. 
				vertsBlockOxideEnd.push_back(NULL);
				break;
			}
			currVert = currVert->NorthVertex;
		}
	}
	
}

void TunnelSolver::initialize()
{
	//the tunneling start vertex is assigned here, because it only depends on the device structure and
	//maintains unchanged during the following calculation.
	this->vertsTunnelOxideStart.clear();
	this->vertsTunnelOxideEnd.clear();
	this->vertsBlockOxideStart.clear();
	this->vertsBlockOxideEnd.clear();
	//the order of vertices in these for vectors are in correspondence with each other.

	FDVertex *currVert = NULL;

	vector<FDVertex *> channelVerts = this->domain->GetContact("Channel")->GetContactVerts();
	for (size_t iVert = 0; iVert != channelVerts.size(); ++iVert)
	{
		currVert = channelVerts.at(iVert);
		//the start vertices of tunneling oxide (channel vertices) are needed to fill the other three vectors.
		vertsTunnelOxideStart.push_back(currVert);
	}
}

double TunnelSolver::supplyFunction_forCurrDens(double energy)
{
	//in this calculation, it is assumed that the band where tunneling ends is almost empty.
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;
	double EfTunnelTo = this->fermiEnergyTunnelTo;

	//protect the possible wrong assignment of the tunnel direction
	if (EfTunnelFrom <= EfTunnelTo)
	{
		return 0;
	}

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 0; // assume the energy level in CB of trapping layer is always empty

	double energyDiff = 0; // in [eV], energy - Fermi energy tunneling from
	energyDiff = energy - EfTunnelFrom;

	if (energy - EfTunnelFrom > 5 * kB * T / q)
		integralTunnelFrom = kB * T * SctmMath::exp(-q * energyDiff / kB / T);
	else
		integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(-q * energyDiff / kB / T));

	if (energy - EfTunnelTo > 5 * kB * T / q)
		integralTunnelTo = kB * T * SctmMath::exp(-q * (energy - EfTunnelTo) / kB / T);
	else
		integralTunnelTo = kB * T * SctmMath::ln(1 + SctmMath::exp(-q * (energy - EfTunnelTo) / kB / T));


	double supply = integralTunnelFrom - integralTunnelTo;
	//the supply function has a dimension of [J]
	return supply;
}

double TunnelSolver::supplyFunction_forTunCoeff(double energy)
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
	static double prefactor = h * h * h / 8 / pi / (SctmMath::sqrt(pi) / 2) /
		(kB*T) / SctmMath::sqrt(kB * T) /
		SctmMath::sqrt(2 * this->effTunnelMass*m0*this->effTunnelMass*m0*this->effTunnelMass*m0);


	//the supply function has a dimension of m^3 * J
	ret = prefactor * kB * T * SctmMath::exp(-q * (energy - cbedgeTunnelFrom) / kB / T);
	return ret;
}

void TunnelSolver::loadBandStructureForHoles(FDVertex* startVert)
{
	//clear and reset all vectors;
	cbEdge_Tunnel.clear();
	eMass_Tunnel.clear();
	deltaX_Tunnel.clear();

	//only consider DT/FN tunneling currently, so these containers are filled but not used.
	cbEdge_Trap.clear();
	eMass_Trap.clear();
	deltaX_Trap.clear();
	eEnergyLevel_Trap.clear();
	verts_Trap.clear();

	cbEdge_Block.clear();
	eMass_Block.clear();
	deltaX_Block.clear();

	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	Normalization norm = Normalization(this->temperature);
	double dx = 0;
	double emass = 0;
	double cbedge = 0;
	double trapEnergyLevel = 0;
	double trapDepth = 0;
	FDVertex *currVert = startVert;

	using namespace MaterialDB;
	static Mat::Name tunnelMat = domain->GetRegion("Tunnel")->Mat->MatName();
	static Mat::Name trapMat = domain->GetTrapMatName();//trap region name is different in different cases
	static Mat::Name blockMat = domain->GetRegion("Block")->Mat->MatName();

	//set the tunneling oxide
	//use valence band information to fill the electron tunneling parameters (class members)
	//by adopting an offset value, 0, i.e. the following relation is used:
	//(energy value for electron tunneling, as the parameter names indicates) = 0 - (energy value for hole tunneling)
	while (true)
	{
		//load dx
		dx = 0;
		if (currVert->SouthVertex != NULL)
		{
			dx += currVert->SouthLength / 2;
		}
		if (currVert->Trap == NULL)
		{
			dx += currVert->NorthLength / 2;
		}
		dx = norm.PullLength(dx);

		//load vbedge, valence band edge
		if (currVert->Trap == NULL)
		{
			cbedge = - currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy);
		}
		else
		{
			cbedge = - currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, tunnelMat);
		}
		cbedge = norm.PullEnergy(cbedge);

		//load emass
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::hMass);

		deltaX_Tunnel.push_back(dx);
		cbEdge_Tunnel.push_back(cbedge);
		eMass_Tunnel.push_back(emass);

		//check if current vertex is the ending vertex for tunneling layer
		//boundary for eDensity is also boundary for hDensity
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsTunnelOxideEnd.push_back(currVert);
			break;
		}
		currVert = currVert->NorthVertex;
	}

	//set the trapping layer
	double dy = 0;
	while (currVert->Trap != NULL)
	{
		//load dx
		DriftDiffusionSolver::getDeltaXYAtVertex(currVert, dx, dy);
		dx = norm.PullLength(dy); //dy means the tunneling direction is south to north

		//load cbedge
		if (!currVert->Phys->HasMultiPrpty(PhysProperty::ValenceBandEnergy))
		{
			cbedge = - currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy);
		}
		else
		{
			cbedge = - currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, trapMat);
		}
		cbedge = norm.PullEnergy(cbedge);

		//load emass
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::hMass);

		//load hole trap energy level
		trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::EnergyFromValeBand);
		trapEnergyLevel = cbedge - norm.PullEnergy(trapDepth); // cbedge is offset value for vbedge

		cbEdge_Trap.push_back(cbedge);
		deltaX_Trap.push_back(dx);
		eMass_Trap.push_back(emass);
		eEnergyLevel_Trap.push_back(trapEnergyLevel);
		verts_Trap.push_back(currVert);

		currVert = currVert->NorthVertex;
	}
	currVert = currVert->SouthVertex; // compensate the last statement in while loop

	/*
	vertsBlockOxideStart.push_back(currVert);
	if (SctmGlobalControl::Get().LateralTunneling && SlopingTunnelTrapToGate::IsSlopingTunnel(currVert))
	{
		SlopingTunnelTrapToGate slopeTunnel = SlopingTunnelTrapToGate(this->domain, currVert);
		slopeTunnel.LoadBandStructureAlongPath(deltaX_Block, cbEdge_Block, eMass_Block);
		this->tunnelTrapToGateEnable = true;
		vertsBlockOxideEnd.push_back(slopeTunnel.GetGateVertex());
	}
	else
	{
	*/
	while (true)
	{
		//load dx
		dx = 0;
		if (currVert->Trap == NULL)
		{
			dx += currVert->SouthLength / 2;
		}
		if (currVert->NorthVertex != NULL)
		{
			dx += currVert->NorthLength / 2;
		}
		dx = norm.PullLength(dx);

		//load cbedge
		if (currVert->Trap == NULL)
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		}
		else
		{
			cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
		}
		cbedge = norm.PullEnergy(cbedge);

		//load emass
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);

		deltaX_Block.push_back(dx);
		cbEdge_Block.push_back(cbedge);
		eMass_Block.push_back(emass);

		if (currVert->IsAtContact()) //the gate contact
		{
			tunnelTrapToGateEnable = true;
			vertsBlockOxideEnd.push_back(currVert);
			break;
		}
		//the if below should not be reached
		if (currVert->NorthVertex == NULL) //the north boundary, within the isolation region
		{
			tunnelTrapToGateEnable = false;
			//to guarantee that the vertices vectors have same amount of member and the corresponding vertices can be obtained with the same index. 
			vertsBlockOxideEnd.push_back(NULL);
			break;
		}
		currVert = currVert->NorthVertex;
	}
	//}
}


void SubsToTrapElecTunnel::setSolver_DTFN(FDVertex *startVertex)
{
	Normalization norm = Normalization(this->temperature);

	//set the silicon band edge, because the difference is fixed
	//all values of band edge are in real value, so pulling is needed.

	//set the fermi energy of the tunneling-in vertex
	//Pulling of the parameters is done here, because TunnelSolver uses real value internally
	double fermiAbove = norm.PullEnergy(fermiAboveMap[startVertex->GetID()]);

	if (eTunDirection == TunnelDirection::North)
	{
		cbedgeTunnelFrom = cbEdge_Tunnel.front() - this->subsBarrier;

		//set the conduction band edge in the trapping layer where the tunneling ends.
		cbedgeTunnelTo = cbEdge_Trap.front();

		fermiEnergyTunnelFrom = cbedgeTunnelFrom + fermiAbove;
	}
	else if (eTunDirection == TunnelDirection::South)
	{
		cbedgeTunnelFrom = cbEdge_Trap.front();
		cbedgeTunnelTo = cbEdge_Tunnel.front() - this->subsBarrier;
		//TODO: this is useless calculation tunneling coefficient.
		fermiEnergyTunnelFrom = 0;
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10041);
	}
}
	

SubsToTrapElecTunnel::SubsToTrapElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	using namespace MaterialDB;
	//substrate/tunnel oxide barrier
	Normalization norm = Normalization(this->temperature);
	this->subsBarrier = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElectronAffinity) 
						- GetMatPrpty(domain->GetRegion("Tunnel")->Mat, MatProperty::Mat_ElectronAffinity); //in normalized value
	this->subsBarrier = norm.PullPotential(subsBarrier);
}

double SubsToTrapElecTunnel::getSupplyFunction(double energy)
{
	if (eTunDirection == TunnelDirection::North)
	{
		return TunnelSolver::supplyFunction_forCurrDens(energy);
	}
	else if (eTunDirection == TunnelDirection::South)
	{
		return TunnelSolver::supplyFunction_forTunCoeff(energy);
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10041);
	}
	return 0;
}

void SubsToTrapElecTunnel::SolveTunnel()
{
	//the following containers are filled when loading band structure, so they are cleared at the beginning of
	//a new turn to solve tunnel problems.
	//Only vertsTunnelOxideStart exists.
	vertsTunnelOxideEnd.clear();
	vertsBlockOxideStart.clear();
	vertsBlockOxideEnd.clear();

	//clear the results of last calculation
	eCurrDens_DTFN.clear();
	eCurrDens_DTFN.resize(vertsTunnelOxideStart.size());
	eCurrDensMap_MFN.clear();
	eCurrDensMap_B2T.clear();

	double cbedge_max = 0;
	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		setTunnelDirection(vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
		setTunnelTag();

		//setSolver_DTFN(vertsTunnelOxideStart.at(iVert));

		if (eTunDirection == TunnelDirection::North) // tunneling in
		{
			cbedge_max = cbEdge_Tunnel.front();
		}
		else if (eTunDirection == TunnelDirection::South) // tunneling out
		{
			cbedge_max = cbEdge_Tunnel.back();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

		currdens = calcDTFNtunneling(deltaX_Tunnel, eMass_Tunnel, cbEdge_Tunnel, cbedge_max);
		currdens += calcThermalEmission(deltaX_Tunnel, eMass_Tunnel, cbEdge_Tunnel, cbedge_max);
		eCurrDens_DTFN.at(iVert) = currdens;

		setSolver_Trap();

		if (SctmGlobalControl::Get().PhysicsMFN)
		{
			calcCurrDens_MFN();
		}
		
		if (SctmGlobalControl::Get().PhysicsB2T)
		{
			calcCurrDens_B2T();
		}

		if (SctmGlobalControl::Get().PhysicsT2B)
		{
			calcTransCoeff_T2B();
		}
	}
}

void SubsToTrapElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	//Result returning does not care the direction of the tunneling current.
	//So the sign of the result is not determined here, it is set in solver pack.
	FDVertex *currVert = NULL;
	int vertID = 0;
	//for current density
	double currDens = 0;
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;
	//for tunneling coefficient
	double coeff = 0;
	double cm_in_m = SctmPhys::cm_in_m;

	Normalization norm = Normalization(this->temperature);
	FDBoundary::TunnelTag tunnelTag = FDBoundary::noTunnel;

	for (size_t iVert = 0; iVert != vertsTunnelOxideEnd.size(); ++iVert)
	{
		//currVert = vertsEnd.at(iVert);
		currVert = vertsTunnelOxideEnd.at(iVert);
		vertID = currVert->GetID();
		
		tunnelTag = currVert->BndCond.GetBCTunnelTag();
		switch (tunnelTag)
		{
			case FDBoundary::eTunnelOut:
				eTunDirection = TunnelSolver::South;
				break;
			case FDBoundary::eTunnelIn:
				eTunDirection = TunnelSolver::North;
				break;
			case FDBoundary::hTunnelIn:
			case FDBoundary::hTunnelOut:
			case FDBoundary::noTunnel:
			default:
				SCTM_ASSERT(SCTM_ERROR, 10056);
				break;
		}

		if (eTunDirection == TunnelDirection::North)
		{
			//currently the eCurrDens_Tunnel saves the current density, which is in A/m2
			currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
			ret[vertID] = norm.PushCurrDens(currDens);
		}
		else if (eTunDirection == TunnelDirection::South)
		{
			//currently the eCurrDens_DTFN saves the coefficient to calculate current density,
			//which has a dimension of A*m.
			coeff = eCurrDens_DTFN.at(iVert) / cm_in_m;
			ret[vertID] = norm.PushTunCoeff(coeff);
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}
		
	}
}

void SubsToTrapElecTunnel::setSolver_Trap()
{
	//reset the vectors
	cbEdge_TunnelTrap.clear();
	eMass_TunnelTrap.clear();
	deltaX_TunnelTrap.clear();
	//
	cbEdge_TunnelTrap.reserve(this->cbEdge_Tunnel.size() + this->cbEdge_Trap.size());
	cbEdge_TunnelTrap.insert(cbEdge_TunnelTrap.end(), this->cbEdge_Tunnel.begin(), this->cbEdge_Tunnel.end());
	cbEdge_TunnelTrap.insert(cbEdge_TunnelTrap.end(), this->cbEdge_Trap.begin(), this->cbEdge_Trap.end());

	deltaX_TunnelTrap.reserve(this->deltaX_Tunnel.size() + this->deltaX_Trap.size());
	deltaX_TunnelTrap.insert(deltaX_TunnelTrap.end(), this->deltaX_Tunnel.begin(), this->deltaX_Tunnel.end());
	deltaX_TunnelTrap.insert(deltaX_TunnelTrap.end(), this->deltaX_Trap.begin(), this->deltaX_Trap.end());

	eMass_TunnelTrap.reserve(this->eMass_Tunnel.size() + this->eMass_Trap.size());
	eMass_TunnelTrap.insert(eMass_TunnelTrap.end(), this->eMass_Tunnel.begin(), this->eMass_Tunnel.end());
	eMass_TunnelTrap.insert(eMass_TunnelTrap.end(), this->eMass_Trap.begin(), this->eMass_Trap.end());
}

void SubsToTrapElecTunnel::calcCurrDens_MFN()
{
	//TODO: write debug information for MFN
	if (cbedgeTunnelFrom > cbedgeTunnelTo)
	{
		//no Modified Fowler-Nordheim tunneling happens
		return;
	}

	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	double eCurrDens = 0; // in [A/m^2]

	//Emin and Emax are processed with q. (divided by q)
	//Emin and Emax are real values, in [eV]
	double Emin = cbedgeTunnelFrom;
	double Emax = cbedgeTunnelTo;

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double hbar = SctmPhys::hbar;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double TC = 0; // the tunneling coefficient
	double supply = 0; // the value of supply function

	double currEnergy = Emax;
	int vertID = 0;
	int vSize = 0;
	FDVertex *vertToAssign = NULL;

	while (currEnergy >= Emin)
	{
		//calculate MFN tunneling current
		vertToAssign = findTrapVertex_MFN(currEnergy, vSize);
		if (vertToAssign != NULL)
		{
			TC = getTransCoeff(currEnergy, deltaX_TunnelTrap, eMass_TunnelTrap, cbEdge_TunnelTrap, vSize);
			supply = getSupplyFunction(currEnergy);
			//mass is normalized using m0, E is normalized using q
			eCurrDens = 0.5 / pi / pi * this->effTunnelMass * m0 * q / hbar / hbar / hbar
				* TC
				* supply
				* q * dE;

			//[A/m^2] -> [A/cm^2]
			eCurrDens = eCurrDens * per_m2_in_per_cm2;

			//assign the calculated current density to this vertex is done outside this class in SolverPack.
			vertID = vertToAssign->GetID();
			eCurrDensMap_MFN[vertID] += eCurrDens;
		}

		//calculate the tunneling current of each energy level from max to min
		currEnergy -= dE;
	}
}

FDVertex * SubsToTrapElecTunnel::findTrapVertex_MFN(double energy, int &size)
{
	//assign the tunneling to the vertex with lower band edge
	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		if (cbEdge_Trap.at(iVert) <= energy)
		{
			//size = index + 1
			size = cbEdge_Tunnel.size() + iVert + 1;
			return verts_Trap.at(iVert);
		}
	}
	return NULL;
}

void SubsToTrapElecTunnel::ReturnResult_MFN(VertexMapDouble &ret)
{
	double eCurrDens_in_per_cm2 = 0;
	Normalization norm = Normalization(this->temperature);

	for (VertexMapDouble::iterator it = eCurrDensMap_MFN.begin(); it != eCurrDensMap_MFN.end(); ++it)
	{
		eCurrDens_in_per_cm2 = it->second;
		ret[it->first] = norm.PushCurrDens(eCurrDens_in_per_cm2);
	}
}

FDVertex * SubsToTrapElecTunnel::findTrapVertex_B2T(double energy, int &size)
{
	//assign the tunneling to vertex with higher trap energy
	for (size_t iVert = 0; iVert != verts_Trap.size() - 1; ++iVert)
	{
		if ((eEnergyLevel_Trap.at(iVert) >= energy) && (eEnergyLevel_Trap.at(iVert + 1) < energy))
		{
			//size = vertices size of tunneling oxide + index + 1
			size = cbEdge_Tunnel.size() + iVert + 1;
			return verts_Trap.at(iVert);
		}
	}
	return NULL;
}

void SubsToTrapElecTunnel::ReturnResult_B2T(VertexMapDouble &ret)
{
	double eCurrDens_in_per_cm2 = 0;
	Normalization norm = Normalization(this->temperature);
	for (VertexMapDouble::iterator it = eCurrDensMap_B2T.begin(); it != eCurrDensMap_B2T.end(); ++it)
	{
		eCurrDens_in_per_cm2 = it->second;
		ret[it->first] = norm.PushCurrDens(eCurrDens_in_per_cm2);
	}
}

void SubsToTrapElecTunnel::calcCurrDens_B2T()
{
	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	double eCurrDens = 0; // in [A/m^2]

	//Emin and Emax are processed with q. (divided by q)
	//Emin and Emax are real values, in [eV]
	double Emin = cbedgeTunnelFrom;
	double Emax = eEnergyLevel_Trap.front();

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double hbar = SctmPhys::hbar;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double TC = 0; // the tunneling coefficient
	double supply = 0; // the value of supply function

	double currEnergy = Emax;
	int vertID = 0;
	int vSize = 0;
	FDVertex *vertToAssign = NULL;

	while (currEnergy >= Emin)
	{
		//calculate trap-to-band tunneling
		vertToAssign = findTrapVertex_B2T(currEnergy, vSize);
		if (vertToAssign != NULL)
		{
			TC = getTransCoeff(currEnergy, deltaX_TunnelTrap, eMass_TunnelTrap, cbEdge_TunnelTrap, vSize);
			supply = getSupplyFunction(currEnergy);
			eCurrDens = q * 0.5 / pi / pi * this->effTunnelMass * m0 / hbar / hbar / hbar
				* TC
				* supply
				* q * dE;

			//[A/m^2] -> [A/cm^2]
			eCurrDens = eCurrDens * per_m2_in_per_cm2;
			vertID = vertToAssign->GetID();
			eCurrDensMap_B2T[vertID] += eCurrDens;
		}

		//calculate the tunneling current of each energy level from max to min
		currEnergy -= dE;
	}
}

void SubsToTrapElecTunnel::setTunnelTag()
{
	using namespace SctmUtils;
	FDVertex *verts = verts_Trap.front(); // the front element is vertex at tunnelOxide/trappingLayer interface
	switch (eTunDirection)
	{
		case TunnelDirection::North:
		{
			verts->BndCond.SetTunnelTag(FDBoundary::eTunnelIn);
			break;
		}
		case TunnelDirection::South:
		{
			verts->BndCond.SetTunnelTag(FDBoundary::eTunnelOut);
			break;
		}
		case TunnelDirection::East:
		case TunnelDirection::West:
		default:
			SCTM_ASSERT(SCTM_ERROR, 10041);
			break;
	}

	/*
	double elecField = 0;
	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != vertsTunnelOxideEnd.size(); ++iVert)
	{
		currVert = vertsTunnelOxideEnd.at(iVert);
		//TODO: this is a temporary method to set the tunnel tags for vertices at interface.
		elecField = currVert->Phys->GetPhysPrpty(PhysProperty::ElectricField_Y);
		if (elecField < 0)
		{
			currVert->BndCond.SetTunnelTag(FDBoundary::eTunnelIn);
		}
	}
	*/
}

void SubsToTrapElecTunnel::setTunnelDirection(FDVertex *vertSubs, FDVertex *vertTrap)
{
	using namespace MaterialDB;

	Normalization norm = Normalization(this->temperature);
	double per_cm3_in_per_m3 = SctmPhys::per_cm3_in_per_m3;
	double m0 = SctmPhys::m0;
	double k0 = SctmPhys::k0;
	double h = SctmPhys::h;
	double pi = SctmMath::PI;
	double q = SctmPhys::q;

	//two fermi energy are in real value, in [eV]
	double fermiSubs = 0;
	double fermiTrap = 0; //the electron quasi-fermi energy of trapping layer
	
	int vertIdSubs = 0;
	vertIdSubs = vertSubs->GetID();
	SCTM_ASSERT(fermiAboveMap.find(vertIdSubs) != fermiAboveMap.end(), 10051);
	fermiSubs = cbEdge_Tunnel.front() - this->subsBarrier + norm.PullPotential(fermiAboveMap[vertIdSubs]);

	double density = 0;
	density = vertTrap->Phys->GetPhysPrpty(PhysProperty::eDensity);
	density = norm.PullDensity(density); //in [cm^-3] up to now
	density = density * per_cm3_in_per_m3; //in [m^3] up to now
	
	double eMassTrap = SctmPhys::m0 * GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_ElecDOSMass);
	double fermiRel = 0;
	//use bandgap when density equals to 0, indicating that fermi energy is far below conduction band
	fermiRel = (density != 0) ? k0 * this->temperature / q * SctmMath::ln(density * h * h * h / 8 / pi / (SctmMath::sqrt(pi) / 2) /
		SctmMath::pow(k0 * this->temperature, 1.5) / SctmMath::sqrt(2 * eMassTrap * eMassTrap * eMassTrap)) : 
		-norm.PullPotential(GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_Bandgap));
	fermiTrap = fermiRel + cbEdge_Trap.front();

	double effMass = 0;
	if (density == 0 || fermiSubs >= fermiTrap) //tunneling into the trap layer
	{
		eTunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = fermiSubs;
		fermiEnergyTunnelTo = fermiTrap;

		cbedgeTunnelFrom = cbEdge_Tunnel.front() - this->subsBarrier;
		cbedgeTunnelTo = cbEdge_Trap.front();

		//set the effective DOS mass for tunneling
		//the material for substrate is silicon
		effMass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass);
		this->effTunnelMass = effMass; // the effective electron mass, in [m0]
	}
	else
	{
		eTunDirection = TunnelDirection::South;
		fermiEnergyTunnelFrom = fermiTrap; //however, this is useless in the calculation of tunneling-out coefficient
		fermiEnergyTunnelTo = fermiSubs;

		cbedgeTunnelFrom = cbEdge_Trap.front();
		cbedgeTunnelTo = cbEdge_Tunnel.front() - this->subsBarrier;

		//set the effective DOS mass for tunneling
		effMass = GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_ElecDOSMass);
		this->effTunnelMass = effMass;
	}

	/* previous method to judge the tunneling direction
	//CAUTION!! this is currently structure-dependent
	double elecField = 0;
	FDVertex *verts = verts_Trap.front(); // the front element is vertex at tunnelOxide/trappingLayer interface
	elecField = verts->Phys->GetPhysPrpty(PhysProperty::ElectricField_Y);
	if (elecField < 0)
	{
		//for electrons, this means program situation.
		eTunDirection = TunnelDirection::North;
	}
	else
	{
		//for electrons, this means retention situation.
		eTunDirection = TunnelDirection::South;
	}
	*/
}

void SubsToTrapElecTunnel::calcTransCoeff_T2B()
{
	double cbedgeTrap_min = cbEdge_Trap.front();
	double cbedgeSubs = cbedgeTunnelTo;
	double trapEnergyLevel = 0;
	double TC = 0;
	int vertID = 0;
	int size = 0;
	FDVertex *currVert = NULL;

	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		trapEnergyLevel = eEnergyLevel_Trap.at(iVert);
		if (trapEnergyLevel < cbedgeTrap_min && trapEnergyLevel > cbedgeSubs)
		{
			//size = vertices size of tunneling oxide + index + 1
			size = cbEdge_Tunnel.size() + iVert + 1;
			currVert = verts_Trap.at(iVert);
			TC = getTransCoeff(trapEnergyLevel, deltaX_TunnelTrap, eMass_TunnelTrap, cbEdge_TunnelTrap, size);

			vertID = currVert->GetID();
			eTransCoeffMap_T2B[vertID] = TC;
		}
	}
}

void SubsToTrapElecTunnel::ReturnResult_T2B(VertexMapDouble &ret)
{
	for (VertexMapDouble::iterator it = eTransCoeffMap_T2B.begin(); it != eTransCoeffMap_T2B.end(); ++it)
	{
		ret[it->first] = it->second;
	}
}


TrapToGateElecTunnel::TrapToGateElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	//the effective DOS mass for tunneling is set when determine the tunneling direction
}

double TrapToGateElecTunnel::getSupplyFunction(double energy)
{
	if (eTunDirection == TunnelDirection::North)
	{
		return TunnelSolver::supplyFunction_forTunCoeff(energy);
	}
	else if (eTunDirection == TunnelDirection::South)
	{
		return TunnelSolver::supplyFunction_forCurrDens(energy);
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10041);
	}
	return 0;

	/* previous method
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
	*/
}

void TrapToGateElecTunnel::setSolver_DTFN(FDVertex *endVertex)
{
	//TODO: this should be enhanced.
	//CAUTION: this method may fail at some situations.
	cbedgeTunnelFrom = cbEdge_Trap.back();

	//band energy of metal contact, equals to minus applied voltage.
	Normalization norm = Normalization(this->temperature);
	SCTM_ASSERT(endVertex->Contact != NULL, 10050);
	double gateVoltage = norm.PullPotential(endVertex->Contact->Voltage);
	cbedgeTunnelTo = -gateVoltage;
}

void TrapToGateElecTunnel::SolveTunnel()
{
	//the following containers are filled when loading band structure, so they are cleared at the beginning of
	//a new turn to solve tunnel problems.
	//Only vertsTunnelOxideStart exists.
	vertsTunnelOxideEnd.clear();
	vertsBlockOxideStart.clear();
	vertsBlockOxideEnd.clear();

	//clear the results of last calculation
	eCurrDens_DTFN.clear();
	eCurrDens_DTFN.resize(vertsTunnelOxideStart.size());
	eTransCoeffMap_T2B.clear();

	double currdens = 0;
	double cbedge_max = 0;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		if (!this->tunnelTrapToGateEnable)
		{
			continue;
		}

		setTunnelDirection(vertsBlockOxideStart.at(iVert), vertsBlockOxideEnd.at(iVert));
		setTunnelTag();

		if (eTunDirection == TunnelDirection::North) // tunneling out
		{
			cbedge_max = cbEdge_Block.front();
		}
		else if (eTunDirection == TunnelDirection::South) // tunneling in
		{
			cbedge_max = cbEdge_Block.back();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

		currdens = calcDTFNtunneling(deltaX_Block, eMass_Block, cbEdge_Block, cbedge_max);
		currdens += calcThermalEmission(deltaX_Block, eMass_Block, cbEdge_Block, cbedge_max);
		eCurrDens_DTFN.at(iVert) = currdens;

		if (SctmGlobalControl::Get().PhysicsT2B)
		{
			setSolver_Trap();
			calcTransCoeff_T2B();
		}
	}
}

void TrapToGateElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	//for current density
	double currDens = 0;
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;
	//for tunneling coefficient
	double tunCoeff = 0;
	double cm_in_m = SctmPhys::cm_in_m;

	Normalization norm = Normalization(this->temperature);
	FDBoundary::TunnelTag tunnelTag = FDBoundary::noTunnel;

	for (size_t iVert = 0; iVert != vertsBlockOxideStart.size(); ++iVert)
	{
		if (vertsBlockOxideEnd.at(iVert) == NULL)
		{
			//no contact vertex is at the end of the tunneling path, so no tunneling-out occurs.
			//not include this vertex in the returned result
			continue;
		}

		currVert = vertsBlockOxideStart.at(iVert);
		vertID = currVert->GetID();

		tunnelTag = currVert->BndCond.GetBCTunnelTag();
		switch (tunnelTag)
		{
			case FDBoundary::eTunnelOut:
			{
				eTunDirection = TunnelSolver::North;
				break;
			}
			case FDBoundary::eTunnelIn:
			{
				eTunDirection = TunnelSolver::South;
				break;
			}
			case FDBoundary::hTunnelIn:
			case FDBoundary::hTunnelOut:
			case FDBoundary::noTunnel:
			default:
				SCTM_ASSERT(SCTM_ERROR, 10056);
				break;
		}

		if (eTunDirection == TunnelDirection::North)
		{
			//currently the eCurrDens_DTFN stores the coefficient to calculate current density
			//the calculated coefficient has a dimension of A*m, and should be converted to A*cm
			//this value[A*cm] * eDensity[cm^-3] = current density[A/cm^2]
			tunCoeff = eCurrDens_DTFN.at(iVert) / cm_in_m;
			ret[vertID] = norm.PushTunCoeff(tunCoeff);
		}
		else if (eTunDirection == TunnelDirection::South)
		{
			currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
			ret[vertID] = norm.PushCurrDens(currDens);
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}
	}
}

void TrapToGateElecTunnel::setTunnelTag()
{
	using namespace SctmUtils;
	FDVertex *verts = verts_Trap.back(); //the back element is vertex at trappingLayer/blockingLayer interface
	switch (eTunDirection)
	{
		case TunnelDirection::North:
		{
			verts->BndCond.SetTunnelTag(FDBoundary::eTunnelOut);
			break;
		}
		case TunnelDirection::South:
		{
			verts->BndCond.SetTunnelTag(FDBoundary::eTunnelIn);
			break;
		}
		case TunnelDirection::East:
		case TunnelDirection::West:
		default:
		SCTM_ASSERT(SCTM_ERROR, 10041);
		break;
	}
}


void TrapToGateElecTunnel::setTunnelDirection(FDVertex *vertTrap, FDVertex *vertGate)
{
	using namespace MaterialDB;

	Normalization norm = Normalization(this->temperature);
	double per_cm3_in_per_m3 = SctmPhys::per_cm3_in_per_m3;
	double m0 = SctmPhys::m0;
	double k0 = SctmPhys::k0;
	double h = SctmPhys::h;
	double pi = SctmMath::PI;
	double q = SctmPhys::q;

	//two fermi energy are in real value, in [eV]
	double fermiTrap = 0; //the electron quasi-fermi energy of trapping layer
	double fermiGate = 0;


	SCTM_ASSERT(vertGate->Contact != NULL, 10050);
	double gateVoltage = norm.PullPotential(vertGate->Contact->Voltage);
	fermiGate = -gateVoltage;

	double density = 0;
	density = vertTrap->Phys->GetPhysPrpty(PhysProperty::eDensity);
	density = norm.PullDensity(density); //in [cm^-3]
	density = density * per_cm3_in_per_m3; //in [m^-3]

	double eMassTrap = SctmPhys::m0 * GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_ElecDOSMass);
	double fermiRel = 0;
	//use bandgap when density equals to 0, indicating that fermi energy is far below conduction band
	fermiRel = (density != 0) ? k0 * this->temperature / q * SctmMath::ln(density * h * h * h / 8 / pi / (SctmMath::sqrt(pi) / 2) /
		SctmMath::pow(k0 * this->temperature, 1.5) / SctmMath::sqrt(2 * eMassTrap * eMassTrap * eMassTrap)) : 
		-norm.PullPotential(GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_Bandgap));
	fermiTrap = fermiRel + cbEdge_Trap.back(); //the last vertex of trapping layer

	double effMass = 0;
	if (density == 0 || fermiGate > fermiTrap) //tunneling into the trap layer from gate
	{
		eTunDirection = TunnelDirection::South;
		fermiEnergyTunnelFrom = fermiGate;
		fermiEnergyTunnelTo = fermiTrap;

		cbedgeTunnelFrom = fermiGate - 10; //use a very large number (10eV) to indicate that gate has no edge
		cbedgeTunnelTo = cbEdge_Trap.back();

		//set the effective mass for tunneling current
		this->effTunnelMass = 1; //the effective mass in gate is 1
	}
	else //tunneling out of the trap layer
	{
		eTunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = fermiTrap;
		fermiEnergyTunnelTo = fermiGate - 10; //use a very large number (10eV) to indicate that gate has no edge

		cbedgeTunnelFrom = cbEdge_Trap.back();
		cbedgeTunnelTo = fermiGate; //for the gate, cbedgeTunnelTo is the same with the fermi energy of gate

		//set the effective DOS mass for tunneling
		effMass = GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_ElecDOSMass);
		this->effTunnelMass = effMass;
	}

	/* previous method to determine the tunneling direction
	double elecField = 0;
	FDVertex *currVert = NULL;

	currVert = verts_Trap.back(); //the back element is vertex at trappingLayer/blockingLayer interface
	//TODO: this is a temporary method to set the tunnel tags for vertices at interface.
	elecField = currVert->Phys->GetPhysPrpty(PhysProperty::ElectricField_Y);
	if (elecField < 0)
	{
		//for electrons, this means program and retention situation
		eTunDirection = TunnelSolver::North;
	}
	else
	{
		//for electrons, this means erase situation.
		eTunDirection = TunnelDirection::South;
	}
	*/
}


void TrapToGateElecTunnel::calcTransCoeff_T2B()
{
	double cbedge_min = cbEdge_Trap.back();
	double gateEdge = cbedgeTunnelTo;
	double trapEnergyLevel = 0;
	double TC = 0;
	int vertID = 0;
	int startindex = 0;
	FDVertex *currVert = NULL;

	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		trapEnergyLevel = eEnergyLevel_Trap.at(iVert);
		if (trapEnergyLevel < cbedge_min && trapEnergyLevel > gateEdge)
		{
			startindex = iVert;
			currVert = verts_Trap.at(iVert);
			//the index is same for container of Trap vertices and TrapBlock vertices.
			TC = getTransCoeff(trapEnergyLevel, deltaX_TrapBlock, eMass_TrapBlock, cbEdge_TrapBlock, 0, startindex);

			vertID = currVert->GetID();
			eTransCoeffMap_T2B[vertID] = TC;
		}
	}
}

void TrapToGateElecTunnel::setSolver_Trap()
{
	//clear the containers
	cbEdge_TrapBlock.clear();
	eMass_TrapBlock.clear();
	deltaX_TrapBlock.clear();

	cbEdge_TrapBlock.reserve(cbEdge_Trap.size() + cbEdge_Block.size());
	cbEdge_TrapBlock.insert(cbEdge_TrapBlock.end(), cbEdge_Trap.begin(), cbEdge_Trap.end());
	cbEdge_TrapBlock.insert(cbEdge_TrapBlock.end(), cbEdge_Block.begin(), cbEdge_Block.end());

	eMass_TrapBlock.reserve(eMass_Trap.size() + eMass_Block.size());
	eMass_TrapBlock.insert(eMass_TrapBlock.end(), eMass_Trap.begin(), eMass_Trap.end());
	eMass_TrapBlock.insert(eMass_TrapBlock.end(), eMass_Block.begin(), eMass_Block.end());

	deltaX_TrapBlock.reserve(deltaX_Trap.size() + deltaX_Block.size());
	deltaX_TrapBlock.insert(deltaX_TrapBlock.end(), deltaX_Trap.begin(), deltaX_Trap.end());
	deltaX_TrapBlock.insert(deltaX_TrapBlock.end(), deltaX_Block.begin(), deltaX_Block.end());
}

void TrapToGateElecTunnel::ReturnResult_T2B(VertexMapDouble &ret)
{
	for (VertexMapDouble::iterator it = eTransCoeffMap_T2B.begin(); it != eTransCoeffMap_T2B.end(); ++it)
	{
		ret[it->first] = it->second;
	}
}



bool SlopingTunnelTrapToGate::IsSlopingTunnel(FDVertex *vert)
{
	//the vertex is at the block/trap interface
	FDElement *southEastElem = vert->SoutheastElem;
	FDElement *southWestElem = vert->SouthwestElem;


	if (southEastElem != NULL && southWestElem != NULL && 
		southEastElem->Region->RegName == "Trap.Iso2" && southWestElem->Region->RegName == "Trap.Iso2")
	{
		return true;
	}
	
	
	if (southEastElem != NULL && southWestElem != NULL && 
		southEastElem->Region->RegName == "Trap.Iso3" && southWestElem->Region->RegName == "Trap.Iso3")
	{
		return true;
	}
	
	return false;
}

SlopingTunnelTrapToGate::SlopingTunnelTrapToGate(FDDomain *_domain, FDVertex *_vert)
{
	domain = _domain;
	vertStart = _vert;
	vertRegName = vertStart->SoutheastElem->Region->RegName;
	temperature = SctmGlobalControl::Get().Temperature;

	setBoundaryVerts();
	setInterpolatedValues();
}

void SlopingTunnelTrapToGate::LoadBandStructureAlongPath(vector<double> &dx, vector<double> &cbedge, vector<double> &emass)
{
	dx.push_back(this->dSlope.front() / 2);
	cbedge.push_back(this->cbEdge.front());
	emass.push_back(this->elecMass.front());

	for (size_t is = 0; is != this->dSlope.size() - 1; ++is)
	{
		dx.push_back((this->dSlope.at(is) + this->dSlope.at(is + 1)) / 2);
		cbedge.push_back(this->cbEdge.at(is + 1));
		emass.push_back(this->elecMass.at(is + 1));
	}

	dx.push_back(this->dSlope.back() / 2);
	cbedge.push_back(this->cbEdge.back());
	emass.push_back(this->elecMass.back());
	
}

void SlopingTunnelTrapToGate::setBoundaryVerts()
{
	//set the angle
	double xDist = 0;
	double yDist = 0;

	FDVertex *vert = NULL;

	if (vertRegName == "Trap.Iso2")
	{
		//the vertex is in left trap region under isolation 2
		vert = vertStart;
		while (vert->NortheastElem->Region->RegName != "Iso2")
		{
			vert = vert->NorthVertex;
		}
		upmostVert = vert;
		yDist = SctmMath::abs(vert->Y - vertStart->Y);

		vert = vertStart;
		while (vert->SoutheastElem->Region->RegName != "Trap.Gate2")
		{
			vert = vert->EastVertex;
		}
		easternmostVert = vert;
		xDist = SctmMath::abs(vert->X - vertStart->X);
		
		tanAngle = xDist / yDist;
	}
	else if (vertRegName == "Trap.Iso3")
	{
		//the vertex is in right trap region under isolation 3
		vert = vertStart;
		while (vert->NorthwestElem->Region->RegName != "Iso3")
		{
			vert = vert->NorthVertex;
		}
		upmostVert = vert;
		yDist = SctmMath::abs(vert->Y - vertStart->Y);

		vert = vertStart;
		while (vert->SouthwestElem->Region->RegName != "Trap.Gate2")
		{
			vert = vert->WestVertex;
		}
		westernmostVert = vert;
		xDist = SctmMath::abs(vert->X - vertStart->X);

		tanAngle = xDist / yDist;
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10052);
	}

}

void SlopingTunnelTrapToGate::setInterpolatedValues()
{
	static Mat::Name blockMat = domain->GetRegion("Block")->Mat->MatName();
	FDVertex *vert = vertStart;
	Normalization norm = Normalization(this->temperature);
	double deltaSlope = 0;
	double cbedge = 0;
	double mass = 0;
	
	double xDist = 0;
	double yDist = 0;
	double slopeDist = 0;
	double dy = 0;

	double dxLeft = 0;
	double dxRight = 0;
	double valueLeft = 0;
	double valueRight = 0;

	FDVertex *leftVert = NULL;
	FDVertex *rightVert = NULL;

	//set the first vertex
	cbedge = vertStart->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
	cbedge = norm.PullPotential(cbedge);
	this->cbEdge.push_back(cbedge);
	mass = vertStart->Phys->GetPhysPrpty(PhysProperty::eMass);
	this->elecMass.push_back(mass);

	vert = vertStart->NorthVertex;
	while (true) //note the iteration stop vertex
	{
		dy = SctmMath::abs(vert->Y - vert->SouthVertex->Y);
		deltaSlope = dy * SctmMath::sqrt(1 + tanAngle * tanAngle);
		this->dSlope.push_back(norm.PullLength(deltaSlope));

		yDist = SctmMath::abs(vert->Y - vertStart->Y);
		xDist = tanAngle * yDist;

		findLeftRigthVertex(vert, xDist, leftVert, rightVert);

		dxLeft = SctmMath::abs(xDist - SctmMath::abs(leftVert->X - vert->X));
		dxRight = SctmMath::abs(SctmMath::abs(rightVert->X - vert->X) - xDist);

		//set the cbedge
		valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
		valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
		cbedge = (valueLeft * dxRight + valueRight * dxLeft) / (dxLeft + dxRight);
		this->cbEdge.push_back(norm.PullPotential(cbedge));

		//set the emass
		valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		mass = (valueLeft * dxRight + valueRight * dxLeft) / (dxLeft + dxRight);
		this->elecMass.push_back(mass);

		if (vert == upmostVert)
		{
			break;
		}
		vert = vert->NorthVertex;
	}
}

void SlopingTunnelTrapToGate::findLeftRigthVertex(FDVertex *vert, double xDist, FDVertex* &left, FDVertex* &right)
{
	double distance = 0;
	FDVertex *iterVert = vert;


	if (vertRegName == "Trap.Iso2")
	{
		while (true)
		{
			distance = SctmMath::abs(iterVert->X - vert->X);
			if (iterVert->IsAtContact())
			{
				break;
			}
			if (distance > xDist)
			{
				break;
			}
			iterVert = iterVert->EastVertex;
		}
		//up to here, iterVert stores the right vertex
		right = iterVert;
		left = iterVert->WestVertex;
	}
	else if (vertRegName == "Trap.Iso3")
	{
		while (true)
		{
			distance = SctmMath::abs(iterVert->X - vert->X);
			if (iterVert->IsAtContact())
			{
				break;
			}
			if (distance > xDist)
			{
				break;
			}
			iterVert = iterVert->WestVertex;
		}
		//up to here, iterVert stores the right vertex
		left = iterVert;
		right = iterVert->EastVertex;
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10052);
	}


	
}

FDVertex* SlopingTunnelTrapToGate::GetGateVertex()
{
	FDVertex *iterVert = upmostVert;
	
	if (vertRegName == "Trap.Iso2")
	{
		while (true)
		{
			if (iterVert->IsAtContact())
			{
				break;
			}
			iterVert = iterVert->EastVertex;
		}
	}
	else if (vertRegName == "Trap.Iso3")
	{
		while (true)
		{
			if (iterVert->IsAtContact())
			{
				break;
			}
			iterVert = iterVert->WestVertex;
		}
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10052);
	}

	return iterVert;
}
