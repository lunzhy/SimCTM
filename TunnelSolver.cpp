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
#include "DDSolver.h"

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
	double Emin = bandEdgeTunnelFrom > bandEdgeTunnelTo ? bandEdgeTunnelFrom : bandEdgeTunnelTo; // in real value
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
	return DTFNdensity;
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
	return TEdensity;
}

TunnelSolver::TunnelSolver(FDDomain *_domain): domain(_domain)
{
	using SctmUtils::SctmGlobalControl;
	this->temperature = SctmGlobalControl::Get().Temperature;
	initialize();
}

void TunnelSolver::ReadInput(VertexMapDouble &fermi)
{
	int vertID = 0;
	double efa = 0; //electron fermi energy above
	double hfa = 0; //hole fermi energy above
	double bandgap = 0;
	FDVertex *vert = NULL;
	using MaterialDB::MatProperty;

	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		vert = vertsTunnelOxideStart.at(iVert);
		vertID = vert->GetID();

		SCTM_ASSERT(fermi.find(vertID) != fermi.end(), 10020);
		efa = fermi[vertID];
		
		this->efermiAboveMap[vertID] = efa;

		bandgap = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_Bandgap); //bandgap value is in normalized value
		hfa = (-efa) - bandgap;
		this->hfermiAboveMap[vertID] = hfa;
	}
}


void TunnelSolver::loadBandStructure(FDVertex *startVert)
{
	//clear and reset all the vectors
	bandEdge_Tunnel.clear();
	ehmass_Tunnel.clear();
	deltaX_Tunnel.clear();

	bandEdge_Trap.clear();
	ehMass_Trap.clear();
	deltaX_Trap.clear();
	trapEnergyLevel.clear();
	verts_Trap.clear();

	bandEdge_Block.clear();
	ehMass_Block.clear();
	deltaX_Block.clear();

	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// CAUTION! the tunneling direction is north
	Normalization norm = Normalization(this->temperature);
	double dx = 0;
	double ehmass = 0;
	double bandedge = 0;
	double trapLevel = 0;
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
		if (tunMode == TunnelMode::ElecTunnel)
		{
			if (currVert->Trap == NULL)
			{
				bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
			}
			else
			{
				bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, tunnelMat);
			}
		}
		else //tunModel == TunnelMode::HoleTunnel
		{
			if (currVert->Trap == NULL)
			{
				bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy);
			}
			else
			{
				bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, tunnelMat);
			}
		}
		
		bandedge = norm.PullEnergy(bandedge);

		//load emass
		if (tunMode == TunnelMode::ElecTunnel)
			ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		else //tunModel == TunnelMode::HoleTunnel
			ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::hMass);

		deltaX_Tunnel.push_back(dx);
		bandEdge_Tunnel.push_back(bandedge);
		ehmass_Tunnel.push_back(ehmass);

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
		if (tunMode == TunnelMode::ElecTunnel)
		{
			if (!currVert->Phys->HasMultiPrpty(PhysProperty::ConductionBandEnergy))
			{
				bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
			}
			else
			{
				bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, trapMat);
			}
		}
		else //tunModel == TunnelMode::HoleTunnel
		{
			if (!currVert->Phys->HasMultiPrpty(PhysProperty::ValenceBandEnergy))
			{
				bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy);
			}
			else
			{
				bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, trapMat);
			}
		}
		
		bandedge = norm.PullEnergy(bandedge);

		//load emass
		if (tunMode == TunnelMode::ElecTunnel)
			ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		else
			ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::hMass);

		//load trap energy level
		if (tunMode == TunnelMode::ElecTunnel)
			trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::eEnergyFromCondBand);
		else //tunModel == TunnelMode::HoleTunnel
			trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::hEnergyFromValeBand);
		
		trapLevel = bandedge - norm.PullEnergy(trapDepth);

		bandEdge_Trap.push_back(bandedge);
		deltaX_Trap.push_back(dx);
		ehMass_Trap.push_back(ehmass);
		trapEnergyLevel.push_back(trapLevel);
		verts_Trap.push_back(currVert);

		currVert = currVert->NorthVertex;
	}
	currVert = currVert->SouthVertex; // compensate the last statement in while loop

	//set the blocking oxide
	vertsBlockOxideStart.push_back(currVert);
	if (SctmGlobalControl::Get().LateralTunneling && SlopingTunnelTrapToGate::IsSlopingTunnel(currVert))
	{
		SlopingTunnelTrapToGate slopeTunnel = SlopingTunnelTrapToGate(this->domain, currVert, this->tunMode);
		slopeTunnel.LoadBandStructureAlongPath(deltaX_Block, bandEdge_Block, ehMass_Block);
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
			if (tunMode == TunnelMode::ElecTunnel)
			{
				if (currVert->Trap == NULL)
				{
					bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
				}
				else
				{
					bandedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
				}
			}
			else //tunModel == TunnelMode::HoleTunnel
			{
				if (currVert->Trap == NULL)
				{
					bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy);
				}
				else
				{
					bandedge = -currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, blockMat);
				}
			}
			
			bandedge = norm.PullEnergy(bandedge);

			//load emass
			if (tunMode == TunnelMode::ElecTunnel)
				ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
			else //tunModel == TunnelMode::HoleTunnel
				ehmass = currVert->Phys->GetPhysPrpty(PhysProperty::hMass);

			deltaX_Block.push_back(dx);
			bandEdge_Block.push_back(bandedge);
			ehMass_Block.push_back(ehmass);

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
	ret = prefactor * kB * T * SctmMath::exp(-q * (energy - bandEdgeTunnelFrom) / kB / T);
	return ret;
}



void SubsToTrapElecTunnel::setSolver_DTFN(FDVertex *startVertex)
{
	Normalization norm = Normalization(this->temperature);

	//set the silicon band edge, because the difference is fixed
	//all values of band edge are in real value, so pulling is needed.

	//set the fermi energy of the tunneling-in vertex
	//Pulling of the parameters is done here, because TunnelSolver uses real value internally
	double fermiAbove = norm.PullEnergy(efermiAboveMap[startVertex->GetID()]);

	if (tunDirection == TunnelDirection::North)
	{
		bandEdgeTunnelFrom = bandEdge_Tunnel.front() - this->eSubsBarrier;

		//set the conduction band edge in the trapping layer where the tunneling ends.
		bandEdgeTunnelTo = bandEdge_Trap.front();

		fermiEnergyTunnelFrom = bandEdgeTunnelFrom + fermiAbove;
	}
	else if (tunDirection == TunnelDirection::South)
	{
		bandEdgeTunnelFrom = bandEdge_Trap.front();
		bandEdgeTunnelTo = bandEdge_Tunnel.front() - this->eSubsBarrier;
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
	this->tunMode = TunnelSolver::ElecTunnel;

	using namespace MaterialDB;
	//substrate/tunnel oxide barrier
	Normalization norm = Normalization(this->temperature);
	this->eSubsBarrier = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElectronAffinity) 
						- GetMatPrpty(domain->GetRegion("Tunnel")->Mat, MatProperty::Mat_ElectronAffinity); //in normalized value
	this->eSubsBarrier = norm.PullPotential(eSubsBarrier);

	solverTAT = new SubsToTrapElecTAT(this);
}

double SubsToTrapElecTunnel::getSupplyFunction(double energy)
{
	if (tunDirection == TunnelDirection::North)
	{
		return TunnelSolver::supplyFunction_forCurrDens(energy);
	}
	else if (tunDirection == TunnelDirection::South)
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

	Normalization norm = Normalization(this->temperature);
	FDVertex* vert = NULL;
	double per_cm3_in_per_m3 = SctmPhys::per_cm3_in_per_m3;
	double edensity = 0;
	double currdens_TAT = 0;

	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		setTunnelDirection(vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
		setTunnelTag();

		if (tunDirection == TunnelDirection::North) // tunneling in
		{
			cbedge_max = bandEdge_Tunnel.front();
		}
		else if (tunDirection == TunnelDirection::South) // tunneling out
		{
			cbedge_max = bandEdge_Tunnel.back();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

		currdens = calcDTFNtunneling(deltaX_Tunnel, ehmass_Tunnel, bandEdge_Tunnel, cbedge_max);
		currdens += calcThermalEmission(deltaX_Tunnel, ehmass_Tunnel, bandEdge_Tunnel, cbedge_max);

		//for TAT current
		if (tunDirection == TunnelDirection::North) // tunnel from subs to trap layer
		{
			//calculate TAT current density, which is in [A/m^2]
			currdens_TAT = this->solverTAT->SolveTAT(vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
			currdens += currdens_TAT;
		}
		else // tunnel from trap layer to subs
		{
			//calculate TAT current density, which is in [A/m^2]
			currdens_TAT = this->solverTAT->SolveTAT(vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
			vert = this->vertsTunnelOxideEnd.at(iVert);
			edensity = vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
			edensity = norm.PullDensity(edensity) * per_cm3_in_per_m3;
			//current density is the tunneling coefficient of the vertex. currdens_TAT is the TAT current density,
			//so the equivalent tunneling coefficient is needed.
			currdens += currdens_TAT / edensity; // [A/m^2] / [/m^3] = [A.m]
		}
		
		if (SctmGlobalControl::Get().Coordinate == "Cylindrical")
		{
			currdens = updateTunCurrForCylindrical(currdens, vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
		}

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
	//This method is also used in SubsToTrapHoleTunnel, so both electron and hole tunneling are considered.

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
		
		if (tunMode == TunnelSolver::ElecTunnel)
			tunnelTag = currVert->BndCond.GetElecTunnelTag();
		else
			tunnelTag = currVert->BndCond.GetHoleTunnelTag();

		switch (tunnelTag)
		{
			case FDBoundary::eTunnelOut:
			{
				tunDirection = TunnelSolver::South;
				break;
			}
			case FDBoundary::eTunnelIn:
			{
				tunDirection = TunnelSolver::North;
				break;
			}
			case FDBoundary::hTunnelOut:
			{
				tunDirection = TunnelSolver::South;
				break;
			}
			case FDBoundary::hTunnelIn:
			{
				tunDirection = TunnelSolver::North;
				break;
			}
			default:
				SCTM_ASSERT(SCTM_ERROR, 10056);
				break;
		}

		if (tunDirection == TunnelDirection::North)
		{
			//currently the eCurrDens_Tunnel saves the current density, which is in A/m2
			currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
			ret[vertID] = norm.PushCurrDens(currDens);
		}
		else if (tunDirection == TunnelDirection::South)
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
	cbEdge_TunnelTrap.reserve(this->bandEdge_Tunnel.size() + this->bandEdge_Trap.size());
	cbEdge_TunnelTrap.insert(cbEdge_TunnelTrap.end(), this->bandEdge_Tunnel.begin(), this->bandEdge_Tunnel.end());
	cbEdge_TunnelTrap.insert(cbEdge_TunnelTrap.end(), this->bandEdge_Trap.begin(), this->bandEdge_Trap.end());

	deltaX_TunnelTrap.reserve(this->deltaX_Tunnel.size() + this->deltaX_Trap.size());
	deltaX_TunnelTrap.insert(deltaX_TunnelTrap.end(), this->deltaX_Tunnel.begin(), this->deltaX_Tunnel.end());
	deltaX_TunnelTrap.insert(deltaX_TunnelTrap.end(), this->deltaX_Trap.begin(), this->deltaX_Trap.end());

	eMass_TunnelTrap.reserve(this->ehmass_Tunnel.size() + this->ehMass_Trap.size());
	eMass_TunnelTrap.insert(eMass_TunnelTrap.end(), this->ehmass_Tunnel.begin(), this->ehmass_Tunnel.end());
	eMass_TunnelTrap.insert(eMass_TunnelTrap.end(), this->ehMass_Trap.begin(), this->ehMass_Trap.end());
}

void SubsToTrapElecTunnel::calcCurrDens_MFN()
{
	//TODO: write debug information for MFN
	if (bandEdgeTunnelFrom > bandEdgeTunnelTo)
	{
		//no Modified Fowler-Nordheim tunneling happens
		return;
	}

	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	double eCurrDens = 0; // in [A/m^2]

	//Emin and Emax are processed with q. (divided by q)
	//Emin and Emax are real values, in [eV]
	double Emin = bandEdgeTunnelFrom;
	double Emax = bandEdgeTunnelTo;

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
		if (bandEdge_Trap.at(iVert) <= energy)
		{
			//size = index + 1
			size = bandEdge_Tunnel.size() + iVert + 1;
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
		if ((trapEnergyLevel.at(iVert) >= energy) && (trapEnergyLevel.at(iVert + 1) < energy))
		{
			//size = vertices size of tunneling oxide + index + 1
			size = bandEdge_Tunnel.size() + iVert + 1;
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
	double Emin = bandEdgeTunnelFrom;
	double Emax = trapEnergyLevel.front();

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
	switch (tunDirection)
	{
		case TunnelDirection::North:
		{

			verts->BndCond.SetElecTunnelTag(FDBoundary::eTunnelIn);
			break;
		}
		case TunnelDirection::South:
		{
			verts->BndCond.SetElecTunnelTag(FDBoundary::eTunnelOut);
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
	SCTM_ASSERT(efermiAboveMap.find(vertIdSubs) != efermiAboveMap.end(), 10051);
	fermiSubs = bandEdge_Tunnel.front() - this->eSubsBarrier + norm.PullPotential(efermiAboveMap[vertIdSubs]);

	double edensity = 0;
	edensity = vertTrap->Phys->GetPhysPrpty(PhysProperty::eDensity);
	edensity = norm.PullDensity(edensity); //in [cm^-3] up to now
	edensity = edensity * per_cm3_in_per_m3; //in [m^3] up to now
	
	double eMassTrap = SctmPhys::m0 * GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_ElecDOSMass);
	double fermiRel = 0;//Ef below Ec
	//use bandgap when density equals to 0, indicating that fermi energy is far below conduction band
	fermiRel = (edensity != 0) ? k0 * this->temperature / q * SctmMath::ln(edensity * h * h * h / 8 / pi / (SctmMath::sqrt(pi) / 2) /
		SctmMath::pow(k0 * this->temperature, 1.5) / SctmMath::sqrt(2 * eMassTrap * eMassTrap * eMassTrap)) : 
		-norm.PullPotential(GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_Bandgap));
	fermiTrap = fermiRel + bandEdge_Trap.front();

	if (edensity == 0 || fermiSubs >= fermiTrap) //tunneling into the trap layer
	{
		tunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = fermiSubs;
		fermiEnergyTunnelTo = fermiTrap;

		bandEdgeTunnelFrom = bandEdge_Tunnel.front() - this->eSubsBarrier;
		bandEdgeTunnelTo = bandEdge_Trap.front();	
	}
	else
	{
		tunDirection = TunnelDirection::South;
		fermiEnergyTunnelFrom = fermiTrap; //however, this is useless in the calculation of tunneling-out coefficient
		fermiEnergyTunnelTo = fermiSubs;

		bandEdgeTunnelFrom = bandEdge_Trap.front();
		bandEdgeTunnelTo = bandEdge_Tunnel.front() - this->eSubsBarrier;
	}

	//set the effective DOS mass for tunneling
	//the material for substrate is silicon
	double eeffmass = 0;
	eeffmass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass);
	this->effTunnelMass = eeffmass; // the effective electron mass, in [m0]

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
	double cbedgeTrap_min = bandEdge_Trap.front();
	double cbedgeSubs = bandEdgeTunnelTo;
	double trapEnergy = 0;
	double tunnelCoeff = 0;
	int vertID = 0;
	int size = 0;
	FDVertex *currVert = NULL;
	double coeff_TAT2B = 0;

	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		trapEnergy = trapEnergyLevel.at(iVert);
		if (trapEnergy < cbedgeTrap_min && trapEnergy > cbedgeSubs)
		{
			//size = vertices size of tunneling oxide + index + 1
			size = bandEdge_Tunnel.size() + iVert + 1;
			currVert = verts_Trap.at(iVert);
			tunnelCoeff = getTransCoeff(trapEnergy, deltaX_TunnelTrap, eMass_TunnelTrap, cbEdge_TunnelTrap, size);

			vertID = currVert->GetID();
			eTransCoeffMap_T2B[vertID] = tunnelCoeff;

			coeff_TAT2B = this->solverTAT->CalcCoeff_TrapAssistedT2B(currVert);
			eCoeffMap_TAT2B[vertID] = coeff_TAT2B;

			eTransCoeffMap_T2B[vertID] += coeff_TAT2B;
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

double SubsToTrapElecTunnel::updateTunCurrForCylindrical(double currdens, FDVertex* inVert, FDVertex* outVert)
{
	double radiusTunFrom = 0;
	double radiusTunTo = 0;
	double updatedCurrDens = 0;

	if (this->tunDirection == TunnelSolver::North) // electron/hole tunneling in
	{
		radiusTunFrom = inVert->R;
		radiusTunTo = outVert->R;
		updatedCurrDens = currdens * radiusTunFrom / radiusTunTo;
	}
	return updatedCurrDens;
}


TrapToGateElecTunnel::TrapToGateElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	//the effective DOS mass for tunneling is set when determine the tunneling direction
	this->tunMode = TunnelSolver::ElecTunnel;
}

double TrapToGateElecTunnel::getSupplyFunction(double energy)
{
	if (tunDirection == TunnelDirection::North)
	{
		return TunnelSolver::supplyFunction_forTunCoeff(energy);
	}
	else if (tunDirection == TunnelDirection::South)
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
	bandEdgeTunnelFrom = bandEdge_Trap.back();

	//band energy of metal contact, equals to minus applied voltage.
	Normalization norm = Normalization(this->temperature);
	SCTM_ASSERT(endVertex->Contact != NULL, 10050);
	double gateVoltage = norm.PullPotential(endVertex->Contact->Voltage);
	bandEdgeTunnelTo = -gateVoltage;
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

		if (tunDirection == TunnelDirection::North) // tunneling out
		{
			cbedge_max = bandEdge_Block.front();
		}
		else if (tunDirection == TunnelDirection::South) // tunneling in
		{
			cbedge_max = bandEdge_Block.back();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

		currdens = calcDTFNtunneling(deltaX_Block, ehMass_Block, bandEdge_Block, cbedge_max);
		currdens += calcThermalEmission(deltaX_Block, ehMass_Block, bandEdge_Block, cbedge_max);

		if (SctmGlobalControl::Get().Coordinate == "Cylindrical")
		{
			currdens = updateTunCurrForCylindrical(currdens, vertsBlockOxideEnd.at(iVert), vertsBlockOxideStart.at(iVert));
		}

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

		if (tunMode == TunnelSolver::ElecTunnel)
			tunnelTag = currVert->BndCond.GetElecTunnelTag();
		else
			tunnelTag = currVert->BndCond.GetHoleTunnelTag();

		switch (tunnelTag)
		{
			case FDBoundary::eTunnelOut:
			{
				tunDirection = TunnelSolver::North;
				break;
			}
			case FDBoundary::eTunnelIn:
			{
				tunDirection = TunnelSolver::South;
				break;
			}
			case FDBoundary::hTunnelOut:
			{
				tunDirection = TunnelSolver::North;
				break;
			}
			case FDBoundary::hTunnelIn:
			{
				tunDirection = TunnelSolver::South;
				break;
			}
			case FDBoundary::noTunnel:
			{
				tunDirection = TunnelSolver::NoTunnel;
				break;
			}
			default:
				SCTM_ASSERT(SCTM_ERROR, 10056);
				break;
		}

		if (tunDirection == TunnelDirection::North)
		{
			//currently the eCurrDens_DTFN stores the coefficient to calculate current density
			//the calculated coefficient has a dimension of A*m, and should be converted to A*cm
			//this value[A*cm] * eDensity[cm^-3] = current density[A/cm^2]
			tunCoeff = eCurrDens_DTFN.at(iVert) / cm_in_m;
			ret[vertID] = norm.PushTunCoeff(tunCoeff);
		}
		else if (tunDirection == TunnelDirection::South)
		{
			currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
			ret[vertID] = norm.PushCurrDens(currDens);
		}
		else if (tunDirection == TunnelDirection::NoTunnel)
		{
			ret[vertID] = 0;
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
	switch (tunDirection)
	{
		case TunnelDirection::North:
		{
			verts->BndCond.SetElecTunnelTag(FDBoundary::eTunnelOut);
			break;
		}
		case TunnelDirection::South:
		{
			verts->BndCond.SetElecTunnelTag(FDBoundary::eTunnelIn);
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
	fermiTrap = fermiRel + bandEdge_Trap.back(); //the last vertex of trapping layer

	double effMass = 0;
	if (density == 0 || fermiGate > fermiTrap) //tunneling into the trap layer from gate
	{
		tunDirection = TunnelDirection::South;
		fermiEnergyTunnelFrom = fermiGate;
		fermiEnergyTunnelTo = fermiTrap;

		bandEdgeTunnelFrom = fermiGate - 10; //use a very large number (10eV) to indicate that gate has no edge
		bandEdgeTunnelTo = bandEdge_Trap.back();

		//set the effective mass for tunneling current
		this->effTunnelMass = 1; //the effective mass in gate is 1
	}
	else //tunneling out of the trap layer
	{
		tunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = fermiTrap;
		fermiEnergyTunnelTo = fermiGate - 10; //use a very large number (10eV) to indicate that gate has no edge

		bandEdgeTunnelFrom = bandEdge_Trap.back();
		bandEdgeTunnelTo = fermiGate; //for the gate, cbedgeTunnelTo is the same with the fermi energy of gate

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
	double cbedge_min = bandEdge_Trap.back();
	double gateEdge = bandEdgeTunnelTo;
	double trapEnergy = 0;
	double TC = 0;
	int vertID = 0;
	int startindex = 0;
	FDVertex *currVert = NULL;

	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		trapEnergy = this->trapEnergyLevel.at(iVert);
		if (trapEnergy < cbedge_min && trapEnergy > gateEdge)
		{
			startindex = iVert;
			currVert = verts_Trap.at(iVert);
			//the index is same for container of Trap vertices and TrapBlock vertices.
			TC = getTransCoeff(trapEnergy, deltaX_TrapBlock, eMass_TrapBlock, cbEdge_TrapBlock, 0, startindex);

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

	cbEdge_TrapBlock.reserve(bandEdge_Trap.size() + bandEdge_Block.size());
	cbEdge_TrapBlock.insert(cbEdge_TrapBlock.end(), bandEdge_Trap.begin(), bandEdge_Trap.end());
	cbEdge_TrapBlock.insert(cbEdge_TrapBlock.end(), bandEdge_Block.begin(), bandEdge_Block.end());

	eMass_TrapBlock.reserve(ehMass_Trap.size() + ehMass_Block.size());
	eMass_TrapBlock.insert(eMass_TrapBlock.end(), ehMass_Trap.begin(), ehMass_Trap.end());
	eMass_TrapBlock.insert(eMass_TrapBlock.end(), ehMass_Block.begin(), ehMass_Block.end());

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

double TrapToGateElecTunnel::updateTunCurrForCylindrical(double currdens, FDVertex* inVert, FDVertex* outVert)
{
	double radiusTunFrom = 0;
	double radiusTunTo = 0;
	double updatedCurrDens = 0;

	if (this->tunDirection == TunnelSolver::South) // electron/hole tunneling in
	{
		radiusTunFrom = inVert->R;
		radiusTunTo = outVert->R;
		updatedCurrDens = currdens * radiusTunFrom / radiusTunTo;
	}
	return updatedCurrDens;
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

SlopingTunnelTrapToGate::SlopingTunnelTrapToGate(FDDomain *_domain, FDVertex *_vert, TunnelSolver::TunnelMode _tunmode)
{
	domain = _domain;
	vertStart = _vert;
	tunMode = _tunmode;
	vertRegName = vertStart->SoutheastElem->Region->RegName;
	temperature = SctmGlobalControl::Get().Temperature;

	setBoundaryVerts();
	setInterpolatedValues();
}

void SlopingTunnelTrapToGate::LoadBandStructureAlongPath(vector<double> &dx, vector<double> &cbedge, vector<double> &emass)
{
	dx.push_back(this->dSlope.front() / 2);
	cbedge.push_back(this->bandEdge.front());
	emass.push_back(this->ehMass.front());

	for (size_t is = 0; is != this->dSlope.size() - 1; ++is)
	{
		dx.push_back((this->dSlope.at(is) + this->dSlope.at(is + 1)) / 2);
		cbedge.push_back(this->bandEdge.at(is + 1));
		emass.push_back(this->ehMass.at(is + 1));
	}

	dx.push_back(this->dSlope.back() / 2);
	cbedge.push_back(this->bandEdge.back());
	emass.push_back(this->ehMass.back());
	
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
	double bandedge = 0;
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
	if (tunMode == TunnelSolver::ElecTunnel)
		bandedge = vertStart->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
	else
		bandedge = -vertStart->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, blockMat);
	bandedge = norm.PullPotential(bandedge);
	this->bandEdge.push_back(bandedge);

	if (tunMode == TunnelSolver::ElecTunnel)
		mass = vertStart->Phys->GetPhysPrpty(PhysProperty::eMass);
	else
		mass = vertStart->Phys->GetPhysPrpty(PhysProperty::hMass);
	this->ehMass.push_back(mass);

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
		if (tunMode == TunnelSolver::ElecTunnel)
		{
			valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
			valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, blockMat);
		}
		else
		{
			valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, blockMat);
			valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, blockMat);
		}
		
		bandedge = (valueLeft * dxRight + valueRight * dxLeft) / (dxLeft + dxRight);
		this->bandEdge.push_back(norm.PullPotential(bandedge));

		//set the mass
		if (tunMode == TunnelSolver::ElecTunnel)
		{
			valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::eMass);
			valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		}
		else
		{
			valueLeft = leftVert->Phys->GetPhysPrpty(PhysProperty::hMass);
			valueRight = rightVert->Phys->GetPhysPrpty(PhysProperty::hMass);
		}
		mass = (valueLeft * dxRight + valueRight * dxLeft) / (dxLeft + dxRight);
		this->ehMass.push_back(mass);

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

void SubsToTrapHoleTunnel::setTunnelDirection(FDVertex* vertSubs, FDVertex* vertTrap)
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
	double hFermiSubs = 0;
	double hFermiTrap = 0; //the hole quasi-fermi energy of trapping layer

	int vertIdSubs = 0;
	vertIdSubs = vertSubs->GetID();
	
	//use hfermiAboveMap for holes, which is the offset value
	SCTM_ASSERT(hfermiAboveMap.find(vertIdSubs) != hfermiAboveMap.end(), 10051);
	hFermiSubs = bandEdge_Tunnel.front() - this->hSubsBarrier + norm.PullPotential(hfermiAboveMap[vertIdSubs]);

	double hdensity = 0;
	hdensity = vertTrap->Phys->GetPhysPrpty(PhysProperty::hDensity);
	hdensity = norm.PullDensity(hdensity); //in [cm^-3] up to now
	hdensity = hdensity * per_cm3_in_per_m3; //in [m^3] up to now

	double hMassTrap = SctmPhys::m0 * GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_HoleDOSMass);
	double hFermiRel = 0; //relative fermi level with respect to valence band edge
	//use bandgap when density equals to 0, indicating that fermi energy is far below valence band

	//Ev - Ef, Ev below Ef, the fermiRel is a negative value commonly
	hFermiRel = (hdensity != 0) ? k0 * this->temperature / q * SctmMath::ln(hdensity * h * h * h / 2 /
		SctmMath::pow((2 * pi * hMassTrap * k0 * this->temperature), 1.5)) :
		-norm.PullPotential(GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_Bandgap));
	hFermiTrap = hFermiRel + bandEdge_Trap.front(); //fermiRel is a negative value

	if (hdensity == 0 || hFermiSubs >= hFermiTrap) //tunneling into the trap layer
	{
		tunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = hFermiSubs;
		fermiEnergyTunnelTo = hFermiTrap;

		bandEdgeTunnelFrom = bandEdge_Tunnel.front() - this->hSubsBarrier;
		bandEdgeTunnelTo = bandEdge_Trap.front();
	}
	else
	{
		tunDirection = TunnelDirection::South;
		fermiEnergyTunnelFrom = hFermiTrap; //however, this is useless in the calculation of tunneling-out coefficient
		fermiEnergyTunnelTo = hFermiSubs;

		bandEdgeTunnelFrom = bandEdge_Trap.front();
		bandEdgeTunnelTo = bandEdge_Tunnel.front() - this->hSubsBarrier;
	}

	//set the effective DOS mass for tunneling
	//the material for substrate is silicon
	double heffmass = 0;
	heffmass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_HoleDOSMass);
	this->effTunnelMass = heffmass; // the effective electron mass, in [m0]
}

SubsToTrapHoleTunnel::SubsToTrapHoleTunnel(FDDomain* _domain) : SubsToTrapElecTunnel(_domain)
{
	this->tunMode = TunnelSolver::HoleTunnel;

	using namespace MaterialDB;
	Normalization norm = Normalization(this->temperature);
	
	//substrate/tunnel oxide valence barrier, positive value
	//material properties are in normalized value
	this->hSubsBarrier = (GetMatPrpty(domain->GetRegion("Tunnel")->Mat, MatProperty::Mat_ElectronAffinity) + 
		GetMatPrpty(domain->GetRegion("Tunnel")->Mat, MatProperty::Mat_Bandgap)) -
		(GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElectronAffinity) +
		GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_Bandgap));
	//convert to real value
	this->hSubsBarrier = norm.PullPotential(this->hSubsBarrier);
}


void SubsToTrapHoleTunnel::setTunnelTag()
{
	using namespace SctmUtils;
	FDVertex *verts = verts_Trap.front(); // the front element is vertex at tunnelOxide/trappingLayer interface
	switch (tunDirection)
	{
		case TunnelDirection::North:
		{
			verts->BndCond.SetHoleTunnelTag(FDBoundary::hTunnelIn);
			break;
		}
		case TunnelDirection::South:
		{
			verts->BndCond.SetHoleTunnelTag(FDBoundary::hTunnelOut);
			break;
		}
		case TunnelDirection::East:
		case TunnelDirection::West:
		default:
		SCTM_ASSERT(SCTM_ERROR, 10041);
		break;
	}
}

void SubsToTrapHoleTunnel::SolveTunnel()
{
	//the following containers are filled when loading band structure, so they are cleared at the beginning of
	//a new turn to solve tunnel problems.
	//Only vertsTunnelOxideStart exists.
	vertsTunnelOxideEnd.clear();
	vertsBlockOxideStart.clear();
	vertsBlockOxideEnd.clear();

	//clear the results of last calculation
	hCurrDens_DTFN.clear();
	hCurrDens_DTFN.resize(vertsTunnelOxideStart.size());

	double bandedge_max = 0;
	double currdens = 0;

	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		setTunnelDirection(vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
		setTunnelTag();

		if (tunDirection == TunnelDirection::North) // tunneling in
		{
			bandedge_max = bandEdge_Tunnel.front();
		}
		else if (tunDirection == TunnelDirection::South) // tunneling out
		{
			bandedge_max = bandEdge_Tunnel.back();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

		currdens = calcDTFNtunneling(deltaX_Tunnel, ehmass_Tunnel, bandEdge_Tunnel, bandedge_max);
		currdens += calcThermalEmission(deltaX_Tunnel, ehmass_Tunnel, bandEdge_Tunnel, bandedge_max);

		if (SctmGlobalControl::Get().Coordinate == "Cylindrical")
		{
			currdens = updateTunCurrForCylindrical(currdens, vertsTunnelOxideStart.at(iVert), vertsTunnelOxideEnd.at(iVert));
		}

		hCurrDens_DTFN.at(iVert) = currdens;

	}

	pretendToBeElecTun();
}

void SubsToTrapHoleTunnel::pretendToBeElecTun()
{
	//this method is used in order to use the methods in SubsToTrapElecTunnel for electron tunneling problem:
	//ReturnResult
	this->eCurrDens_DTFN = this->hCurrDens_DTFN;
}

TrapToGateHoleTunnel::TrapToGateHoleTunnel(FDDomain* _domain) : TrapToGateElecTunnel(_domain)
{
	this->tunMode = TunnelSolver::HoleTunnel;
}

void TrapToGateHoleTunnel::setTunnelDirection(FDVertex *vertTrap, FDVertex *vertGate)
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
	fermiGate = gateVoltage;

	double density = 0;
	density = vertTrap->Phys->GetPhysPrpty(PhysProperty::hDensity);
	density = norm.PullDensity(density); //in [cm^-3]
	density = density * per_cm3_in_per_m3; //in [m^-3]

	double hMassTrap = SctmPhys::m0 * GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_HoleDOSMass);
	double fermiRel = 0;
	//use bandgap when density equals to 0, indicating that fermi energy is far below conduction band
	fermiRel = (density != 0) ? k0 * this->temperature / q * SctmMath::ln(density * h * h * h / 8 / pi / (SctmMath::sqrt(pi) / 2) /
		SctmMath::pow(k0 * this->temperature, 1.5) / SctmMath::sqrt(2 * hMassTrap * hMassTrap * hMassTrap)) :
		-norm.PullPotential(GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_Bandgap));
	fermiTrap = fermiRel + bandEdge_Trap.back(); //the last vertex of trapping layer

	double effMass = 0;
	if (density == 0 || fermiGate > fermiTrap) //tunneling into the trap layer from gate
	{
		// there is no hole tunneling from gate to trap layer, the no calculation will be applied when South is met.
		tunDirection = TunnelDirection::NoTunnel;
	}
	else //tunneling out of the trap layer
	{
		tunDirection = TunnelDirection::North;
		fermiEnergyTunnelFrom = fermiTrap;
		fermiEnergyTunnelTo = fermiGate - 10; //use a very large number (10eV) to indicate that gate has no edge

		bandEdgeTunnelFrom = bandEdge_Trap.back();
		bandEdgeTunnelTo = fermiGate; //for the gate, cbedgeTunnelTo is the same with the fermi energy of gate

		//set the effective DOS mass for tunneling
		effMass = GetMatPrpty(GetMaterial(domain->GetTrapMatName()), MatProperty::Mat_HoleDOSMass);
		this->effTunnelMass = effMass;
	}
}

void TrapToGateHoleTunnel::setTunnelTag()
{
	using namespace SctmUtils;
	FDVertex *verts = verts_Trap.back(); //the back element is vertex at trappingLayer/blockingLayer interface
	switch (tunDirection)
	{
		case TunnelDirection::North:
		{
			verts->BndCond.SetHoleTunnelTag(FDBoundary::hTunnelOut);
			break;
		}
		case TunnelDirection::NoTunnel:
		{
			verts->BndCond.SetHoleTunnelTag(FDBoundary::noTunnel);
			break;
		}
		case TunnelDirection::South:
		case TunnelDirection::East:
		case TunnelDirection::West:
		default:
			SCTM_ASSERT(SCTM_ERROR, 10041);
			break;
	}
}

void TrapToGateHoleTunnel::SolveTunnel()
{
	//the following containers are filled when loading band structure, so they are cleared at the beginning of
	//a new turn to solve tunnel problems.
	//Only vertsTunnelOxideStart exists.
	vertsTunnelOxideEnd.clear();
	vertsBlockOxideStart.clear();
	vertsBlockOxideEnd.clear();

	//clear the results of last calculation
	hCurrDens_DTFN.clear();
	hCurrDens_DTFN.resize(vertsTunnelOxideStart.size());

	double currdens = 0;
	double vbedge_max = 0;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		if (!this->tunnelTrapToGateEnable)
		{
			continue;
		}

		setTunnelDirection(vertsBlockOxideStart.at(iVert), vertsBlockOxideEnd.at(iVert));
		setTunnelTag();

		if (tunDirection == TunnelDirection::North) // tunneling out
		{
			vbedge_max = bandEdge_Block.front();
			currdens = calcDTFNtunneling(deltaX_Block, ehMass_Block, bandEdge_Block, vbedge_max);
			currdens += calcThermalEmission(deltaX_Block, ehMass_Block, bandEdge_Block, vbedge_max);
			hCurrDens_DTFN.at(iVert) = currdens;
		}
		else if (tunDirection == TunnelDirection::NoTunnel) // tunneling in
		{
			//hole tunneling from gate to trap layer does not exist
			hCurrDens_DTFN.at(iVert) = 0;
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10041);
		}

	}
	pretendToBeElecTun();
}

void TrapToGateHoleTunnel::pretendToBeElecTun()
{
	eCurrDens_DTFN = hCurrDens_DTFN;
}

SubsToTrapElecTAT::SubsToTrapElecTAT(SubsToTrapElecTunnel* _tunSolver):
tunnelSolver(_tunSolver)
{
	setOxideTrapParam();
}

bool SubsToTrapElecTAT::calcTimeConstant(FDVertex* oxideVert, double& ctime, double& etime)
{
	using namespace MaterialDB;

	Normalization norm = Normalization(this->tunnelSolver->temperature);
	double cm_in_m = SctmPhys::cm_in_m;
	double cm2_in_m2 = cm_in_m * cm_in_m;
	double m0 = SctmPhys::m0;
	double k0 = SctmPhys::k0;
	double h = SctmPhys::h;
	double pi = SctmMath::PI;
	double q = SctmPhys::q;
	double kT = k0 * tunnelSolver->temperature;

	double crosssection = 0;
	double tunnelcoeff = 0;
	double activationEnergy = 0;
	double evth = 0;
	
	double trapEnergyLevel = 0;
	double bandedgeEnergyLevel = 0;
	double energy = 0;

	//double cbandEnergy = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
	//cbandEnergy = norm.PullEnergy(cbandEnergy);
	//if (this->tunnelSolver->bandEdgeTunnelFrom > cbandEnergy)
	//{
	//	return false;
	//}

	if (this->tunnelDirection == TunnelSolver::North) // tunnel from substrate to trap layer
	{
		//calculation of the capture time constant, use [m]
		tunnelcoeff = getTunCoeffTAT(this->tunStartVert, oxideVert, "in");
		double edensity = 2.0 *
			SctmMath::pow((2 * pi * m0 * GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * kT), 1.5)
			/ h / h / h *
			SctmMath::exp((tunnelSolver->fermiEnergyTunnelFrom - tunnelSolver->bandEdgeTunnelFrom) /
			(kT / q)); //density already in [1/m^3]

		double dosmass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * m0;
		evth = SctmMath::sqrt(3 * kT / dosmass); // in [m/s]

		crosssection = this->oxideTrapCrossSection * cm2_in_m2; // in [m^2]

		trapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) -
			norm.PushEnergy(this->oxideTrapEnergyFromCond); // in normalized value
		bandedgeEnergyLevel = norm.PushEnergy(this->tunnelSolver->bandEdgeTunnelFrom);

		energy = SctmMath::max_val(0, trapEnergyLevel - bandedgeEnergyLevel);

		/*energy = SctmMath::square(this->phononEnergy +
		(trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom)) / 4 / this->phononEnergy;
		energy = (trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom) * ((trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom)/4/this->phononEnergy + 0.5);*/

		ctime = 1 / (crosssection * tunnelcoeff * evth * edensity *
			SctmMath::exp(-energy));

		//ctime = SctmMath::exp(-this->activeEnergy) / (crosssection * tunnelcoeff * evth * density);

		//calculation of the emission time constant, use [cm]
		tunnelcoeff = getTunCoeffTAT(oxideVert, this->tunEndVert, "out");
		double dos = tunnelSolver->verts_Trap.front()->Phys->GetPhysPrpty(PhysProperty::eEffDOS); // in normalized value, in [1/cm^3]
		dos = norm.PullDensity(dos); // in [1/cm^3]

		evth = tunnelSolver->verts_Trap.front()->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity); // in normalized value, in [cm/s]
		evth = norm.PullVelocity(evth); // in [cm/s]

		crosssection = this->oxideTrapCrossSection; // in [cm^2]

		trapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) - 
			norm.PushEnergy(this->oxideTrapEnergyFromCond);
		bandedgeEnergyLevel = norm.PushEnergy(this->tunnelSolver->bandEdgeTunnelTo);

		energy = SctmMath::max_val(0, bandedgeEnergyLevel - trapEnergyLevel);

		/*energy = SctmMath::square(this->phononEnergy -
		(trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom)) / 4 / this->phononEnergy;
		energy = (trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom) * ((trapEnergyLevel - this->tunnelSolver->bandEdgeTunnelFrom) / 4 / this->phononEnergy - 0.5);*/

		etime = 1 / (crosssection * tunnelcoeff * evth * dos * SctmMath::exp(-energy));

		/*double trapEnergy = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) -
		this->trapEnergyFromCond;
		etime = SctmMath::exp(this->activeEnergy) / (crosssection * tunnelcoeff * evth * density) *
		SctmMath::exp(-(this->tunnelSolver->fermiEnergyTunnelTo));*/
	}
	else // tunnelDirection == TunnelSolver::South, tunnel from trap layer to substrate
	{
		//calculation of the capture time constant, use [cm]
		tunnelcoeff = getTunCoeffTAT(oxideVert, this->tunEndVert, "in");
		double edensity = tunnelSolver->verts_Trap.front()->Phys->GetPhysPrpty(PhysProperty::eDensity); // in normalized value, in [1/cm^3]
		edensity = norm.PullDensity(edensity); // in [1/cm^3]

		evth = tunnelSolver->verts_Trap.front()->Phys->GetPhysPrpty(PhysProperty::eThermalVelocity); // in normalized value, in [cm/s]
		evth = norm.PullVelocity(evth); // in [cm/s]
		crosssection = this->oxideTrapCrossSection; // in [cm^2]
		trapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) - 
			norm.PushEnergy(this->oxideTrapEnergyFromCond); // in normalized value, in [eV]
		bandedgeEnergyLevel = norm.PushEnergy(this->tunnelSolver->bandEdgeTunnelFrom);
		energy = SctmMath::max_val(0, trapEnergyLevel - bandedgeEnergyLevel);

		ctime = 1 / (crosssection * tunnelcoeff * evth * edensity * SctmMath::exp(-energy));

		//calculation of the emission time constant, use [m]
		tunnelcoeff = getTunCoeffTAT(this->tunStartVert, oxideVert, "out");
		double dos = 2.0 *
			SctmMath::pow((2 * pi * m0 * GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * kT), 1.5)
			/ h / h / h; //density already in [1/m^3]

		double dosmass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * m0;
		evth = SctmMath::sqrt(3 * kT / dosmass); // in [m/s]

		crosssection = this->oxideTrapCrossSection * cm2_in_m2; // in [m^2]
		trapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) - 
			norm.PushEnergy(this->oxideTrapEnergyFromCond); // in normalized value, in [eV]
		bandedgeEnergyLevel = norm.PushEnergy(this->tunnelSolver->bandEdgeTunnelTo);
		energy = SctmMath::max_val(0, bandedgeEnergyLevel - trapEnergyLevel);
		etime = 1 / (crosssection * tunnelcoeff * evth * dos * SctmMath::exp(-energy));
	}
	

	return true;
}

double SubsToTrapElecTAT::SolveTAT(FDVertex* tunStart, FDVertex* tunEnd)
{
	Normalization norm = Normalization(this->tunnelSolver->temperature);
	this->tunStartVert = tunStart;
	this->tunEndVert = tunEnd;
	loadOxideVertices();

	FDVertex* oxidevert = NULL;
	int vertID = 0;

	double ctime = 0;
	double etime = 0;
	double currDens = 0;
	double dx = 0; //in real value, in [m]
	double q = SctmPhys::q;
	double cm_in_m = SctmPhys::cm_in_m;
	double per_cm3_in_per_m3 = SctmPhys::per_cm3_in_per_m3;
	double oxideTrapDensity = 0; //in real value, in [m^-3]

	//set the tunneling direction before the calculation of time constant
	this->tunnelDirection = this->tunnelSolver->tunDirection;

	for (size_t iVert = 0; iVert != this->oxideVertices.size(); ++iVert)
	{
		oxidevert = oxideVertices.at(iVert);
		vertID = oxidevert->GetID();

		calcTimeConstant(oxidevert, ctime, etime);
		oxidevert->Phys->SetPhysPrpty(PhysProperty::eCaptureTime, norm.PushTime(ctime));
		oxidevert->Phys->SetPhysPrpty(PhysProperty::eEmissionTime, norm.PushTime(etime));

		dx = getDeltaX(oxidevert) * cm_in_m; // get the real value in [m]
		oxideTrapDensity = this->getOxideTrapDensity(oxidevert) * per_cm3_in_per_m3;; //in real value, in [cm^-3]
		currDens += q * oxideTrapDensity / (ctime + etime) * dx; //in real value, in [A/m^2]
	}

	return currDens;
}

void SubsToTrapElecTAT::loadOxideVertices()
{
	FDVertex* currVert = this->tunStartVert;

	this->oxideVertices.clear();
	while (true)
	{
		oxideVertices.push_back(currVert);
		if (currVert == this->tunEndVert)
		{
			break;
		}
		else
		{
			//previous knowledge of the structure has to be known.
			currVert = currVert->NorthVertex;
		}
	}
}

double SubsToTrapElecTAT::getTunCoeffTAT(FDVertex* startVert, FDVertex* endVert, string inout)
{
	int size = 0;
	int startindex = 0;
	double retTunCoeff = 0;
	double energy = 0;

	if (this->tunnelDirection == TunnelSolver::North) // tunnel from substrate to trap layer
	{
		if (inout == "in") //for tunneling from substrate into the trap
		{
			startindex = 0;
			size = std::find(oxideVertices.begin(), oxideVertices.end(), endVert) - oxideVertices.begin() + 1;
			if ((size_t)startindex > oxideVertices.size())
			{
				SCTM_ASSERT(SCTM_ERROR, 10059);
			}
			energy = this->tunnelSolver->bandEdgeTunnelFrom;
		}
		else // inout == "out", for tunneling from the trap into trapping layer
		{
			startindex = std::find(oxideVertices.begin(), oxideVertices.end(), startVert) - oxideVertices.begin();
			size = oxideVertices.size();
			energy = tunnelSolver->bandEdgeTunnelTo;
		}
		retTunCoeff = tunnelSolver->getTransCoeff(energy, tunnelSolver->deltaX_Tunnel, tunnelSolver->ehmass_Tunnel,
			tunnelSolver->bandEdge_Tunnel, size, startindex);
	}
	else // tunnel from trap layer to substrate
	{
		if (inout == "in")
		{
			startindex = std::find(oxideVertices.begin(), oxideVertices.end(), startVert) - oxideVertices.begin();
			size = oxideVertices.size();
			energy = tunnelSolver->bandEdgeTunnelFrom;
		}
		else // inout == "out"
		{
			startindex = 0;
			size = std::find(oxideVertices.begin(), oxideVertices.end(), endVert) - oxideVertices.begin() + 1;
			if ((size_t)startindex > oxideVertices.size())
			{
				SCTM_ASSERT(SCTM_ERROR, 10059);
			}
			energy = this->tunnelSolver->bandEdgeTunnelTo;
		}
		retTunCoeff = tunnelSolver->getTransCoeff(energy, tunnelSolver->deltaX_Tunnel, tunnelSolver->ehmass_Tunnel,
			tunnelSolver->bandEdge_Tunnel, size, startindex);
	}

	return retTunCoeff;
}

double SubsToTrapElecTAT::getDeltaX(FDVertex* vert)
{
	int pos = std::find(oxideVertices.begin(), oxideVertices.end(), vert) - oxideVertices.begin();
	double deltax =  tunnelSolver->deltaX_Tunnel.at(pos);
	// deltax is already in real value, but in cm
	return deltax;
}

void SubsToTrapElecTAT::setOxideTrapParam()
{
	using SctmUtils::ParamBase;
	using SctmUtils::SctmParameterParser;
	using SctmUtils::Param;

	Normalization norm = Normalization(this->tunnelSolver->temperature);
	ParamBase *parBase = NULL;

	//the calculation of TAT current is in real value,
	//so, all the parameters are in real values.

	//this->oxideTrapMaxDensity = 0; // in [1/cm^3]
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_trap_maxdens);
	this->oxideTrapMaxDensity = dynamic_cast<Param<double> *>(parBase)->Value();

	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_trap_pos);
	this->oxideTrapPosition = dynamic_cast<Param<double> *>(parBase)->Value();

	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_trap_sig);
	this->oxideTrapSigma = dynamic_cast<Param<double> *>(parBase)->Value();

	//this->oxideTrapCrossSection = 1e-13; // in  [1/cm^2]
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_trap_xsection);
	this->oxideTrapCrossSection = dynamic_cast<Param<double> *>(parBase)->Value();
	
	//this->oxideTrapEnergyFromCond = 1.6; // in real value, in [eV]
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_trap_energy);
	this->oxideTrapEnergyFromCond = dynamic_cast<Param<double> *>(parBase)->Value();
	
	//this->trapAssistedT2BTunnelFrequency = 1e6; // in real value, in [1/s]
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tat_t2b_frequency);
	this->trapAssistedT2BTunnelFrequency = dynamic_cast<Param<double> *>(parBase)->Value();
	
}

double SubsToTrapElecTAT::CalcCoeff_TrapAssistedT2B(FDVertex* trapvert)
{
	Normalization norm = Normalization(this->tunnelSolver->temperature);
	FDVertex* oxidevert = NULL;
	double q = SctmPhys::q;
	int vertID = 0;
	double ctime = 0;
	double etime = 0;
	double dx = 0; //in [cm]
	double dummy = 0;
	double oxideTrapDensity = 0; //in [cm^-3]
	double currdens = 0; //in [A/cm^2]
	double eTrappedDensity = 0;
	double frequency_T2B = 0;
	double tunCoeff_TAT2B = 0;

	for (size_t iVert = 0; iVert != this->oxideVertices.size(); ++iVert)
	{
		oxidevert = oxideVertices.at(iVert);
		vertID = oxidevert->GetID();

		calcTimeConstant_TrapAssistedT2B(oxidevert, trapvert, ctime, etime);
		dx = this->getDeltaX(oxidevert); //in real value, in [cm]
		oxideTrapDensity = this->getOxideTrapDensity(oxidevert); //in real value, in [cm^-3]

		currdens += q * oxideTrapDensity / (ctime + etime) * dx; //in real value, in [A/cm^2]
	}

	DriftDiffusionSolver::getDeltaXYAtVertex(trapvert, dummy, dx); //dx is the length in y direction
	dx = norm.PullLength(dx); //in [cm]
	eTrappedDensity = norm.PullDensity(trapvert->Trap->GetTrapPrpty(TrapProperty::eTrapped)); //in real value, in [cm^-3]

	frequency_T2B = trapvert->Trap->GetTrapPrpty(TrapProperty::eFrequency_T2B);
	frequency_T2B = norm.PullFrequency(frequency_T2B);

	//obtain the equivalent tunneling coefficient for trap-assisted trap-to-band tunneling current under the same
	//frequency of trap-to-band tunneling
	tunCoeff_TAT2B = currdens / q / dx / eTrappedDensity / frequency_T2B; //equivalent tunneling coefficient

	return tunCoeff_TAT2B;
}

void SubsToTrapElecTAT::calcTimeConstant_TrapAssistedT2B(FDVertex* oxideVert, FDVertex* trapVert, double& ctime, double& etime)
{
	using namespace MaterialDB;
	Normalization norm = Normalization(this->tunnelSolver->temperature);
	double cm_in_m = SctmPhys::cm_in_m;
	double cm2_in_m2 = cm_in_m * cm_in_m;
	double m0 = SctmPhys::m0;
	double k0 = SctmPhys::k0;
	double h = SctmPhys::h;
	double pi = SctmMath::PI;
	double q = SctmPhys::q;
	double kT = k0 * tunnelSolver->temperature;

	double crosssection = 0;
	double tunnelCoeff = 0;
	double frequency = 0;
	double dummy = 0;
	double dy = 0;
	double deltaX = 0;
	double eTrappedDensity = 0;
	double energy = 0;
	double oxideTrapEnergyLevel = 0;
	double trapSiteEnergyLevel = 0;
	double bandedgeEnergyLevel = 0;

	double dos = 0;
	double evth = 0;
	
	//capture time constant
	crosssection = this->oxideTrapCrossSection; //in [cm^2]
	tunnelCoeff = this->getTunCoeff_TrapAssistedT2B(oxideVert, trapVert);
	frequency = this->trapAssistedT2BTunnelFrequency;
	DriftDiffusionSolver::getDeltaXYAtVertex(trapVert, dummy, dy);
	deltaX = norm.PullLength(dy); //in [cm]
	eTrappedDensity = norm.PullDensity(trapVert->Trap->GetTrapPrpty(TrapProperty::eTrapped)); //in [cm^-3]
	oxideTrapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) - 
		norm.PushEnergy(this->oxideTrapEnergyFromCond); //in normalize value, in [eV]
	trapSiteEnergyLevel = trapVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) -
		trapVert->Trap->GetTrapPrpty(TrapProperty::eEnergyFromCondBand); //in normalized value, in [eV]
	energy = SctmMath::max_val(0, oxideTrapEnergyLevel - trapSiteEnergyLevel);
	ctime = 1 / (crosssection * tunnelCoeff * frequency * eTrappedDensity * deltaX * SctmMath::exp(-energy));
	
	//emission time constant
	tunnelCoeff = getTunCoeffTAT(this->tunStartVert, oxideVert, "out");
	dos = 2.0 *
		SctmMath::pow((2 * pi * m0 * GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * kT), 1.5)
		/ h / h / h; //density already in [1/m^3]

	double dosmass = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_ElecDOSMass) * m0;
	evth = SctmMath::sqrt(3 * kT / dosmass); // in [m/s]

	crosssection = this->oxideTrapCrossSection * cm2_in_m2; // in [m^2]
	oxideTrapEnergyLevel = oxideVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy) -
		norm.PushEnergy(this->oxideTrapEnergyFromCond); // in normalized value, in [eV]
	bandedgeEnergyLevel = norm.PushEnergy(this->tunnelSolver->bandEdgeTunnelTo);
	energy = SctmMath::max_val(0, bandedgeEnergyLevel - oxideTrapEnergyLevel);
	etime = 1 / (crosssection * tunnelCoeff * evth * dos * SctmMath::exp(-energy));
}

double SubsToTrapElecTAT::getTunCoeff_TrapAssistedT2B(FDVertex* oxideVert, FDVertex* trapVert)
{
	double energy = 0;
	int startindex = 0;
	int trapindex = 0;
	int size = 0;
	double tunCoeff = 0;

	startindex = std::find(this->oxideVertices.begin(), this->oxideVertices.end(), oxideVert) - this->oxideVertices.begin();
	trapindex = std::find(this->tunnelSolver->verts_Trap.begin(), this->tunnelSolver->verts_Trap.end(), trapVert) -
		this->tunnelSolver->verts_Trap.begin();
	size = this->tunnelSolver->bandEdge_Tunnel.size() + trapindex + 1;

	energy = this->tunnelSolver->trapEnergyLevel.at(trapindex); // in real value, in [eV]

	tunCoeff = this->tunnelSolver->getTransCoeff(energy, this->tunnelSolver->deltaX_TunnelTrap,
		this->tunnelSolver->eMass_TunnelTrap, this->tunnelSolver->cbEdge_TunnelTrap,
		size, startindex);
	return tunCoeff;
}

double SubsToTrapElecTAT::getOxideTrapDensity(FDVertex* vert)
{
	Normalization norm = Normalization(this->tunnelSolver->temperature);
	double trapDens = 0;
	double yCoord = 0;
	double cm_in_nm = SctmPhys::cm_in_nm;

	yCoord = norm.PullLength(vert->Y) * cm_in_nm;
	trapDens = this->oxideTrapMaxDensity *
		SctmMath::exp(-SctmMath::square(yCoord - this->oxideTrapPosition) / SctmMath::square(this->oxideTrapSigma) / 2);

	return trapDens;
}
