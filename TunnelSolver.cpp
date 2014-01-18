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

double TunnelSolver::getTransCoeff(double energy, vector<double> &deltax, vector<double> &emass, vector<double> &cbedge, int size /* = 0 */)
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

	for (vector<double>::size_type ix = 0; ix != size; ++ ix)
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

double TunnelSolver::calcDTFNtunneling(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge)
{
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	//Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	//Emin and Emax are real values, in [eV]
	//TODO: the choice of Emin determines the tunneling direction
	double Emin = cbedgeTunnelFrom > cbedgeTunnelTo ? cbedgeTunnelFrom : cbedgeTunnelTo; // in normalized value
	double Emax = cbedge.front();

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

double TunnelSolver::calcThermalEmission(vector<double> &deltaX, vector<double> &emass, vector<double> &cbedge)
{
	double TEdensity = 0; // in [A/m^2]
	double Emin = cbedge.front(); // in [eV]
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
	cbEdge_Oxide.clear();
	eMass_Oxide.clear();
	deltaX_Oxide.clear();

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
	static Mat::Name tunnelMat = domain->GetRegion(FDRegion::Tunneling)->Mat->MatName();
	static Mat::Name trapMat = domain->GetRegion(FDRegion::Trapping)->Mat->MatName();
	static Mat::Name blockMat = domain->GetRegion(FDRegion::Blocking)->Mat->MatName();

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

		deltaX_Oxide.push_back(dx);
		cbEdge_Oxide.push_back(cbedge);
		eMass_Oxide.push_back(emass);

		//check if current vertex is the ending vertex for tunneling layer
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsTunnelOxideEnd.push_back(currVert);
			//this tag is important in drift-diffusion solver.
			currVert->BndCond.SetTunnelTag(FDBoundary::eTunnelIn);
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
			vertsBlockOxideEnd.push_back(currVert);
			break;
		}
		currVert = currVert->NorthVertex;
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


void SubsToTrapElecTunnel::setSolver_DTFN(FDVertex *startVertex)
{
	Normalization norm = Normalization(this->temperature);

	//set the silicon band edge, because the difference is fixed
	//all values of band edge are in real value, so pulling is needed.
	using namespace MaterialDB;
	double barrier = 0;
	barrier = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(domain->GetRegion(FDRegion::Tunneling)->Mat, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelFrom = cbEdge_Oxide.front() - barrier;

	//set the conduction band edge in the trapping layer where the tunneling ends.
	cbedgeTunnelTo = cbEdge_Trap.front();
	
	//set the fermi energy of the tunneling-in vertex
	//Pulling of the parameters is done here, because TunnelSolver uses real value internally
	double fermiAbove = norm.PullEnergy(fermiAboveMap[startVertex->GetID()]);
	fermiEnergyTunnelFrom = cbedgeTunnelFrom + fermiAbove;
}

SubsToTrapElecTunnel::SubsToTrapElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	//set the effective tunneling mass
	using namespace MaterialDB;
	//the material for substrate is silicon
	double effSiMass = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_ElectronMass);
	this->effTunnelMass = effSiMass; // the effective electron mass, in [m0]
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
	//the following containers are filled when loading band structure, so they are cleared at the beginning of
	//a new turn to solve tunnel problems.
	//Only vertsTunnelOxideStart exists.
	vertsTunnelOxideEnd.clear();
	vertsBlockOxideStart.clear();
	vertsBlockOxideEnd.clear();

	//clear the results of last calculation
	eCurrDens_DTFN.resize(vertsTunnelOxideStart.size());
	eCurrDensMap_MFN.clear();
	eCurrDensMap_B2T.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		setTunnelTag();

		setSolver_DTFN(vertsTunnelOxideStart.at(iVert));

		currdens = calcDTFNtunneling(deltaX_Oxide, eMass_Oxide, cbEdge_Oxide);
		currdens += calcThermalEmission(deltaX_Oxide, eMass_Oxide, cbEdge_Oxide);
		eCurrDens_DTFN.at(iVert) = currdens;

		setSolver_Trap(vertsTunnelOxideStart.at(iVert));

		if (cbedgeTunnelFrom > cbedgeTunnelTo)
		{
			//no Modified Fowler-Nordheim tunneling happens
			continue;
		}
		//TODO: write debug information for MFN

		if (SctmGlobalControl::Get().PhysicsMFN)
		{
			//calcCurrDens_MFN();
		}
		
		//calcCurrDens_B2T();
	}
}

void SubsToTrapElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double currDens = 0;
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	Normalization norm = Normalization(this->temperature);
	for (size_t iVert = 0; iVert != vertsTunnelOxideEnd.size(); ++iVert)
	{
		//currVert = vertsEnd.at(iVert);
		currVert = vertsTunnelOxideEnd.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel saves the current density, which is in A/m2
		currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
		ret[vertID] = norm.PushCurrDens(currDens);
	}
}

void SubsToTrapElecTunnel::setSolver_Trap(FDVertex *startVertex)
{
	//reset the vectors
	cbEdge_Total.clear();
	eMass_Total.clear();
	deltaX_Total.clear();
	//
	cbEdge_Total.reserve(this->cbEdge_Oxide.size() + this->cbEdge_Trap.size());
	cbEdge_Total.insert(cbEdge_Total.end(), this->cbEdge_Oxide.begin(), this->cbEdge_Oxide.end());
	cbEdge_Total.insert(cbEdge_Total.end(), this->cbEdge_Trap.begin(), this->cbEdge_Trap.end());

	deltaX_Total.reserve(this->deltaX_Oxide.size() + this->deltaX_Trap.size());
	deltaX_Total.insert(deltaX_Total.end(), this->deltaX_Oxide.begin(), this->deltaX_Oxide.end());
	deltaX_Total.insert(deltaX_Total.end(), this->deltaX_Trap.begin(), this->deltaX_Trap.end());

	eMass_Total.reserve(this->eMass_Oxide.size() + this->eMass_Trap.size());
	eMass_Total.insert(eMass_Total.end(), this->eMass_Oxide.begin(), this->eMass_Oxide.end());
	eMass_Total.insert(eMass_Total.end(), this->eMass_Trap.begin(), this->eMass_Trap.end());
}

void SubsToTrapElecTunnel::calcCurrDens_MFN()
{
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
			TC = getTransCoeff(currEnergy, deltaX_Total, eMass_Total, cbEdge_Total, vSize);
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
	for (size_t iVert = 0; iVert != verts_Trap.size(); ++iVert)
	{
		if (cbEdge_Trap.at(iVert) <= energy)
		{
			//size = index + 1
			size = cbEdge_Oxide.size() + iVert + 1;
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
	for (size_t iVert = 0; iVert != verts_Trap.size() - 1; ++iVert)
	{
		if ((eEnergyLevel_Trap.at(iVert) <= energy) && (energy >= eEnergyLevel_Trap.at(iVert + 1)))
		{
			//size = index + 1
			size = cbEdge_Oxide.size() + iVert + 1;
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
			TC = getTransCoeff(currEnergy, deltaX_Total, eMass_Total, cbEdge_Total, vSize);
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
}


TrapToGateElecTunnel::TrapToGateElecTunnel(FDDomain *_domain): TunnelSolver(_domain)
{
	//set the effective tunneling mass
	using namespace MaterialDB;
	double effMass = GetMatPrpty(domain->GetRegion(FDRegion::Trapping)->Mat, MatProperty::Mat_ElectronMass);
	this->effTunnelMass = effMass; // the effective electron mass, in [m0]
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


void TrapToGateElecTunnel::setSolver_DTFN(FDVertex *endVertex)
{
	cbedgeTunnelFrom = cbEdge_Trap.back();


	//TODO: read this from simulation configuration
	//band energy of metal contact, equals to minus applied voltage.
	//double gateAppliedVoltage = -16; // in [V]
	using SctmUtils::SctmGlobalControl;
	//the value stored in SctmGlobalControl is in real value, no conversion is needed here.
	double gateAppliedVoltage = -SctmGlobalControl::Get().GateVoltage;
	cbedgeTunnelTo = gateAppliedVoltage;
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
	eCurrDens_DTFN.resize(vertsTunnelOxideStart.size());
	eCurrDensMap_T2B.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsTunnelOxideStart.size(); ++iVert)
	{
		loadBandStructure(vertsTunnelOxideStart.at(iVert));
		setTunnelTag();

		setSolver_DTFN(vertsTunnelOxideStart.at(iVert));

		currdens = calcDTFNtunneling(deltaX_Block, eMass_Block, cbEdge_Block);
		currdens += calcThermalEmission(deltaX_Block, eMass_Block, cbEdge_Block);
		eCurrDens_DTFN.at(iVert) = currdens;
	}
}

void TrapToGateElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double tunCoeff = 0;
	double cm_in_m = SctmPhys::cm_in_m;

	Normalization norm = Normalization(this->temperature);
	for (size_t iVert = 0; iVert != vertsBlockOxideStart.size(); ++iVert)
	{
		//currVert = vertsStart.at(iVert);
		currVert = vertsBlockOxideStart.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel stores the coefficient to calculate current density
		//the calculated coefficient has a dimension of A*m, and should be converted to A*cm
		//this value[A*cm] * eDensity[cm^-3] = current density[A/cm^2]
		tunCoeff = eCurrDens_DTFN.at(iVert) / cm_in_m;
		ret[vertID] = norm.PushTunCoeff(tunCoeff);
	}
}

void TrapToGateElecTunnel::setTunnelTag()
{
	double elecField = 0;
	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != vertsBlockOxideStart.size(); ++iVert)
	{
		currVert = vertsBlockOxideStart.at(iVert);
		//TODO: this is a temporary method to set the tunnel tags for vertices at interface.
		elecField = currVert->Phys->GetPhysPrpty(PhysProperty::ElectricField_Y);
		if (elecField < 0)
		{
			currVert->BndCond.SetTunnelTag(FDBoundary::eTunnelOut);
		}
	}
}
