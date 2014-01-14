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
#include "DDSolver.h"
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

double TunnelSolver::calcDTFNtunneling()
{
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC

	//the unit of calculated current density is [A/m^2]
	//Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	//Emin and Emax are real values, in [eV]
	double Emin = cbedgeTunnelFrom > cbedgeTunnelTo ? cbedgeTunnelFrom : cbedgeTunnelTo; // in normalized value
	double Emax = cbEdge_Oxide.front();

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
		TC = getTransCoeff(currEnergy, this->deltaX_Oxide, this->eMass_Oxide, this->cbEdge_Oxide);
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

double TunnelSolver::calcThermalEmission()
{
	double TEdensity = 0; // in [A/m^2]
	double Emin = cbEdge_Oxide.front(); // in [eV]
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
	this->eCurrDens = 0;
	this->temperature = SctmGlobalControl::Get().Temperature;
}

void TunnelSolver::ReadInput(VertexMapDouble &fermi)
{
	int vertID = 0;
	double val = 0;
	FDVertex *vert = NULL;
	for (size_t iVert = 0; iVert != vertsStart.size(); ++iVert)
	{
		vert = vertsStart.at(iVert);
		vertID = vert->GetID();

		SCTM_ASSERT(fermi.find(vertID) != fermi.end(), 10020);
		val = fermi[vertID];
		fermiAboveMap[vertID] = val;
	}
}

double TunnelSolver::calcCurrDens_DTFN()
{
	double ret = 0;
	ret = calcDTFNtunneling();
	ret += calcThermalEmission();
	return ret;
}


void SubsToTrapElecTunnel::initialize()
{
	//the tunneling start vertex is assigned here, because it only depends on the device structure and
	//maintains unchanged during the following calculation.
	this->vertsStart.clear();
	//this->vertsEnd_Trap.clear();

	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->domain->GetVertices().size(); ++iVert)
	{
		//the sequence in the vertices vector is its corresponding vertex index
		currVert = this->domain->GetVertex(iVert);
		if ( (currVert->IsAtContact()) && (currVert->Contact->ContactName == "Channel") )
		{
			vertsStart.push_back(currVert);
			//vertsEnd_Trap.push_back(currVert);
		}
	}
	//set the effective tunneling mass
	using namespace MaterialDB;
	//the material for substrate is silicon
	double effSiMass = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_ElectronMass);
	this->effTunnelMass = effSiMass; // the effective electron mass, in [m0]
}

void SubsToTrapElecTunnel::setSolver_DTFN(FDVertex *startVertex)
{
	//reset the vectors used in the calculation
	cbEdge_Oxide.clear();
	eMass_Oxide.clear();
	deltaX_Oxide.clear();
	//
	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// CAUTION! the tunneling direction is north
	Normalization norm = Normalization(this->temperature);
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
		//for boundary vertex of trapping layer, it always at the boundary of eDensity
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			vertsEnd.push_back(currVert);
			//this tag is important in drift-diffusion solver.
			currVert->BndCond.SetTunnelTag(FDBoundary::eTunnelIn);
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

		this->deltaX_Oxide.push_back(dx);
		this->eMass_Oxide.push_back(emass);
		this->cbEdge_Oxide.push_back(cbedge);

		currVert = currVert->NorthVertex; // to run the iteration
	}
	//CAUTION this is a temporary method
	//set the silicon band edge, because the difference is fixed
	//all values of band edge are in real value, so pulling is needed.
	using namespace MaterialDB;
	double barrier = 0;
	barrier = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(domain->GetRegion(FDRegion::Tunneling)->Mat, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelFrom = cbEdge_Oxide.front() - barrier;

	//set the conduction band edge in the trapping layer where the tunneling ends.
	barrier = GetMatPrpty(domain->GetRegion(FDRegion::Trapping)->Mat, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(domain->GetRegion(FDRegion::Tunneling)->Mat, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelTo = cbEdge_Oxide.back() - barrier;

	//new method
	//currVert = currVert->NorthVertex;
	//cbedgeTunnelTo = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
	//cbedgeTunnelTo = norm.PullEnergy(cbedgeTunnelTo);
	
	//set the fermi energy of the tunneling-in vertex
	//Pulling of the parameters is done here, because TunnelSolver uses real value internally
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
	//clear the results of last calculation
	vertsEnd.clear();
	eCurrDens_DTFN.resize(vertsStart.size());
	eCurrDensMap_MFN.clear();
	eCurrDensMap_B2T.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsStart.size(); ++iVert)
	{
		setSolver_DTFN(vertsStart.at(iVert));
		currdens = calcCurrDens_DTFN();
		eCurrDens_DTFN.at(iVert) = currdens;

		if (cbedgeTunnelFrom > cbedgeTunnelTo)
		{
			//no Modified Fowler-Nordheim tunneling happens
			continue;
		}

		//TODO: write debug information for MFN

		setSolver_Trap(vertsStart.at(iVert));
		calcCurrDens_MFN();

		calcCurrDens_B2T();
	}
}

void SubsToTrapElecTunnel::ReturnResult(VertexMapDouble &ret)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double currDens = 0;
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	Normalization norm = Normalization(this->temperature);
	for (size_t iVert = 0; iVert != vertsEnd.size(); ++iVert)
	{
		currVert = vertsEnd.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel saves the current density, which is in A/m2
		currDens = eCurrDens_DTFN.at(iVert) * per_m2_in_per_cm2;
		ret[vertID] = norm.PushCurrDens(currDens);
	}
}

void SubsToTrapElecTunnel::setSolver_Trap(FDVertex *startVertex)
{
	//the tunneling direction is north

	//reset the vectors used in the calculation
	cbEdge_Trap.clear();
	eMass_Trap.clear();
	deltaX_Trap.clear();
	verts_Trap.clear();
	eEnergyLevel_Trap.clear();
	//

	//cbedgeTunnelTo has to be set before this method, i.e. setSolver_Tunnel is called.
	Normalization norm = Normalization(this->temperature);
	double dx = 0;
	double emass = 0;
	double cbedge = 0;
	FDVertex *currVert = startVertex;

	double deltaX = 0;
	double deltaY = 0;

	double trapEnergyLevel = 0;
	double trapDepth = 0;

	while (true)
	{
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			break;
		}
		currVert = currVert->NorthVertex;
	}
	//currVert points to the vertex at trapping layer boundary now

	//cbedge = cbedgeTunnelTo; //the cbedge of first vertex along the trapping layer
	//emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
	//dx = norm.PullLength(deltaY);
	
	//move to next vertex in the trapping layer
	//currVert = currVert->NorthVertex;

	
	emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
	cbedge = cbedgeTunnelTo;
	trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::EnergyFromCondBand);

	DriftDiffusionSolver::getDeltaXYAtVertex(currVert, deltaX, deltaY);
	
	//the tunneling direction is from south to north, in y direction
	dx = norm.PullLength(deltaY);
	cbedge = norm.PullEnergy(cbedge);
	trapEnergyLevel = cbedge - norm.PullEnergy(trapDepth);

	verts_Trap.push_back(currVert);
	eMass_Trap.push_back(emass);
	cbEdge_Trap.push_back(cbedge);
	deltaX_Trap.push_back(dx);
	eEnergyLevel_Trap.push_back(trapEnergyLevel);

	currVert = currVert->NorthVertex;
	

	while (true)
	{
		if (currVert->Trap == NULL)
		{
			break;
		}
		emass = currVert->Phys->GetPhysPrpty(PhysProperty::eMass);
		cbedge = currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		trapDepth = currVert->Trap->GetTrapPrpty(TrapProperty::EnergyFromCondBand);

		DriftDiffusionSolver::getDeltaXYAtVertex(currVert, deltaX, deltaY);

		//the tunneling direction is from south to north, in y direction
		dx = norm.PullLength(deltaY);
		cbedge = norm.PullEnergy(cbedge);
		trapEnergyLevel = cbedge - norm.PullEnergy(trapDepth);

		verts_Trap.push_back(currVert);
		eMass_Trap.push_back(emass);
		cbEdge_Trap.push_back(cbedge);
		deltaX_Trap.push_back(dx);
		eEnergyLevel_Trap.push_back(trapEnergyLevel);

		if ((cbedge < cbedgeTunnelFrom && trapEnergyLevel < cbedgeTunnelFrom))
		{
			//this is a temporary method to find the last vertex along the tunneling path
			//break when the iteration meets the end of trapping layer
			break;
		}

		currVert = currVert->NorthVertex;
	}

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
	vertsEnd.clear();

	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->domain->GetVertices().size(); ++iVert)
	{
		//the sequence in the vertices vector is its corresponding vertex index
		currVert = this->domain->GetVertex(iVert);
		if (currVert->IsAtContact() && currVert->Contact->ContactName == "Gate")
		{
			vertsEnd.push_back(currVert);
			//vertsEnd_Trap.push_back(currVert);
		}
	}
	//set the effective tunneling mass
	//TODO: this should be obtained from user input
	using namespace MaterialDB;
	double effMass = GetMatPrpty(domain->GetRegion(FDRegion::Trapping)->Mat, MatProperty::Mat_ElectronMass);
	this->effTunnelMass = effMass; // the effective electron mass, in [m0]
}

void TrapToGateElecTunnel::setSolver_DTFN(FDVertex *endVertex)
{
	//reset the vectors used in the calculation
	cbEdge_Oxide.clear();
	eMass_Oxide.clear();
	deltaX_Oxide.clear();
	//
	// IMPORTANT! the parameters are in normalization values. They are converted before the calculation. (Not here) !
	// the tunneling direction is north, but the vertex find direction is south
	Normalization norm = Normalization(this->temperature);
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
			vertsStart.push_back(currVert);
			currVert->BndCond.SetTunnelTag(FDBoundary::eTunnelOut);
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

		this->deltaX_Oxide.push_back(dx);
		this->eMass_Oxide.push_back(emass);
		this->cbEdge_Oxide.push_back(cbedge);

		currVert = currVert->SouthVertex; // to run the iteration
	}

	std::reverse(deltaX_Oxide.begin(), deltaX_Oxide.end());
	std::reverse(eMass_Oxide.begin(), eMass_Oxide.end());
	std::reverse(cbEdge_Oxide.begin(), cbEdge_Oxide.end());

	//CAUTION this is a temporary method
	using namespace MaterialDB;
	double barrier = 0;
	barrier = GetMatPrpty(domain->GetRegion(FDRegion::Trapping)->Mat, MatProperty::Mat_ElectronAffinity)
		- GetMatPrpty(domain->GetRegion(FDRegion::Blocking)->Mat, MatProperty::Mat_ElectronAffinity);
	barrier = norm.PullEnergy(barrier);
	cbedgeTunnelFrom = cbEdge_Oxide.front() - barrier;

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
	eCurrDens_DTFN.resize(vertsEnd.size());
	vertsStart.clear();

	double currdens = 0;
	for (size_t iVert = 0; iVert != vertsEnd.size(); ++iVert)
	{
		setSolver_DTFN(vertsEnd.at(iVert));
		currdens = calcCurrDens_DTFN();
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
	for (size_t iVert = 0; iVert != vertsEnd.size(); ++iVert)
	{
		currVert = vertsStart.at(iVert);
		vertID = currVert->GetID();
		//currently the eCurrDens_Tunnel stores the coefficient to calculate current density
		//the calculated coefficient has a dimension of A*m, and should be converted to A*cm
		//this value[A*cm] * eDensity[cm^-3] = current density[A/cm^2]
		tunCoeff = eCurrDens_DTFN.at(iVert) / cm_in_m;
		ret[vertID] = norm.PushTunCoeff(tunCoeff);
	}
}
