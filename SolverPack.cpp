/**
* @file SolverPack.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-11-19   16:23
* @note
* @todo
*/
#include "FDDomain.h"
#include "SolverPack.h"
#include "PoissonSolver.h"
#include "TunnelSolver.h"
#include "DDSolver.h"
#include "SctmPhys.h"
#include "SctmUtils.h"
#include "DomainDetails.h"
#include "Normalization.h"
#include "TrapSolver.h"
#include "SubstrateSolver.h"

using namespace SctmUtils;

SolverPack::SolverPack(FDDomain *_domain): domain(_domain)
{
	this->temperature = SctmGlobalControl::Get().Temperature;
	initialize();
	//fakeFermiEnergy();
}

void SolverPack::initialize()
{
	subsSolver = new OneDimSubsSolver(domain);
	poissonSolver = new TwoDimPoissonSolver(domain);
	tunnelOxideSolver = new SubsToTrapElecTunnel(domain);
	blockOxideSolver = new TrapToGateElecTunnel(domain);
	ddSolver = new DriftDiffusionSolver(domain);
	trappingSolver = new TrapSolver(domain);

	mapPotential.clear();
	mapCurrDens_Tunnel.clear();
	mapSiFermiAboveCBedge.clear();
	mapCurrDensCoeff_Block.clear();
	mapChannelPotential.clear();
}

void SolverPack::callIteration()
{
	while (!SctmTimeStep::Get().End())
	{
		SctmTimeStep::Get().GenerateNext();
		SctmTimer::Get().Set();

		//solve substrate
		subsSolver->SolveSurfacePot();
		fetchSubstrateResult();
		SctmData::Get().WriteSubstrateResult(subsSolver);

		//solver Poisson equation
		poissonSolver->ReadChannelPotential(mapChannelPotential);
		poissonSolver->SolvePotential();
		fetchPoissonResult();
		SctmData::Get().WritePotential(domain->GetVertices());
		SctmData::Get().WriteBandInfo(domain->GetVertices());
		SctmData::Get().WriteElecField(domain->GetVertices());

		//solve tunneling problem in tunneling oxide
		tunnelOxideSolver->ReadInput(mapSiFermiAboveCBedge);
		tunnelOxideSolver->SolveTunnel();
		fetchTunnelOxideResult();
		SctmData::Get().WriteTunnelCurrentFromSubs(domain, mapCurrDens_Tunnel);

		//solve tunneling problem in blocking oxide
		blockOxideSolver->SolveTunnel();
		fetchBlockOxideResult();

		//solve trapping
		trappingSolver->SolveTrap();
		fetchTrappingResult();
		SctmData::Get().WriteTrapOccupation(domain->GetDDVerts());

		//solver drift-diffusion equation
		ddSolver->SolveDD(mapCurrDens_Tunnel, mapCurrDensCoeff_Block);
		fetchDDResult();
		SctmData::Get().WriteTunnelCoeff(domain, mapCurrDens_Tunnel, mapCurrDensCoeff_Block);
		SctmData::Get().WriteElecDens(domain->GetDDVerts());
		SctmData::Get().WriteElecCurrDens(domain->GetDDVerts());

		//write the final result
		SctmData::Get().WriteTotalElecDens(domain->GetDDVerts());
		SctmData::Get().WriteFlatBandVoltageShift(domain);

		SctmMessaging::Get().PrintTimeElapsed(SctmTimer::Get().SinceLastSet());
	}
}

void SolverPack::fetchPoissonResult()
{
	poissonSolver->UpdatePotential();
}

void SolverPack::Run()
{
	callIteration();

	//SctmDebug::GetInstance().WritePoisson(domain);
	//SctmDebug::GetInstance().WriteBandInfo(domain);
	//SctmDebug::GetInstance().WriteDensity(domain);
}

void SolverPack::fetchTunnelOxideResult()
{
	int vertID = 0;
	FDVertex *vert = NULL;

	//deal with Direct or FN tunneling result
	//it is critical to clear the map
	this->mapCurrDens_Tunnel.clear();
	tunnelOxideSolver->ReturnResult(mapCurrDens_Tunnel);
	//set the sign of boundary current for dd solver.
	for (VertexMapDouble::iterator it = mapCurrDens_Tunnel.begin(); it != mapCurrDens_Tunnel.end(); ++it)
	{	
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		if (vert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn)
		{
			//does not need to change
			//it->second = it->second;	
		}
		else // eTunnelOut electron tunnel out of the trapping layer
		{
			it->second = -it->second;
		}
	}


	//deal with MFN tunneling result
	//Assign the calculated current density to the specific vertex.
	double eCurrDens_MFN = 0; // in normalized value, in [A/cm^2]
	this->mapCurrDensMFN.clear();
	tunnelOxideSolver->ReturnResult_MFN(mapCurrDensMFN);
	for (VertexMapDouble::iterator it = mapCurrDensMFN.begin(); it != mapCurrDensMFN.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eCurrDens_MFN = -it->second;
		vert->Phys->SetPhysPrpty(PhysProperty::eCurrDensMFN_Y, eCurrDens_MFN);
	}

	//deal with Band-to-Trap tunneling result
	double eSubsCurrDens_B2T = 0; // in normalized value, in [A/cm^2]
	this->mapCurrDensB2T.clear();
	tunnelOxideSolver->ReturnResult_B2T(mapCurrDensB2T);
	for (VertexMapDouble::iterator it = mapCurrDensB2T.begin(); it != mapCurrDensB2T.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eSubsCurrDens_B2T = -it->second;
		vert->Phys->SetPhysPrpty(PhysProperty::eSubsCurrDensB2T, eSubsCurrDens_B2T);
	}

	//deal with Trap-to-Band tunneling out result
	double eTransCoeff_T2B = 0;
	this->mapTransCoeffT2B_Tunnel.clear();
	tunnelOxideSolver->ReturnResult_T2B(mapTransCoeffT2B_Tunnel);
	for (VertexMapDouble::iterator it = mapTransCoeffT2B_Tunnel.begin(); it != mapTransCoeffT2B_Tunnel.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eTransCoeff_T2B = vert->Trap->GetTrapPrpty(TrapProperty::eTransCoeffT2B);
		eTransCoeff_T2B = eTransCoeff_T2B + it->second;
		vert->Trap->SetTrapPrpty(TrapProperty::eTransCoeffT2B, eTransCoeff_T2B);
	}
}

void SolverPack::fetchDDResult()
{
	ddSolver->UpdateElecDens();
}

void SolverPack::fetchBlockOxideResult()
{
	//it is critical to clear the map
	this->mapCurrDensCoeff_Block.clear();
	blockOxideSolver->ReturnResult(mapCurrDensCoeff_Block);
	//the current(or tunnCoeff) should be same with the direction of the boundary condition
	//so, reversed value is used to calculate current density
	FDVertex *vert = NULL;
	int vertID = 0;
	
	//set the tunneling coefficient in DT/FN tunneling out of electrons in conduction band
	for (VertexMapDouble::iterator it = mapCurrDensCoeff_Block.begin(); it != mapCurrDensCoeff_Block.end(); ++it)
	{
		it->second = - it->second;
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		//save the tunneling-out coefficient in the physics property, to be used in calculating the tunneling out current
		//because, in the boundary condition of the vertex, the tunneling coefficient is not stored.
		vert->Phys->SetPhysPrpty(PhysProperty::TunnelCoeff, it->second);
	}

	double eTransCoeff_T2B = 0;
	this->mapTransCoeffT2B_Block.clear();
	blockOxideSolver->ReturnResult_T2B(mapTransCoeffT2B_Block);
	//set the trap-to-band tunneling out coefficient
	for (VertexMapDouble::iterator it = mapTransCoeffT2B_Block.begin(); it != mapTransCoeffT2B_Block.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eTransCoeff_T2B = vert->Trap->GetTrapPrpty(TrapProperty::eTransCoeffT2B);
		eTransCoeff_T2B = eTransCoeff_T2B + it->second;
		vert->Trap->SetTrapPrpty(TrapProperty::eTransCoeffT2B, eTransCoeff_T2B);
	}
}

void SolverPack::fetchTrappingResult()
{
	trappingSolver->UpdateTrapped();
}

void SolverPack::fetchSubstrateResult()
{
	double fermiAbove = 0;
	double channelPot = 0;
	subsSolver->ReturnResult(fermiAbove, channelPot);
	
	FDContact *channelContact = NULL;
	int vertID = 0;

	channelContact = domain->GetContact("Channel");
	std::vector<FDVertex *> &channelVerts = channelContact->GetContactVerts();
	
	for (size_t iVert = 0; iVert != channelVerts.size(); ++iVert)
	{
		vertID = channelVerts.at(iVert)->GetID();
		mapChannelPotential[vertID] = channelPot;
		mapSiFermiAboveCBedge[vertID] = fermiAbove;
	}
}
