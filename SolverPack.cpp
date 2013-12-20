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
	mapCurrDensFromTunnelLayer.clear();
	mapSiFermiAboveCBedge.clear();
	mapCurrDensCoeff.clear();
	mapChannelPotential.clear();
}

void SolverPack::callIteration()
{
	while (!UtilsTimeStep.End())
	{
		UtilsTimeStep.GenerateNext();
		UtilsTimer.Set();

		//solve substrate
		subsSolver->SolveSurfacePot();
		fetchSubstrateResult();
		UtilsData.WriteSubstrateResult(subsSolver);

		//solver Poisson equation
		poissonSolver->ReadChannelPotential(mapChannelPotential);
		poissonSolver->SolvePotential();
		fetchPoissonResult();
		UtilsData.WritePotential(domain->GetVertices());
		UtilsData.WriteBandInfo(domain->GetVertices());
		UtilsData.WriteElecField(domain->GetVertices());

		//solve tunneling problem in tunneling oxide
		tunnelOxideSolver->ReadInput(mapSiFermiAboveCBedge);
		tunnelOxideSolver->SolveTunnel();
		fetchTunnelOxideResult();
		UtilsData.WriteTunnelCurrentFromSubs(domain, mapCurrDensFromTunnelLayer);

		//solve tunneling problem in blocking oxide
		blockOxideSolver->SolveTunnel();
		fetchBlockOxideResult();

		//solve trapping
		trappingSolver->SolveTrap();
		fetchTrappingResult();
		UtilsData.WriteTrapOccupation(domain->GetDDVerts());

		//solver drift-diffusion equation
		ddSolver->SolveDD(mapCurrDensFromTunnelLayer, mapCurrDensCoeff);
		fetchDDResult();
		UtilsData.WriteTunnelCoeff(domain, mapCurrDensFromTunnelLayer, mapCurrDensCoeff);
		UtilsData.WriteElecDens(domain->GetDDVerts());
		UtilsData.WriteElecCurrDens(domain->GetDDVerts());

		//write the final result
		UtilsData.WriteTotalElecDens(domain->GetDDVerts());
		UtilsData.WriteFlatBandVoltageShift(domain);

		UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());
	}
}

void SolverPack::fetchPoissonResult()
{
	poissonSolver->UpdatePotential();
}

void SolverPack::Run()
{
	callIteration();

	//UtilsDebug.WritePoisson(domain);
	//UtilsDebug.WriteBandInfo(domain);
	//UtilsDebug.WriteDensity(domain);
}

void SolverPack::fetchTunnelOxideResult()
{
	//it is critical to clear the map
	this->mapCurrDensFromTunnelLayer.clear();
	tunnelOxideSolver->ReturnResult(mapCurrDensFromTunnelLayer);
	//set the sign of boundary current for dd solver.
	for (VertexMapDouble::iterator it = mapCurrDensFromTunnelLayer.begin(); it != mapCurrDensFromTunnelLayer.end(); ++it)
	{
		//does not need to change
		//it->second = it->second;
	}
}

void SolverPack::fetchDDResult()
{
	ddSolver->UpdateElecDens();
}

void SolverPack::fetchBlockOxideResult()
{
	//it is critical to clear the map
	this->mapCurrDensCoeff.clear();
	blockOxideSolver->ReturnResult(mapCurrDensCoeff);
	//the current(or tunnCoeff) should be same with the direction of the boundary condition
	//so, reversed value is used to calculate current density
	FDVertex *currVert = NULL;
	int vertID = 0;
	for (VertexMapDouble::iterator it = mapCurrDensCoeff.begin(); it != mapCurrDensCoeff.end(); ++it)
	{
		it->second = - it->second;
		vertID = it->first;
		currVert = domain->GetVertex(vertID);
		//save the tunneling-out coefficient in this physics property, to be used in calculating the tunneling out current
		//because, in the boundary condition of the vertex, the tunneling coefficient is not stored.
		currVert->Phys->SetPhysPrpty(PhysProperty::TunnelCoeff, it->second);
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
