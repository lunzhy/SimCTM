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

using namespace SctmUtils;

SolverPack::SolverPack(FDDomain *_domain): domain(_domain)
{
	initialize();
	fakeFermiEnergy();
}

void SolverPack::initialize()
{
	poissonSolver = new TwoDimPoissonSolver(domain);
	tunnelOxideSolver = new SubsToTrapElecTunnel(domain);
	blockOxideSolver = new TrapToGateElecTunnel(domain);
	ddSolver = new DriftDiffusionSolver(domain);
	trappingSolver = new TrapSolver(domain);

	mapPotential.clear();
	mapCurrDensFromTunnelLayer.clear();
	mapSiFermiAboveCBedge.clear();
}

void SolverPack::callIteration()
{
	while (!UtilsTimeStep.End())
	{
		UtilsTimeStep.GenerateNext();

		poissonSolver->SolvePotential();
		fetchPoissonResult();
		UtilsData.WritePotential(domain->GetVertices());
		UtilsData.WriteBandInfo(domain->GetVertices());
		UtilsData.WriteElecField(domain->GetVertices());

		tunnelOxideSolver->ReadInput(mapSiFermiAboveCBedge);
		tunnelOxideSolver->SolveTunnel();
		fetchTunnelOxideResult();
		UtilsData.WriteTunnelCurrentFromSubs(domain, mapCurrDensFromTunnelLayer);

		blockOxideSolver->SolveTunnel();
		fetchBlockOxideResult();

		trappingSolver->SolveTrap();
		fetchTrappingResult();
		UtilsData.WriteTrapOccupation(domain->GetDDVerts());

		ddSolver->SolveDD(mapCurrDensFromTunnelLayer, mapCurrDensCoeff);
		fetchDDResult();
		UtilsData.WriteTunnelCoeff(domain, mapCurrDensFromTunnelLayer, mapCurrDensCoeff);
		UtilsData.WriteElecDens(domain->GetDDVerts());
		UtilsData.WriteElecCurrDens(domain->GetDDVerts());

		UtilsData.WriteTotalElecDens(domain->GetDDVerts());
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

void SolverPack::fakeFermiEnergy()
{
	mapSiFermiAboveCBedge.clear();
	Normalization norm = Normalization();
	FDVertex *currVert = NULL;
	int vertID = 0;
	double val = 0.05;
	val = norm.PushEnergy(val);
	for (size_t iVert = 0; iVert != domain->GetVertices().size(); ++iVert)
	{
		currVert = domain->GetVertices().at(iVert);
		if (currVert->IsAtContact() && currVert->Contact->ContactName == "Channel")
		{
			vertID = currVert->GetID();
			mapSiFermiAboveCBedge[vertID] = val;
		}
	}
}

void SolverPack::fetchTunnelOxideResult()
{
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
