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
using namespace SctmUtils;

SolverPack::SolverPack(FDDomain *_domain): domain(_domain)
{
	initialize();
	fakeFermiEnergy();
}

void SolverPack::initialize()
{
	poissonSolver = new TwoDimPoissonSolver(domain);
	TunnelOxideSolver = new SubsToTrapElecTunnel(domain);
	BlockOxideSolver = new TrapToGateElecTunnel(domain);
	ddSolver = new DriftDiffusionSolver(domain);

	mapPotential.clear();
	mapCurrDensFromTunnelLayer.clear();
	mapSiFermiAboveCBedge.clear();
}

void SolverPack::callIteration()
{
	for (size_t it = 0; it != 200; ++it)
	{
		UtilsTimeStep.GenerateNext();

		poissonSolver->SolvePotential();
		fetchPoissonResult();
		UtilsData.WritePotential(domain->GetVertices());
		UtilsData.WriteBandInfo(domain->GetVertices());

		TunnelOxideSolver->ReadInput(mapSiFermiAboveCBedge);
		TunnelOxideSolver->SolveTunnel();
		fetchTunnelOxideResult();
		UtilsData.WriteTunnelCurrentFromSubs(domain, mapCurrDensFromTunnelLayer);
		UtilsData.WriteElecField(domain->GetVertices());
		UtilsData.WriteTotalElecDens(domain->GetDDVerts());

		BlockOxideSolver->SolveTunnel();
		fetchBlockOxideResult();

		ddSolver->ReadCurrDensBC_in(mapCurrDensFromTunnelLayer);
		ddSolver->ReadCurrDensBC_out(mapCurrDensCoeff);
		ddSolver->SolveDD();
		fetchDDResult();
		UtilsData.WriteTunnelCoeff(domain, mapCurrDensFromTunnelLayer, mapCurrDensCoeff);
		UtilsData.WriteElecDens(domain->GetDDVerts());
		UtilsData.WriteElecCurrDens(domain->GetDDVerts());
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
	TunnelOxideSolver->ReturnResult(mapCurrDensFromTunnelLayer);
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
	BlockOxideSolver->ReturnResult(mapCurrDensCoeff);
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
