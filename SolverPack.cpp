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
	this->poissonSolver = new TwoDimPoissonSolver(domain);
	this->tunnelSolver = new SubsToGateEletronTunnel(domain);
	this->ddSolver = new DriftDiffusionSolver(domain);

	this->mapPotential.clear();
}

void SolverPack::callIteration()
{
	UtilsTimeStep.GenerateNext();

	poissonSolver->SolvePotential();
	fetchPoissonResult();

	tunnelSolver->ReadInput(siFermiAboveCBedge);
	tunnelSolver->SolveTunnel_Interface();
	fetchTunnelResult();
	
	ddSolver->ReadInputCurrentBC(mapCurrDensFromTunnelLayer);
	ddSolver->SolveDD();
	fetchDDResult();

}

void SolverPack::fetchPoissonResult()
{
	poissonSolver->UpdatePotential();
}

void SolverPack::Run()
{
	callIteration();
	UtilsDebug.WritePoisson(domain);
	UtilsDebug.WriteBandInfo(domain);
	UtilsDebug.WriteDensity(domain);
}

void SolverPack::fakeFermiEnergy()
{
	siFermiAboveCBedge.clear();
	Normalization norm = Normalization();
	FDVertex *currVert = NULL;
	int vertID = 0;
	double val = 0.1;
	val = norm.PushEnergy(val);
	for (size_t iVert = 0; iVert != domain->GetVertices().size(); ++iVert)
	{
		currVert = domain->GetVertices().at(iVert);
		if (currVert->IsAtContact() && currVert->Contact->ContactName == "Channel")
		{
			vertID = currVert->GetID();
			siFermiAboveCBedge[vertID] = val;
		}
	}
}

void SolverPack::fetchTunnelResult()
{
	this->mapCurrDensFromTunnelLayer.clear();
	tunnelSolver->ReturnResult(mapCurrDensFromTunnelLayer);
}

void SolverPack::fetchDDResult()
{
	ddSolver->UpdateElecDens();
}
