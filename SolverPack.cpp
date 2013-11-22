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
using namespace SctmUtils;

SolverPack::SolverPack(FDDomain *_domain): domain(_domain)
{
	initialize();
}

void SolverPack::initialize()
{
	this->poissonSolver = new TwoDimPoissonSolver(domain);
	//this->tunnelSolver = new SubsToGateEletronTunnel(domain);
	this->ddSolver = NULL;

	this->retPotential.clear();
}

void SolverPack::callIteration()
{
	poissonSolver->SolvePotential();
	updateWithPoissonResult();
}

void SolverPack::updateWithPoissonResult()
{
	poissonSolver->ReturnResult(this->retPotential);

	using SctmPhys::PhysProperty;
	FDVertex *vert = NULL;
	int vertID = 0;
	double val = 0;
	
	for (VertexMapDouble::iterator it = retPotential.begin(); it != retPotential.end(); ++it )
	{
		vertID = it->first;
		val = it->second;
		vert = domain->GetVertex(vertID);
		vert->Phys->UpdateValue(PhysProperty::ElectrostaticPotential, val);
	}
}

void SolverPack::Run()
{
	callIteration();
	UtilsDebug.WritePoisson(domain);
}
