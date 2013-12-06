/**
* @file TrapSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-12-6   15:40
* @note
* @todo
*/
#include "TrapSolver.h"
#include "FDDomain.h"
#include "DomainDetails.h"
#include "SctmPhys.h"
#include "SctmUtils.h"

using SctmPhys::TrapProperty;

TrapSolver::TrapSolver(FDDomain *_domain): domain(_domain), vertices(_domain->GetDDVerts())
{

}

void TrapSolver::initializeSolver()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		eMobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eMobility)));
		eTrapDensMap.insert(VertexMapDouble::value_type(vertID, currVert->Trap->GetTrapPrpty(TrapProperty::eTrapDensity)));
		eXsectionMap.insert(VertexMapDouble::value_type(vertID, currVert->Trap->GetTrapPrpty(TrapProperty::eCrossSection)));

		//initialize the maps
		coeffMap.insert(VertexMapDouble::value_type(vertID, 0));
		rhsMap.insert(VertexMapDouble::value_type(vertID, 0));
		eTrappedMap.insert(VertexMapDouble::value_type(vertID, 0));
	}
}

void TrapSolver::setSolverTrapping()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double captureCoeff = 0;
	double eTrapDens = 0;
	double eDens = 0; // the electron density in conduction band
	
	double timeStep = 0;
	double coeff = 0;
	double rhs = 0;

	timeStep = SctmUtils::UtilsTimeStep.TimeStep();
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff);
		eTrapDens = eTrapDensMap[vertID];

		eDens = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);

		coeff = captureCoeff * timeStep * eDens;
		// += is used here
		coeffMap[vertID] += coeff;

		rhs = captureCoeff * timeStep * eDens * eTrapDens;
		rhsMap[vertID] += rhs;
	}
}

void TrapSolver::refreshSolver()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double eTrappedLastTime = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();
		eTrappedLastTime = currVert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
		//refresh the coefficient map and right-hand side map
		coeffMap[vertID] = 1;
		rhsMap[vertID] = eTrappedLastTime;
	}
}

void TrapSolver::SolveTrap()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double eTrapped = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();
		eTrapped = rhsMap[vertID] / coeffMap[vertID];
		eTrapDensMap[vertID] = eTrapped;
	}
}

void TrapSolver::UpdateTrappe()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double eTrapped = 0;
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		eTrapped = eTrappedMap[vertID];
		currVert->Trap->SetTrapPrpty(TrapProperty::eTrapped, eTrapped);
	}
}


