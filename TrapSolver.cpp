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
using SctmUtils::SctmGlobalControl;
using SctmUtils::SctmTimeStep;

TrapSolver::TrapSolver(FDDomain *_domain): domain(_domain), vertices(_domain->GetDDVerts())
{
	this->temperature = SctmGlobalControl::Get().Temperature;
	initializeSolver();
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
	double eFreeDens = 0; // the electron density in conduction band
	
	double timeStep = 0;
	double coeff = 0;
	double rhs = 0;

	static string captureModel = SctmGlobalControl::Get().TrapCaptureModel;
	//static string captureModel = "J-Model";

	timeStep = SctmTimeStep::Get().TimeStep();

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();


		if (captureModel == "J-Model")
		{
			captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_J_Model);
		}
		else if (captureModel == "V-Model")
		{
			captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_V_Model);
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10036);
		}
		
		eTrapDens = eTrapDensMap[vertID];

		eFreeDens = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);

		//coeff of nt due this trapping = capture coeff * delta_time * nf
		coeff = captureCoeff * timeStep * eFreeDens;
		coeffMap[vertID] += coeff; // += is used here

		//rhs of the equation = capture coeff * delta_time * nf * eTrapDens
		rhs = captureCoeff * timeStep * eFreeDens * eTrapDens;
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
	refreshSolver();
	setSolverTrapping();
	setSolverDetrapping_BasicSRH();
	solveEachVertex();
}

void TrapSolver::UpdateTrapped()
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

void TrapSolver::solveEachVertex()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double eTrapped = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();
		eTrapped = rhsMap[vertID] / coeffMap[vertID];
		eTrappedMap[vertID] = eTrapped;
	}
}

void TrapSolver::setSolverDetrapping_BasicSRH()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double timeStep = 0;	

	double coeff_detrapping = 0;
	double eEmission = 0;

	timeStep = SctmTimeStep::Get().TimeStep();
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		eEmission = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_BasicSRH);
		
		//coeff of nt due to detrapping = emission rate * dalta_time
		coeff_detrapping = timeStep * eEmission;
		coeffMap[vertID] += coeff_detrapping;

		//rhs the equation due to detrapping = 0; none effect
		rhsMap[vertID] += 0;
	}
}


