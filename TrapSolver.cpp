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

ElecTrapSolver::ElecTrapSolver(FDDomain *_domain): domain(_domain), vertices(_domain->GetDDVerts())
{
	this->temperature = SctmGlobalControl::Get().Temperature;
	initializeSolver();
}

void ElecTrapSolver::initializeSolver()
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

void ElecTrapSolver::setSolverTrapping()
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

void ElecTrapSolver::refreshSolver()
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

void ElecTrapSolver::SolveTrap()
{
	refreshSolver();

	setSolverTrapping();
	setSolverDetrapping_BasicSRH();
	if (SctmGlobalControl::Get().PhysicsB2T)
	{
		setSolverBandToTrap();
	}
	if (SctmGlobalControl::Get().PhysicsT2B)
	{
		setSolverTrapToBand();
	}
	if (SctmGlobalControl::Get().PhysicsPFModel == "Frequency")
	{
		setSolverPooleFrenkel_Frequency();
	}
	solveEachVertex();
}

void ElecTrapSolver::UpdateTrapped()
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

void ElecTrapSolver::solveEachVertex()
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

void ElecTrapSolver::setSolverDetrapping_BasicSRH()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double timeStep = 0;	

	double coeff = 0;
	double rhs = 0;
	
	double eEmission = 0;

	timeStep = SctmTimeStep::Get().TimeStep();
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		//if PhysicsPFModel == EtDecrease, PF effect is included in the emission coefficient
		eEmission = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_BasicSRH);
		
		//coeff of nt due to detrapping = emission rate * dalta_time
		coeff = timeStep * eEmission;
		coeffMap[vertID] += coeff;

		//rhs of the equation due to detrapping = 0; none effect
		rhsMap[vertID] += 0;
	}
}

void ElecTrapSolver::setSolverBandToTrap()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double timeStep = SctmTimeStep::Get().TimeStep();
	double coeff_B2T = 0;

	double coeff = 0;
	double rhs = 0;

	double eTrapDens = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		//if the electron current density is negative, the electrons flow into the vertex
		//TODO: A temporary method is used here, because the sign is related to the current flow direction.
		coeff_B2T = - currVert->Trap->GetTrapPrpty(TrapProperty::eCoeff_B2T);
		
		//coeff of nt
		coeff = coeff_B2T * timeStep;
		coeffMap[vertID] += coeff;

		//rhs of the equation
		eTrapDens = eTrapDensMap[vertID];
		rhs = coeff_B2T * timeStep * eTrapDens;
		rhsMap[vertID] += rhs;
	}
}

void ElecTrapSolver::setSolverTrapToBand()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double coeff_T2B = 0;
	double coeff = 0;
	double rhs = 0;
	double timeStep = SctmTimeStep::Get().TimeStep();

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		//Trap-to-Band tunneling out always leads to the decrease of trapped charge
		coeff_T2B = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_T2B);

		//coeff of nt due to Trap-to-Band tunneling
		coeff = coeff_T2B * timeStep;
		coeffMap[vertID] += coeff;

		//rhs of the equation
		rhsMap[vertID] += 0;
	}
}

void ElecTrapSolver::setSolverPooleFrenkel_Frequency()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double timeStep = 0;

	double coeff = 0;
	double rhs = 0;

	double eEmission_pf = 0;

	timeStep = SctmTimeStep::Get().TimeStep();
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();
		eEmission_pf = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_PF);

		//coeff of nt due to detrapping = emission rate * dalta_time
		coeff = timeStep * eEmission_pf;
		coeffMap[vertID] += coeff;

		//rhs of the equation due to detrapping = 0; none effect
		rhs = 0;
		rhsMap[vertID] += rhs;
	}
}



HoleConserveTrapSolver::HoleConserveTrapSolver(FDDomain *_domain):
	domain(_domain), 
	vertices(_domain->GetDDVerts()), 
	temperature(SctmGlobalControl::Get().Temperature)
{
	FDVertex *currVert = NULL;
	int vertID = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		maphTrapDens.insert(VertexMapDouble::value_type(vertID, currVert->Trap->GetTrapPrpty(TrapProperty::hTrapDensity)));
		mapQuadEqu_A.insert(VertexMapDouble::value_type(vertID, 0));
		mapQuadEqu_B.insert(VertexMapDouble::value_type(vertID, 0));
		mapQuadEqu_C.insert(VertexMapDouble::value_type(vertID, 0));
	}
}

void HoleConserveTrapSolver::SolveTrap()
{
	this->timestep = SctmTimeStep::Get().TimeStep();

	setSolverTrap();
	solveEachVertex();

}

void HoleConserveTrapSolver::refreshSolver()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double val = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		val = currVert->Trap->GetTrapPrpty(TrapProperty::hTrapped);
		maphTrappedLasttime[vertID] = val;

		val = currVert->Phys->GetPhysPrpty(PhysProperty::hDensity);
		maphDensLasttime[vertID] = val;

		mapQuadEqu_A[vertID] = 0;
		mapQuadEqu_B[vertID] = -1;
		mapQuadEqu_C[vertID] = maphTrappedLasttime[vertID];
	}
}

void HoleConserveTrapSolver::setSolverTrap()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double val = 0;
	double captureCoeff = 0;
	double cap_time_product = 0;
	
	static string captureModel = SctmGlobalControl::Get().TrapCaptureModel;

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

		cap_time_product = this->timestep * captureCoeff;
		mapQuadEqu_A[vertID] += cap_time_product;
		mapQuadEqu_B[vertID] += -cap_time_product * (maphDensLasttime[vertID] + maphTrappedLasttime[vertID] + maphTrapDens[vertID]);
		mapQuadEqu_C[vertID] += cap_time_product * maphTrapDens[vertID] * (maphDensLasttime[vertID] + maphTrappedLasttime[vertID]);
	}
}

void HoleConserveTrapSolver::solveEachVertex()
{
	FDVertex *vert = NULL;
	int vertID = 0;
	double htrappedSolved;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		vert = vertices.at(iVert);
		vertID = vert->GetID();

		htrappedSolved = SctmMath::SolveQuadEquation(mapQuadEqu_A[vertID], mapQuadEqu_B[vertID], mapQuadEqu_C[vertID]);
		maphTrappedSolved[vertID] = htrappedSolved;
	}
}
