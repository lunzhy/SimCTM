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

TrapSolver::TrapSolver(FDDomain *_domain, TrapType _traptype): 
	domain(_domain), 
	vertices(_domain->GetDDVerts()),
	trapType(_traptype)
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

		if (this->trapType == TrapType::eTrap)
		{
			mapTrapDensity.insert(VertexMapDouble::value_type(vertID, currVert->Trap->GetTrapPrpty(TrapProperty::eTrapDensity)));
		}
		else // TrapType::hTrap
		{
			mapTrapDensity.insert(VertexMapDouble::value_type(vertID, currVert->Trap->GetTrapPrpty(TrapProperty::hTrapDensity)));
		}

		//initialize the maps
		coeffMap.insert(VertexMapDouble::value_type(vertID, 0));
		rhsMap.insert(VertexMapDouble::value_type(vertID, 0));
		mapTrappedSolved.insert(VertexMapDouble::value_type(vertID, 0));
	}
}

void TrapSolver::setSolverTrapping()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double captureCoeff = 0;
	double trapDens = 0;
	double freeDens = 0; // the electron density in conduction band
	
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

		if (this->trapType == TrapType::eTrap)
		{
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
			freeDens = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);
		}
		else
		{
			if (captureModel == "J-Model")
			{
				captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::hCaptureCoeff_J_Model);
			}
			else if (captureModel == "V-Model")
			{
				captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::hCaptureCoeff_V_Model);
			}
			else
			{
				SCTM_ASSERT(SCTM_ERROR, 10036);
			}
			freeDens = currVert->Phys->GetPhysPrpty(PhysProperty::hDensity);
		}

		
		
		trapDens = mapTrapDensity[vertID];

		//coeff of nt due this trapping = capture coeff * delta_time * nf
		coeff = captureCoeff * timeStep * freeDens;
		coeffMap[vertID] += coeff; // += is used here

		//rhs of the equation = capture coeff * delta_time * nf * eTrapDens
		rhs = captureCoeff * timeStep * freeDens * trapDens;
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
	if (this->trapType == TrapType::eTrap)
		eSolveTrap();
	else // TrapType::hTrap
		hSolveTrap();
}

void TrapSolver::UpdateTrapped()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double trappedSolved = 0;
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		trappedSolved = mapTrappedSolved[vertID];
		if (this->trapType == TrapType::eTrap)
			currVert->Trap->SetTrapPrpty(TrapProperty::eTrapped, trappedSolved);
		else // TrapType::hTrap
			currVert->Trap->SetTrapPrpty(TrapProperty::hTrapped, trappedSolved);
		
	}
}

void TrapSolver::solveEachVertex()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double trappedSolved = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();
		trappedSolved = rhsMap[vertID] / coeffMap[vertID];
		mapTrappedSolved[vertID] = trappedSolved;
	}
}

void TrapSolver::setSolverDetrapping_BasicSRH()
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

void TrapSolver::setSolverBandToTrap()
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
		eTrapDens = mapTrapDensity[vertID];
		rhs = coeff_B2T * timeStep * eTrapDens;
		rhsMap[vertID] += rhs;
	}
}

void TrapSolver::setSolverTrapToBand()
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

void TrapSolver::setSolverPooleFrenkel_Frequency()
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

void TrapSolver::eSolveTrap()
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

void TrapSolver::hSolveTrap()
{
	refreshSolver();

	setSolverTrapping();
// 	setSolverDetrapping_BasicSRH();
// 	if (SctmGlobalControl::Get().PhysicsB2T)
// 	{
// 		setSolverBandToTrap();
// 	}
// 	if (SctmGlobalControl::Get().PhysicsT2B)
// 	{
// 		setSolverTrapToBand();
// 	}
// 	if (SctmGlobalControl::Get().PhysicsPFModel == "Frequency")
// 	{
// 		setSolverPooleFrenkel_Frequency();
// 	}
	solveEachVertex();
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

	refreshSolver();
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
			captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::hCaptureCoeff_J_Model);
		}
		else if (captureModel == "V-Model")
		{
			captureCoeff = currVert->Trap->GetTrapPrpty(TrapProperty::hCaptureCoeff_V_Model);
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
	double htrappedSolved = 0;

	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		vert = vertices.at(iVert);
		vertID = vert->GetID();

		htrappedSolved = SctmMath::SolveQuadEquation(mapQuadEqu_A[vertID], mapQuadEqu_B[vertID], mapQuadEqu_C[vertID]);
		maphTrappedSolved[vertID] = htrappedSolved;
	}
}

void HoleConserveTrapSolver::UpdateTrapped()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	double trappedDens = 0;
	for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		vertID = currVert->GetID();

		trappedDens = maphTrappedSolved[vertID];
		currVert->Trap->SetTrapPrpty(TrapProperty::hTrapped, trappedDens);
	}
}
