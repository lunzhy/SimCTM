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
	simStructure = SctmGlobalControl::Get().Structure;
	initialize();
}

void SolverPack::initialize()
{
	subsSolver = new OneDimSubsSolver(domain);
	poissonSolver = new TwoDimPoissonSolver(domain);
	tunnelOxideElecSolver = new SubsToTrapElecTunnel(domain);
	tunnelOxideHoleSolver = new SubsToTrapHoleTunnel(domain);
	blockOxideElecSolver = new TrapToGateElecTunnel(domain);
	ddSolver = new DriftDiffusionSolver(domain);
	trappingSolver = new TrapSolver(domain);

	mapPotential.clear();
	mapElecCurrDensOrCoeff_Tunnel.clear();
	mapSiFermiAboveCBedge.clear();
	mapElecCurrDensOrCoeff_Block.clear();
	mapChannelPotential.clear();
}

void SolverPack::callIteration()
{
	if (SctmEnv::Get().IsLinux())
	{
		SctmPyCaller::PyBuildStructure();
	}
	SctmMessaging::Get().PrintHeader("Start to solve iterations.");
	SctmTimer::Get().Set();

	while (!SctmTimeStep::Get().End())
	{
		SctmTimeStep::Get().GenerateNext(); //the simulation starts with step 1
		SctmTimer::Get().Set();

		domain->RefreshGateVoltage();

		if (simStructure == "Single")
		{
			//solve substrate, no matter under Windows or Linux environment
			subsSolver->SolveSurfacePot();
			fetchSubstrateResult();
			SctmData::Get().WriteSubstrateResult(subsSolver);
		}
		else if (simStructure == "Triple" || simStructure == "TripleFull")
		{
			//write the flatband voltage shift of each slice
			SctmData::Get().WriteVfbShiftEachInterface(domain);
			//call Pytaurus to run Sentaurus to write substrate.in
			//Pytaurus will read the charge.in file
			if (SctmEnv::IsLinux())
			{
				SctmPyCaller::PySolve();
			}
			//if SimCTM is in running on Windows, read the same file, temporarily.
			readSubstrateFromFile();
			SctmData::Get().WriteSubstrateFromInput();
		}

		if (SctmTimeStep::Get().IsGateVoltageChanged() && SctmGlobalControl::Get().ClearCarrier)
		{
			domain->ClearCarrier();
		}

		//solver Poisson equation
		SctmTimer::Get().Set();
		poissonSolver->ReadChannelPotential(mapChannelPotential);
		poissonSolver->SolvePotential();
		SctmTimer::Get().Timeit("Poisson", SctmTimer::Get().PopLastSet());
		fetchPoissonResult();
		SctmData::Get().WritePotential(domain->GetVertices());
		SctmData::Get().WriteBandInfo(domain->GetVertices());
		SctmData::Get().WriteElecField(domain->GetVertices());
		
		//solve tunneling problem in tunneling oxide
		tunnelOxideElecSolver->ReadInput(mapSiFermiAboveCBedge);
		SctmTimer::Get().Set();
		tunnelOxideElecSolver->SolveTunnel();
		SctmTimer::Get().Timeit("Transport", SctmTimer::Get().PopLastSet());

		tunnelOxideHoleSolver->ReadInput(mapSiFermiAboveCBedge);
		SctmTimer::Get().Set();
		tunnelOxideHoleSolver->SolveTunnel();
		SctmTimer::Get().Timeit("Transport", SctmTimer::Get().PopLastSet());
		fetchTunnelOxideResult();
		SctmData::Get().WriteTunnelCurrentFromSubs(domain, mapElecCurrDensOrCoeff_Tunnel);

		//solve tunneling problem in blocking oxide
		SctmTimer::Get().Set();
		blockOxideElecSolver->SolveTunnel();
		SctmTimer::Get().Timeit("Transport", SctmTimer::Get().PopLastSet());
		fetchBlockOxideResult();

		//solver drift-diffusion equation
		SctmTimer::Get().Set();
		ddSolver->SolveDD(mapElecCurrDensOrCoeff_Tunnel, mapElecCurrDensOrCoeff_Block);
		SctmTimer::Get().Timeit("Transport", SctmTimer::Get().PopLastSet());
		fetchDDResult();
		SctmData::Get().WriteTunnelCoeff(domain, mapElecCurrDensOrCoeff_Tunnel, mapElecCurrDensOrCoeff_Block);
		SctmData::Get().WriteElecDens(domain->GetDDVerts());
		SctmData::Get().WriteElecCurrDens(domain->GetDDVerts());

		//solve trapping
		SctmTimer::Get().Set();
		trappingSolver->SolveTrap();
		SctmTimer::Get().Timeit("Transport", SctmTimer::Get().PopLastSet());
		fetchTrappingResult();
		SctmData::Get().WriteTrappedInfo(domain->GetDDVerts());
		
		//write the final result
		SctmData::Get().WriteTunnelOutDensity(domain, mapElecCurrDensOrCoeff_Tunnel, mapElecTransCoeffT2B_Tunnel, mapElecCurrDensOrCoeff_Block, mapElecTransCoeffT2B_Block);
		SctmData::Get().WriteTotalElecDens(domain->GetDDVerts());
		SctmData::Get().WriteTrappedDensRegionwise(domain);
		if (simStructure == "Single")
		{
			SctmData::Get().WriteFlatBandVoltageShift(domain);
		}

		SctmMessaging::Get().PrintTimeElapsed(SctmTimer::Get().PopLastSet());
	}

	SctmTimer::Get().Timeit("Total", SctmTimer::Get().PopLastSet());
	SctmData::Get().WriteTimerInfo(SctmTimer::Get());

	if (simStructure == "Triple" || simStructure == "TripleFull")
	{
		SctmPyCaller::PyParseAvgVfb();
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
	this->mapElecCurrDensOrCoeff_Tunnel.clear();
	tunnelOxideElecSolver->ReturnResult(mapElecCurrDensOrCoeff_Tunnel);
	//set the sign of boundary current for dd solver
	//the sign of the boundary condition is set according to the boundary direction
	for (VertexMapDouble::iterator it = mapElecCurrDensOrCoeff_Tunnel.begin(); it != mapElecCurrDensOrCoeff_Tunnel.end(); ++it)
	{	
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		if (vert->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelIn)
		{
			//does not need to change
			//it->second = it->second;	
		}
		else // eTunnelOut electron tunnel out of the trapping layer
		{
			it->second = -it->second;
			vert->Phys->SetPhysPrpty(PhysProperty::eTunnelCoeff, it->second);
		}
	}

	this->mapHoleCurrDensOrCoeff_Tunnel.clear();
	tunnelOxideHoleSolver->ReturnResult(mapHoleCurrDensOrCoeff_Tunnel);

	//

	//deal with MFN tunneling result
	//Assign the calculated current density to the specific vertex.
	double eCurrDens_MFN = 0; // in normalized value, in [A/cm^2]
	this->mapElecCurrDensMFN.clear();
	tunnelOxideElecSolver->ReturnResult_MFN(mapElecCurrDensMFN);
	for (VertexMapDouble::iterator it = mapElecCurrDensMFN.begin(); it != mapElecCurrDensMFN.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eCurrDens_MFN = -it->second;
		vert->Phys->SetPhysPrpty(PhysProperty::eCurrDensMFN_Y, eCurrDens_MFN);
	}

	//deal with Band-to-Trap tunneling result
	double eSubsCurrDens_B2T = 0; // in normalized value, in [A/cm^2]
	this->mapElecCurrDensB2T.clear();
	tunnelOxideElecSolver->ReturnResult_B2T(mapElecCurrDensB2T);
	for (VertexMapDouble::iterator it = mapElecCurrDensB2T.begin(); it != mapElecCurrDensB2T.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eSubsCurrDens_B2T = -it->second;
		vert->Phys->SetPhysPrpty(PhysProperty::eSubsCurrDensB2T, eSubsCurrDens_B2T);
	}

	//deal with Trap-to-Band tunneling out result
	double eTransCoeff_T2B = 0;
	this->mapElecTransCoeffT2B_Tunnel.clear();
	tunnelOxideElecSolver->ReturnResult_T2B(mapElecTransCoeffT2B_Tunnel);
	for (VertexMapDouble::iterator it = mapElecTransCoeffT2B_Tunnel.begin(); it != mapElecTransCoeffT2B_Tunnel.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		eTransCoeff_T2B = vert->Trap->GetTrapPrpty(TrapProperty::eTransCoeffT2B);
		eTransCoeff_T2B = eTransCoeff_T2B + it->second;//T2B occurs in two directions
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
	this->mapElecCurrDensOrCoeff_Block.clear();
	blockOxideElecSolver->ReturnResult(mapElecCurrDensOrCoeff_Block);
	//the current(or tunnCoeff) should be same with the direction of the boundary condition
	//so, reversed value is used to calculate current density
	FDVertex *vert = NULL;
	int vertID = 0;
	
	//set the sign of boundary current for dd solver
	//the sign of the boundary condition is set according to the boundary direction
	for (VertexMapDouble::iterator it = mapElecCurrDensOrCoeff_Block.begin(); it != mapElecCurrDensOrCoeff_Block.end(); ++it)
	{
		vertID = it->first;
		vert = domain->GetVertex(vertID);
		if (vert->BndCond.GetElecTunnelTag() == FDBoundary::eTunnelOut)
		{
			//set the tunneling coefficient in DT/FN tunneling out of electrons in conduction band
			//save the tunneling-out coefficient in the physics property, to be used in calculating the tunneling out current
			//because, in the boundary condition of the vertex, the tunneling coefficient is not stored.
			it->second = -it->second;
			vert->Phys->SetPhysPrpty(PhysProperty::eTunnelCoeff, it->second);
		}
		else
		{
			//does not need to change, because current direction is the same with the boundary direction
		}
	}

	//deal with trap-to-band tunneling out result
	double eTransCoeff_T2B = 0;
	this->mapElecTransCoeffT2B_Block.clear();
	blockOxideElecSolver->ReturnResult_T2B(mapElecTransCoeffT2B_Block);
	//set the trap-to-band tunneling out coefficient
	for (VertexMapDouble::iterator it = mapElecTransCoeffT2B_Block.begin(); it != mapElecTransCoeffT2B_Block.end(); ++it)
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
	subsSolver->ReturnResult(mapSiFermiAboveCBedge, mapChannelPotential);
}

void SolverPack::readSubstrateFromFile()
{
	SctmData::Get().ReadSubsInfoFromFile(mapSiFermiAboveCBedge, mapChannelPotential);
}
