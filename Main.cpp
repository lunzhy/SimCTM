#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"
#include "DDSolver.h"
#include "SctmPhys.h"
#include "SolverPack.h"
#include "Normalization.h"
#include "SubstrateSolver.h"
#include "TripleCells.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>

using namespace SctmUtils;

void initialize(const char *prjdir ="", const char *defaulParFile = "")
{
	SctmMessaging::Get().PrintWelcomingInformation();

	string prj(prjdir);
	string defaultParam(defaulParFile);

	if (prj.empty())
	{
		SctmMessaging::Get().PrintHeader("Cleaning debug folder.");
		prj = SctmEnv::Get().DebugPrjPath;
		SctmPyCaller::PyClean(prj);
	}
	else
	{
		if (SctmEnv::IsWindows())
		{
			//clean the project in Windows
			SctmPyCaller::PyClean(prj);
		}
		if (SctmEnv::IsLinux())
		{
			SctmMessaging::Get().PrintHeader("Preparing project folder.");
			SctmPyCaller::PyPrepare(prj);//currently the preparation of project folders only exists in Linux
		}
	}

	if (defaultParam.empty())
	{
		defaultParam = SctmEnv::Get().DefaultParamPath;
	}

	SctmGlobalControl::SetGlobalControl(defaultParam, prj);

	SctmMessaging::Get().PrintHeader("Initializing the simulator.");
	SctmTimer::Get().Start();
	
	//the initialization of the simulation goes here
	//MaterialDB::SetMaterials_Directly();
	MaterialDB::SetMaterial_FromParFile();
	SctmPhys::SetPhysConstant();
}

void DomainTest()
{
	SctmMessaging::Get().PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	SctmDebug::Get().PrintDomainDetails(aTest);
}

void TripleCellsDomainTest()
{
	FDDomain *aTriple = new TripleCellsFull();
	aTriple->BuildDomain();
	SctmDebug::Get().PrintDomainDetails(aTriple);
}

void PoissonTest()
{
	SctmMessaging::Get().PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	SctmDebug::Get().PrintDomainDetails(aTest);

	SctmMessaging::Get().PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aTest);
	poisson.SolvePotential();
	//SctmDebug::GetInstance().PrintDomainDetails(*aTest);
}

void DDSolverTest()
{
	SctmMessaging::Get().PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);

	SctmMessaging::Get().PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aDomain);
	//poisson.SolvePotential();

	SctmMessaging::Get().PrintHeader("Testing the drift diffusion solver.");
	DDTest *ddSolver = new DDTest(aDomain);

	int i = 22;
	while ( i-->0)
	{
		SctmTimeStep::Get().GenerateNext();
		ddSolver->SolveDD();
	}
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);
	//SctmTimeStep::GetInstance().GenerateNext();
	//ddSolver->SolveDD();
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);
}

void RunSolverPack()
{
	FDDomain *aDomain = NULL;
	if (SctmGlobalControl::Get().Structure == "TripleFull")
	{
		SctmMessaging::Get().PrintHeader("Building full triple-cell domain.");
		aDomain = new TripleCellsFull();
	}
	if (SctmGlobalControl::Get().Structure == "Single")
	{
		SctmMessaging::Get().PrintHeader("Building single-cell domain.");
		aDomain = new SimpleONO();
	}
	aDomain->BuildDomain();
	SolverPack aPack = SolverPack(aDomain);
	aPack.Run();
}

void TimeStepTest()
{
	Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
	while (!SctmTimeStep::Get().End())
	{
		SctmTimeStep::Get().GenerateNext();
		SctmDebug::Get().PrintValue(SctmTimeStep::Get().StepNumber());
		SctmDebug::Get().PrintValue(norm.PullTime(SctmTimeStep::Get().TimeStep()));
		SctmDebug::Get().PrintValue(norm.PullTime(SctmTimeStep::Get().ElapsedTime()));
		SctmDebug::Get().PrintNewLine();
	}
}

void SubsSolverTest()
{
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	OneDimSubsSolver subsSolver = OneDimSubsSolver(aDomain);
	subsSolver.SolveSurfacePot();
}

void ParaFileTest()
{
	SctmParameterParser parser = SctmParameterParser();
	double tem = SctmGlobalControl::Get().Temperature;
}

void PardisoTest()
{
	SctmMath::SctmSparseMatrixSolver solver;
	solver.PardisoTest();
}

int main(int argc, char* argv[])
{
	switch (argc)
	{
	case 1:
		initialize();
		break;
	case 2:
		initialize(argv[1]);
		break;
	case 3:
		initialize(argv[1], argv[2]);
		break;
	default:
		SctmMessaging::Get().PrintHeader("Argument Error");
		exit(0);
		break;
	}
	RunSolverPack();
	//TripleCellsDomainTest();
	//ParaFileTest();
	//SubsSolverTest();
	//DomainTest();
	//TimeStepTest();
	//DDSolverTest();
	//TunnelSolverTest();
	//PoissonTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
	return 0;
}