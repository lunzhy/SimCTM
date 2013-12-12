#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"
#include "DDSolver.h"
#include "SctmPhys.h"
#include "SolverPack.h"
#include "Normalization.h"
#include <stdlib.h>

using namespace SctmUtils;

void initialize()
{
	std::system("E:\\\"PhD Study\"\\SimCTM\\PySimFig\\DeleteData.py");
	UtilsMsg.PrintWelcomingInformation();
	UtilsMsg.PrintHeader("Initializing the simulator.");
	UtilsTimer.Start();
	
	//the initialization of the simulation goes here
	MaterialDB::SetMaterials();
	SctmPhys::SetPhysConstant();
}

void DomainTest()
{
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	UtilsDebug.PrintDomainDetails(aTest);
}

void PoissonTest()
{
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	UtilsDebug.PrintDomainDetails(aTest);

	UtilsMsg.PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aTest);
	poisson.SolvePotential();
	//UtilsDebug.PrintDomainDetails(*aTest);
}

void DDSolverTest()
{
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	//UtilsDebug.PrintDomainDetails(*aDomain);

	UtilsMsg.PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aDomain);
	//poisson.SolvePotential();

	UtilsMsg.PrintHeader("Testing the drift diffusion solver.");
	DDTest *ddSolver = new DDTest(aDomain);

	int i = 22;
	while ( i-->0)
	{
		UtilsTimeStep.GenerateNext();
		ddSolver->SolveDD();
	}
	//UtilsDebug.PrintDomainDetails(*aDomain);
	//UtilsTimeStep.GenerateNext();
	//ddSolver->SolveDD();
	//UtilsDebug.PrintDomainDetails(*aDomain);
}

void SolverPackTest()
{
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	SolverPack aPack = SolverPack(aDomain);
	aPack.Run();
}

void TimeStepTest()
{
	Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
	while (!UtilsTimeStep.End())
	{
		UtilsTimeStep.GenerateNext();
		UtilsDebug.PrintValue(UtilsTimeStep.StepNumber());
		UtilsDebug.PrintValue(norm.PullTime(UtilsTimeStep.TimeStep()));
		UtilsDebug.PrintValue(norm.PullTime(UtilsTimeStep.ElapsedTime()));
		UtilsDebug.PrintNewLine();
	}
}

int main()
{
	initialize();
	//DomainTest();
	SolverPackTest();
	//TimeStepTest();
	//DDSolverTest();
	//TunnelSolverTest();
	//PoissonTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}