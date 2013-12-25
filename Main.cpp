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
#include <stdlib.h>

#include <iostream>
using namespace std;

using namespace SctmUtils;

void initialize()
{
	std::system("E:\\\"PhD Study\"\\SimCTM\\PySimFig\\DeleteData.py");
	SctmMessaging::GetInstance().PrintWelcomingInformation();
	SctmMessaging::GetInstance().PrintHeader("Initializing the simulator.");
	SctmTimer::GetInstance().Start();
	
	//the initialization of the simulation goes here
	//MaterialDB::SetMaterials_Directly();
	MaterialDB::SetMaterial_FromParFile();
	SctmPhys::SetPhysConstant();
}

void DomainTest()
{
	SctmMessaging::GetInstance().PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	SctmDebug::GetInstance().PrintDomainDetails(aTest);
}

void PoissonTest()
{
	SctmMessaging::GetInstance().PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	SctmDebug::GetInstance().PrintDomainDetails(aTest);

	SctmMessaging::GetInstance().PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aTest);
	poisson.SolvePotential();
	//SctmDebug::GetInstance().PrintDomainDetails(*aTest);
}

void DDSolverTest()
{
	SctmMessaging::GetInstance().PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);

	SctmMessaging::GetInstance().PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aDomain);
	//poisson.SolvePotential();

	SctmMessaging::GetInstance().PrintHeader("Testing the drift diffusion solver.");
	DDTest *ddSolver = new DDTest(aDomain);

	int i = 22;
	while ( i-->0)
	{
		SctmTimeStep::GetInstance().GenerateNext();
		ddSolver->SolveDD();
	}
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);
	//SctmTimeStep::GetInstance().GenerateNext();
	//ddSolver->SolveDD();
	//SctmDebug::GetInstance().PrintDomainDetails(*aDomain);
}

void SolverPackTest()
{
	SctmMessaging::GetInstance().PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	SolverPack aPack = SolverPack(aDomain);
	aPack.Run();
}

void TimeStepTest()
{
	Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
	while (!SctmTimeStep::GetInstance().End())
	{
		SctmTimeStep::GetInstance().GenerateNext();
		SctmDebug::GetInstance().PrintValue(SctmTimeStep::GetInstance().StepNumber());
		SctmDebug::GetInstance().PrintValue(norm.PullTime(SctmTimeStep::GetInstance().TimeStep()));
		SctmDebug::GetInstance().PrintValue(norm.PullTime(SctmTimeStep::GetInstance().ElapsedTime()));
		SctmDebug::GetInstance().PrintNewLine();
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

int main(int argc, char* argv[])
{
	cout << argv[0] << endl;
	//initialize();
	//ParaFileTest();
	//SubsSolverTest();
	//DomainTest();
	//SolverPackTest();
	//TimeStepTest();
	//DDSolverTest();
	//TunnelSolverTest();
	//PoissonTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}