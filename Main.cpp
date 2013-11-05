#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"
#include "DDSolver.h"
#include "SctmPhys.h"

using namespace SctmUtils;
void initialize()
{
	UtilsMsg.PrintWelcomingInformation();
	UtilsMsg.PrintHeader("Initializing the simulator.");
	UtilsTimer.Start();
	
	//the initialization of the simulation goes here
	MaterialDB::SetMaterials();
	SctmPhys::SetPhysConstant();
}

void PoissonTest()
{
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	UtilsDebug.PrintDomainDetails(*aTest);

	UtilsMsg.PrintHeader("Solving potential using initial value.");
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aTest);
	poisson.SolvePotential();
	//UtilsDebug.PrintDomainDetails(*aTest);
}

void TunnelSolverTest()
{
	UtilsMsg.PrintHeader("Testing the tunneling solver.");
	TunnelTest *tunnelDemo = new TunnelTest();
	tunnelDemo->SolveParamterSet();
	tunnelDemo->SolveCalibrate();
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
	ddSolver->SolveDD();
	//UtilsDebug.PrintDomainDetails(*aDomain);
	ddSolver->SolveDD();
	//UtilsDebug.PrintDomainDetails(*aDomain);
	ddSolver->SolveDD();
}

int main()
{
	initialize();
	DDSolverTest();
	//TunnelSolverTest();
	//PoissonTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}