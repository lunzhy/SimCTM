#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"
#include "DDSolver.h"

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

void SctmTest()
{
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	UtilsDebug.PrintDomainDetails(*aTest);

	UtilsTimer.Set();
	TwoDimPoissonSolver poisson = TwoDimPoissonSolver(aTest);
	poisson.SolvePotential();

	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());
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
	UtilsTimer.Set();
	UtilsMsg.PrintHeader("Building a simple ONO domain.");
	FDDomain *aDomain = new SimpleONO();
	aDomain->BuildDomain();
	UtilsDebug.PrintDomainDetails(*aDomain);
	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());

	UtilsTimer.Set();
	UtilsMsg.PrintHeader("Testing the drift diffusion solver.");
	DDTest *ddSolver = new DDTest(aDomain);

	ddSolver->SolveDD();
	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());
}

int main()
{
	initialize();
	DDSolverTest();
	//TunnelSolverTest();
	//SctmTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}