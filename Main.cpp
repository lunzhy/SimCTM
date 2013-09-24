#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"

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
	TunnelSolver *tunnelDemo = new TunnelTest();
	tunnelDemo->PrepareProblem(NULL);
	tunnelDemo->SolveTunneling();
}

int main()
{
	initialize();
	TunnelSolverTest();
	//SctmTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}