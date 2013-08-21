#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"
#include "PoissonSolver.h"

using namespace SctmUtils;
void initialize()
{
	MaterialDB::SetMaterials();
	SctmPhys::SetPhysConstant();
}

void SctmTest()
{
	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();
	SctmDebug::PrintDomainDetails(*aTest);

	TwoDimPoisson poisson = TwoDimPoisson(aTest);
	poisson.SolvePotential();
}
int main()
{
	initialize();
	SctmTest();
	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}