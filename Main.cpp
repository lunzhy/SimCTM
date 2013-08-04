#include "FDDomainTest.h"
#include "TunnelSolver.h"
int main()
{
	MaterialDB::SetMaterials();

	FDDomainTest aTest = FDDomainTest();
	aTest.BuildDomain();
	aTest.StuffPotential();

	SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	tunnelDemo.PrepareProblem(aTest.getVertex(0));
}