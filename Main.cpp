#include "FDDomainTest.h"
#include "TunnelSolver.h"
int main()
{
	MaterialDB::SetMaterials();

	SimpleONO aTest = SimpleONO();
	aTest.BuildDomain();
	aTest.StuffPotential();

	SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	tunnelDemo.PrepareProblem(aTest.getVertex(0));
}