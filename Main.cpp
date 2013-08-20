#include "SimpleONO.h"
#include "TunnelSolver.h"
void initialize()
{
	MaterialDB::SetMaterials();
	SctmPhys::SetPhysConstant();
}
int main()
{
	initialize();

	SimpleONO aTest = SimpleONO();
	aTest.BuildDomain();

	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}