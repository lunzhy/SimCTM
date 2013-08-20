#include "SimpleONO.h"
#include "FDDomain.h"
#include "TunnelSolver.h"
#include "SctmUtils.h"

void initialize()
{
	MaterialDB::SetMaterials();
	SctmPhys::SetPhysConstant();
}
int main()
{
	initialize();

	FDDomain *aTest = new SimpleONO();
	aTest->BuildDomain();

	SctmUtils::SctmDebug::PrintDomainDetails(*aTest);

	//SubsToGateEletronTunnel tunnelDemo = SubsToGateEletronTunnel();
	//tunnelDemo.PrepareProblem(aTest.getVertex(0));
}