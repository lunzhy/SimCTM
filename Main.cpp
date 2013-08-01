#include "FDDomainTest.h"
int main()
{
	MaterialDB::SetMaterials();

	FDDomainTest aTest = FDDomainTest();
	aTest.BuildDomain();
	aTest.StuffPotential();
}