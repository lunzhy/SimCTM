/**
* @file SubstrateSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-12-18   21:10
* @note
* @todo
*/

#include "SubstrateSolver.h"
#include "FDDomain.h"
#include "DomainDetails.h"
#include "SctmPhys.h"
#include "Material.h"
#include "SctmUtils.h"
#include "SctmMath.h"

void SubstrateSolver::initializeSolver()
{
	using SctmUtils::SctmGlobalControl;
	this->temperature = SctmGlobalControl::Get().Temperature;
	//////////calculate the equilibrium electron and hole density
	double subsDop = SctmGlobalControl::Get().SubstrateDoping;
	double ni = SctmPhys::IntrinsicConcentration;
	if ( subsDop > 0)
	{
		//N-type
		eDensEqui = subsDop;
		hDensEqui = ni*ni / eDensEqui;
	}
	else
	{
		//P-type
		hDensEqui = SctmMath::abs(subsDop);
		eDensEqui = ni*ni / hDensEqui;
	}

	//////////calculate the gate capacitance
	FDContact *channelCont = NULL;
	std::string channelContName = "Channel";
	channelCont = domain->GetContact(channelContName);

	//TODO: this is an temporary method to get the required vertex
	FDVertex *startVert = channelCont->vertices.at(0);
	FDVertex *currVert = startVert;

	double epsilon = 0;
	double delta_d = 0;
	double cap_reciprocal = 0;

	while (currVert != NULL)
	{
		epsilon = currVert->Phys->GetPhysPrpty(PhysProperty::DielectricConstant);
		delta_d = (currVert->NorthLength + currVert->SouthLength) / 2;
		cap_reciprocal += delta_d / epsilon;

		currVert = currVert->NorthVertex;
	}
	this->gateCapacitance = 1 / cap_reciprocal;
}

void SubstrateSolver::calcFuncAndItsDeriv(double surfpot)
{
	using namespace MaterialDB;
	double eps_Si = GetMatPrpty(MaterialMap[Mat::Silicon], MatProperty::Mat_DielectricConstant);

	FDContact *gateCont = NULL;
	this->gateVoltage = domain->GetContact("Gate")->Voltage;
	//calculate flatband voltage
	this->flatbandVoltage = 0;

	//double kT_div_q = SctmPhys::k0 * this->temperature / SctmPhys::q;
	double item_in_square_bracket = 0;

	item_in_square_bracket = SctmMath::sqrt(hDensEqui * (SctmMath::exp(-surfpot) + surfpot - 1) +
		eDensEqui * (SctmMath::exp(surfpot) - surfpot - 1));

	func_SurfPot = SctmMath::sqrt(2.0 * eps_Si) / gateCapacitance * item_in_square_bracket + surfpot - (gateVoltage - flatbandVoltage);
	double numerator = hDensEqui * (-SctmMath::exp(-surfpot) + 1) + eDensEqui*(SctmMath::exp(surfpot) - 1);
	funcDeriv_SurfPot = SctmMath::sqrt(eps_Si / 2.0) / gateCapacitance * numerator / item_in_square_bracket + 1;
}
