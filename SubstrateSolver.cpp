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
#include "Normalization.h"

using SctmUtils::Normalization;

void OneDimSubsSolver::initializeSolver()
{
	using SctmUtils::SctmGlobalControl;

	this->temperature = SctmGlobalControl::Get().Temperature;
	Normalization norm = Normalization(temperature);
	
	//////////calculate the equilibrium electron and hole density
	subsDopConc = norm.PushDensity(SctmGlobalControl::Get().SubstrateDoping);
	if ( subsDopConc > 0)
	{
		//N-type
		subsType = NType;
		eDensEqui = subsDopConc;
		hDensEqui = 1 / eDensEqui; //in normalized value, i.e. in real value, it is ni*ni/subsDop
	}
	else
	{
		//P-type
		subsType = PType;
		subsDopConc = SctmMath::abs(subsDopConc);
		hDensEqui = SctmMath::abs(subsDopConc);
		eDensEqui = 1 / hDensEqui;
	}

	//////////calculate the gate capacitance
	FDContact *channelCont = NULL;
	channelCont = domain->GetContact("Channel");

	//TODO: this is an temporary method to get the required vertex
	FDVertex *startVert = channelCont->GetContactVerts().at(0);
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
	
	//double gateCap = norm.PullCapacitancePerArea(gateCapacitance);
}

void OneDimSubsSolver::calcFuncAndItsDeriv(double surfpot)
{
	using namespace MaterialDB;
	static double eps_Si = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_DielectricConstant);

	//double kT_div_q = SctmPhys::k0 * this->temperature / SctmPhys::q;
	double item_in_square_bracket = SctmMath::sqrt(hDensEqui * (SctmMath::exp(-surfpot) + surfpot - 1) + 
		eDensEqui * (SctmMath::exp(surfpot) - surfpot - 1));

	double numerator = hDensEqui * (-SctmMath::exp(-surfpot) + 1) + eDensEqui * (SctmMath::exp(surfpot) - 1);

	if (gateVoltage - flatbandVoltage > 0)
	{
		func_SurfPot = SctmMath::sqrt(2.0 * eps_Si) / gateCapacitance * item_in_square_bracket + surfpot - (gateVoltage - flatbandVoltage);
		funcDeriv_SurfPot = SctmMath::sqrt(eps_Si / 2.0) / gateCapacitance * numerator / item_in_square_bracket + 1;
	}
	else
	{
		func_SurfPot = -SctmMath::sqrt(2.0 * eps_Si) / gateCapacitance * item_in_square_bracket + surfpot - (gateVoltage - flatbandVoltage);
		funcDeriv_SurfPot = -SctmMath::sqrt(eps_Si / 2.0) / gateCapacitance * numerator / item_in_square_bracket + 1;
	}
	
}

double OneDimSubsSolver::solve_NewtonMethod()
{
	static double tolerance = 1e-7;
	static double eps = 1e-50;
	static int maxIterations = 1000;

	double guessSurfPot = 0; // in [V]
	if (gateVoltage - flatbandVoltage > 0)
	{
		guessSurfPot = 0.5;
	}
	else
	{
		guessSurfPot = -0.5;
	}

	Normalization norm = Normalization(temperature);

	double guessRoot = norm.PushPotential(guessSurfPot);
	double currRoot = guessRoot;
	double nextRoot = 0;
	int iteration = 0;
	double rootNotFound = true;

	while (iteration++ < maxIterations && rootNotFound)
	{
		calcFuncAndItsDeriv(currRoot);
		if (SctmMath::abs(funcDeriv_SurfPot) < eps)
		{
			// the denominator is too small
			SCTM_ASSERT(SCTM_ERROR, 10031);
		}
		nextRoot = currRoot - 0.5 * func_SurfPot / funcDeriv_SurfPot;
		if (SctmMath::abs((nextRoot - currRoot) / nextRoot) < tolerance)
		{
			rootNotFound = false;
		}
		else
		{
			currRoot = nextRoot;
		}
	}

	if (rootNotFound)
	{
		SCTM_ASSERT(SCTM_ERROR, 10032);
	}

	return nextRoot;
}

void OneDimSubsSolver::SolveSurfacePot()
{
	gateVoltage = domain->GetContact("Gate")->Voltage;
	calcFlatbandVoltage();
	this->surfacePotBend = solve_NewtonMethod();
	Normalization norm = Normalization(temperature);
	//double pot = norm.PullPotential(surfacePot);
	calcFermiAboveCB();
	calcChannelPotential();
}

OneDimSubsSolver::OneDimSubsSolver(FDDomain *_domain) : domain(_domain)
{
	initializeSolver();
}

void OneDimSubsSolver::calcFermiAboveCB()
{
	using namespace MaterialDB;
	double bandgap = GetMatPrpty(MaterialMap(Mat::Silicon), MatProperty::Mat_Bandgap);
	//double kT = SctmPhys::k0 * temperature;
	
	if (subsType == PType)
	{
		fermiAbove = -(bandgap / 2 + SctmMath::ln(subsDopConc)) + surfacePotBend;
	}
	else
	{
		fermiAbove = -(bandgap / 2 - SctmMath::ln(subsDopConc)) + surfacePotBend;
	}
	// todo with the normalization problem
}

void OneDimSubsSolver::calcChannelPotential()
{
	double subsPot = 0; // the potential at the substrate contact
	if (subsType = PType)
	{
		subsPot = SctmMath::asinh(-subsDopConc / 2.0);
	}
	else
	{
		subsPot = SctmMath::asinh(subsDopConc / 2);
	}
	channelPot = subsPot + surfacePotBend;
}

void OneDimSubsSolver::calcFlatbandVoltage()
{
	using namespace SctmUtils;
	Normalization norm = Normalization(temperature);
	double gateWorkFunction = norm.PushPotential(SctmGlobalControl::Get().GateWorkFunction);
	double workFuncDifference = gateWorkFunction - SctmPhys::ReferencePotential;

	double vfbShift_charge = SctmPhys::CalculateFlatbandShift_domain(domain);

	this->flatbandVoltage = vfbShift_charge + workFuncDifference;
}

void OneDimSubsSolver::ReturnResult(double &fermiAbove, double &channelPot)
{
	fermiAbove = this->fermiAbove;
	channelPot = this->channelPot;
}
