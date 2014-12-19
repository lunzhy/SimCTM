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

}

void OneDimSubsSolver::calcFuncAndItsDeriv(double surfpot)
{
	using namespace MaterialDB;
	static double eps_Si = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_DielectricConstant);

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
	//clear previous result
	fermiAboveMap.clear();
	channelPotMap.clear();
	
	//CAUTION! TODO: this should be revised according to the change of gate voltage during the simulation 
	gateVoltage = domain->GetContact("Gate")->Voltage;

	static FDContact *subsContact = domain->GetContact("Channel");
	static vector<FDVertex *> &channelVerts = subsContact->GetContactVerts();

	FDVertex *vert = NULL;
	bool toSolve = false;

	for (size_t iVert = 0; iVert != channelVerts.size(); ++iVert)
	{
		vert = channelVerts.at(iVert);

		toSolve = this->isSubsUnderGate(vert);
		if (!toSolve)
		{
			continue;
		}

		this->gateCapacitance = calcGateCapacitance(vert);
		this->flatbandVoltage = calcFlatbandVoltage(vert);
		this->surfacePotBend = solve_NewtonMethod();
		setFermiAboveCB(vert);
		setChannelPotential(vert);
	}
}

OneDimSubsSolver::OneDimSubsSolver(FDDomain *_domain) : domain(_domain)
{
	initializeSolver();
}

void OneDimSubsSolver::setFermiAboveCB(FDVertex *channelVert)
{
	using namespace MaterialDB;
	double bandgap = GetMatPrpty(GetMaterial(Mat::Silicon), MatProperty::Mat_Bandgap);
	//double kT = SctmPhys::k0 * temperature;
	double fermi_above = 0;

	if (subsType == PType)
	{
		fermi_above = -(bandgap / 2 + SctmMath::ln(subsDopConc)) + surfacePotBend;
	}
	else
	{
		fermi_above = -(bandgap / 2 - SctmMath::ln(subsDopConc)) + surfacePotBend;
	}
	// todo with the normalization problem
	fermiAboveMap[channelVert->GetID()] = fermi_above;
}

void OneDimSubsSolver::setChannelPotential(FDVertex *channelVert)
{
	double subsPot = 0; // the potential at the substrate contact
	double channel_pot = 0;
	if (subsType = PType)
	{
		subsPot = SctmMath::asinh(-subsDopConc / 2.0);
	}
	else
	{
		subsPot = SctmMath::asinh(subsDopConc / 2);
	}
	channel_pot = subsPot + this->surfacePotBend;
	channelPotMap[channelVert->GetID()] = channel_pot;
}

double OneDimSubsSolver::calcFlatbandVoltage(FDVertex *channelVert)
{
	using namespace SctmUtils;
	Normalization norm = Normalization(temperature);
	double gateWorkFunction = norm.PushPotential(SctmGlobalControl::Get().GateWorkFunction);
	double workFuncDifference = gateWorkFunction - SctmPhys::ReferencePotential;

	double vfbShift_charge = SctmPhys::CalculateFlatbandShift_slice_for1D(channelVert);

	return vfbShift_charge + workFuncDifference;
	//this->flatbandVoltage = vfbShift_charge + workFuncDifference;
}

void OneDimSubsSolver::ReturnResult(VertexMapDouble &_fermiAboveMap, VertexMapDouble &_channelPotMap)
{
	for (VertexMapDouble::iterator it = this->fermiAboveMap.begin(); it != this->fermiAboveMap.end(); ++it)
	{
		_fermiAboveMap[it->first] = it->second;
	}
	for (VertexMapDouble::iterator it = this->channelPotMap.begin(); it != this->channelPotMap.end(); ++it)
	{
		_channelPotMap[it->first] = it->second;
	}

}

double OneDimSubsSolver::calcGateCapacitance(FDVertex *channelVert)
{
	using SctmUtils::SctmGlobalControl;
	//////////calculate the gate capacitance in the slice
	FDVertex *currVert = channelVert;

	double epsilon = 0;
	double delta_d = 0;
	double cap_reciprocal = 0;

	while (currVert != NULL)
	{
		epsilon = currVert->Phys->GetPhysPrpty(PhysProperty::DielectricConstant);

		if (SctmGlobalControl::Get().Coordinate == "Cylindrical")
		{
			delta_d = currVert->R * SctmMath::ln((currVert->R + currVert->NorthLength / 2) /
				(currVert->R - currVert->SouthLength / 2));
		}
		else // SctmGlobalControl::Get().Coordinate == "Cartesian"
		{
			delta_d = (currVert->SouthLength + currVert->NorthLength) / 2;
		}

		cap_reciprocal += delta_d / epsilon;

		currVert = currVert->NorthVertex;
	}

	return 1.0 / cap_reciprocal;

	//this->gateCapacitance = 1 / cap_reciprocal;
	//double gateCap = norm.PullCapacitancePerArea(gateCapacitance);
}

bool OneDimSubsSolver::isSubsUnderGate(FDVertex* vert)
{
	//this is temporarily related to the specific structure of the cell.
	string contactName = "";
	size_t found = 0;
	while (vert != NULL)
	{
		if (vert->IsAtContact())
		{
			contactName = vert->Contact->ContactName;
			found = contactName.find("Gate");
			if (found != std::string::npos)
			{
				return true;
			}
		}
		vert = vert->NorthVertex;
	}
	return false;
}
