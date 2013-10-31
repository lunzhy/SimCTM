/**
* @file TunnelSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-19   10:31
* @note
* @todo
*/

#include "TunnelSolver.h"
#include "SctmMath.h"
#include "SctmUtils.h"
#include "SctmPhys.h"

using SctmPhys::PhysProperty;
using SctmUtils::SctmFileStream;

double TunnelSolver::getSupplyFunction(double energy)
{
	double T = this->temperature;
	double EfTunnelFrom = this->fermiEnergyTunnelFrom;
	double EfTunnelTo = this->fermiEnergyTunnelTo;

	double kB = SctmPhys::BoltzmanConstant;
	double q = SctmPhys::ElementaryCharge;
	double integralTunnelFrom = 1;
	double integralTunnelTo = 1;

	if (energy - EfTunnelFrom > 5 * kB * T / q)
		integralTunnelFrom =  kB * T * SctmMath::exp(- q * (energy - EfTunnelFrom) / kB / T);
	else
		integralTunnelFrom = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelFrom) / kB / T));

	if (energy - EfTunnelTo > 5 * kB * T / q)
		integralTunnelTo = kB * T * SctmMath::exp(- q * (energy - EfTunnelTo) / kB / T);
	else
		integralTunnelTo = kB * T * SctmMath::ln(1 + SctmMath::exp(- q * (energy - EfTunnelTo) / kB / T));
	
	double supply = integralTunnelFrom - integralTunnelTo;
	return supply;
}

double TunnelSolver::getTransmissionCoefficient(double energy)
{
	double tunnelFactor = 1;
	// The transmission coefficient is proportional to the integral used in the calculation.
	// The difference/factor is used to describe the reflection of the tunneling wave.
	// Here we assume that the reflection is tiny in this problem.

	double m0 = SctmPhys::ElectronMass;
	double q = SctmPhys::ElementaryCharge;
	double hbar = SctmPhys::hbar;
	double cm_in_m = SctmPhys::cm_in_m;
	double integral = 0;
	//direct tunneling and Fowler-Nordheim tunneling are included.
	//the international unit (I.U.) is used in calculating TC
	for (vector<double>::size_type ix = 0; ix != this->deltaX.size(); ++ ix)
	{
		if (cbegde.at(ix) >= energy)
		{
			integral += SctmMath::sqrt(2 * emass.at(ix) * m0 * q * (cbegde.at(ix) - energy)) 
						* deltaX.at(ix) * cm_in_m;
		}
	}
	double tc = tunnelFactor * SctmMath::exp(-2 / hbar * integral);

	return tc;
}

void TunnelSolver::calcDTFNtunneling()
{
	// the unit of calculated current density is [A/m^2]
	// Emin and Emax are processed with q. (divided by q)
	double DTFNdensity = 0; // in [A/m^2]
	double Emin = cbedgeTunnelFrom; // in [eV]
	double Emax = cbegde.front(); // in [eV]

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double h = SctmPhys::h;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double TC = 0;
	double supply = 0;

	double currEnergy = Emin;
	while (currEnergy <= Emax)
	{
		TC = getTransmissionCoefficient(currEnergy);
		supply = getSupplyFunction(currEnergy);
		//mass is normalized using m0, E is normalized using q
		DTFNdensity += 4 * pi * this->effTunnelMass * m0 * q / h / h / h
			* TC
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/cm^2]
	this->currentDensity += DTFNdensity * per_m2_in_per_cm2;
}

void TunnelSolver::calcThermalEmission()
{
	double TEdensity = 0; // in [A/m^2]
	double Emin = cbegde.front(); // in [eV]
	double Emax = Emin + 5; // only calculate the energy within 5eV large than the smallest energy to surpass the barrier

	double T = this->temperature;
	double m0 = SctmPhys::ElectronMass;
	double pi = SctmMath::PI;
	double h = SctmPhys::h;
	double q = SctmPhys::ElementaryCharge;
	double dE = SctmPhys::BoltzmanConstant * T / q; // dE = kT/q, here E is normalized using q
	double per_m2_in_per_cm2 = SctmPhys::per_sqr_m_in_per_sqr_cm;

	double supply = 0;
	double currEnergy = Emin;
	//mass is normalized using m0, E is normalized using q
	while (currEnergy <= Emax)
	{
		supply = getSupplyFunction(currEnergy);
		TEdensity += 4 * pi * this->effTunnelMass * m0 * q / h / h / h
			* 1  //for transmission coefficient
			* supply
			* q * dE;
		currEnergy += dE;
	}
	//the unit of currentDensity is [A/cm^2]
	this->currentDensity += TEdensity * per_m2_in_per_cm2;
}

TunnelSolver::TunnelSolver()
{
	this->currentDensity = 0;
}

double TunnelSolver::GetCurrentDensity()
{
	return this->currentDensity;
}

void TunnelSolver::SolveTunneling()
{
	this->currentDensity = 0;
	calcDTFNtunneling();
	calcThermalEmission();
	//SctmUtils::UtilsDebug.PrintValue(this->currentDensity);
}

void SubsToGateEletronTunnel::PrepareProblem(FDVertex *startVertex)
{
	// IMPORTANT! the parameters are in normalization values. They must be converted !
	// the tunneling direction is north
	double val = 0;
	this->areaFactor = startVertex->EastLength / 2 + startVertex->WestLength / 2;
	FDVertex *currentVertex = startVertex;
	while (currentVertex != NULL)
	{
		this->deltaX.push_back((currentVertex->NorthLength + currentVertex->SouthLength) / 2);
		
		//the method to get physical property is changed
		//this->cbegde.push_back(currentVertex->Phys.conductionBandEnergy);
		val = currentVertex->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy);
		this->cbegde.push_back(val);
		//this->emass.push_back(currentVertex->Phys->electronMass);
		val = currentVertex->Phys->GetPhysPrpty(PhysProperty::eMass);
		this->emass.push_back(val);
		
		currentVertex = currentVertex->NorthVertex; // move the vertex north
	}
}

TunnelTest::TunnelTest()
{
	this->oxideEmass = 0;
	this->siliconBandEdge = 0;
}

void TunnelTest::PrepareProblem(FDVertex *startVertex)
{
	//clear existing vectors and prepare for the following calculation
	this->cbegde.clear();
	this->emass.clear();
	this->deltaX.clear();

	//double oxideThickness = 3; // in [nm]
	double oxide_thickness = this->oxideThickness;

	int gridNumber = 100;
	this->temperature = 300;
	this->effTunnelMass = 0.5; // in m0

	double bandedgeDifference = 3.15; // in [eV], Silicon and oxide band edge difference
	
	//double siliconBandegde = -0.112; // in [eV]
	double silicon_bandedge = this->siliconBandEdge;

	//double elecField = 12.9e6; // in [V/cm]
	double elec_field = this->elecField;

	//double gateVoltage = 4; // in [V]
	double gate_voltage = this->gateVoltage;

	this->cbedgeTunnelFrom = silicon_bandedge;
	this->fermiEnergyTunnelFrom = 0;
	this->fermiEnergyTunnelTo = -gate_voltage;

	double deltax = 0;
	deltax = oxide_thickness / gridNumber * SctmPhys::nm_in_cm;

	double currentX = 0;
	double val = 0;
	for (size_t iVert = 0; iVert != gridNumber + 1; ++iVert)
	{
		currentX = iVert * deltax;

		val = (iVert == 0 || iVert == gridNumber) ? deltax / 2 : deltax;
		this->deltaX.push_back(val);

		val = this->cbedgeTunnelFrom + bandedgeDifference - currentX * elec_field;
		this->cbegde.push_back(val);

		this->emass.push_back(this->oxideEmass);
	}
}

void TunnelTest::SolveParamterSet()
{
	vector<vector<double>> currentMatrix;
	vector<double> oxideEmassSet;
	vector<double> siliconBandEdgeSet;

	this->gateVoltage = 4; // [V]
	this->elecField = 12.9e6; // [V/cm]
	this->oxideThickness = 3; // [nm]

	double val = 0;
	for (int ix = 30; ix <= 210; ix += 10) // totally 19 values
	{
		val = - ix * 0.001;
		siliconBandEdgeSet.push_back(val);
	}

	for (int ix = 38; ix <= 58; ix += 1) // totally 21 values
	{
		val = ix * 0.01;
		oxideEmassSet.push_back(val);
	}

	for (size_t imass = 0 ; imass != oxideEmassSet.size(); ++imass)
	{
		this->oxideEmass = oxideEmassSet.at(imass);
		vector<double> v;
		for (size_t iband = 0 ; iband != siliconBandEdgeSet.size(); ++iband)
		{
			this->siliconBandEdge = siliconBandEdgeSet.at(iband);
			this->PrepareProblem(NULL);
			this->SolveTunneling();
			v.push_back(this->currentDensity);
		}
		currentMatrix.push_back(v);
	}
	SctmFileStream writeFile = SctmFileStream("C:\\Users\\Lunzhy\\Desktop\\TunnelTest.txt", SctmFileStream::Write);
	writeFile.Write2DVectorForOrigin(oxideEmassSet, siliconBandEdgeSet, currentMatrix, "Tunneling current -- x: oxide emass -- y: silicon band edge");
}

void TunnelTest::SolveCalibrate()
{
	this->oxideEmass = 0.42; // [m0]

	vector<double> tunnelCurrent;
	vector<double> voltageSet;

	vector<double> cbedges;
	vector<double> elecFields;

	for (int ix = 1; ix < 11; ++ix)
	{
		voltageSet.push_back(ix);
	}

	// read-in from the parameter file
	SctmFileStream writeFile = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\TunnelCalibrate\\result.txt", SctmFileStream::Write);
	
	//calculation of tunneling current for 2nm
	this->oxideThickness = 2;
	SctmFileStream readFile_2nm = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\TunnelCalibrate\\2nm.txt", SctmFileStream::Read);
	readFile_2nm.ReadTunnelParameter(cbedges, elecFields);
	tunnelCurrent.clear();
	for (int ix = 0; ix != 10; ++ix)
	{
		this->gateVoltage = ix + 1;
		this->siliconBandEdge = cbedges.at(ix);
		this->elecField = elecFields.at(ix) * 1e6;
		this->PrepareProblem(NULL);
		this->SolveTunneling();
		tunnelCurrent.push_back(this->currentDensity);
	}
	writeFile.WriteVector(tunnelCurrent, "2nm current");

	//calculation of tunneling current for 3nm
	this->oxideThickness = 3;
	SctmFileStream readFile_3nm = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\TunnelCalibrate\\3nm.txt", SctmFileStream::Read);
	readFile_3nm.ReadTunnelParameter(cbedges, elecFields);
	tunnelCurrent.clear();
	for (int ix = 0; ix != 10; ++ix)
	{
		this->gateVoltage = ix + 1;
		this->siliconBandEdge = cbedges.at(ix);
		this->elecField = elecFields.at(ix) * 1e6;
		this->PrepareProblem(NULL);
		this->SolveTunneling();
		tunnelCurrent.push_back(this->currentDensity);
	}
	writeFile.WriteVector(tunnelCurrent, "3nm current");
	
	//calculation of tunneling current for 4nm
	this->oxideThickness = 4;
	SctmFileStream readFile_4nm = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\TunnelCalibrate\\4nm.txt", SctmFileStream::Read);
	readFile_4nm.ReadTunnelParameter(cbedges, elecFields);
	tunnelCurrent.clear();
	for (int ix = 0; ix != 10; ++ix)
	{
		this->gateVoltage = ix + 1;
		this->siliconBandEdge = cbedges.at(ix);
		this->elecField = elecFields.at(ix) * 1e6;
		this->PrepareProblem(NULL);
		this->SolveTunneling();
		tunnelCurrent.push_back(this->currentDensity);
	}
	writeFile.WriteVector(tunnelCurrent, "4nm current");

	//calculation of tunneling current for 5nm
	this->oxideThickness = 5;
	SctmFileStream readFile_5nm = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\TunnelCalibrate\\5nm.txt", SctmFileStream::Read);
	readFile_5nm.ReadTunnelParameter(cbedges, elecFields);
	tunnelCurrent.clear();
	for (int ix = 0; ix != 10; ++ix)
	{
		this->gateVoltage = ix + 1;
		this->siliconBandEdge = cbedges.at(ix);
		this->elecField = elecFields.at(ix) * 1e6;
		this->PrepareProblem(NULL);
		this->SolveTunneling();
		tunnelCurrent.push_back(this->currentDensity);
	}
	writeFile.WriteVector(tunnelCurrent, "5nm current");
}
