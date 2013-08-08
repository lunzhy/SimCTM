/**
* @file SctmPhys.cpp
* @brief This file contains the detailed implementation of the physics problems in the simulation
*
*
*
* @author
* @version 
* @date 2013-8-5   15:22
* @note
* @todo
*/

#pragma once
#include "SctmPhys.h"
#include "Material.h"
#include "SctmUtils.h"

namespace SctmPhys
{
	const double &k0 = BoltzmanConstant;
	const double &h = PlanckConstant;
	const double &eps = VacuumDielectricConstant;
	const double &T0 = RoomTemperature;
	const double &q = ElementaryCharge;
	
	double ReferencePotential;

	void PhysProperty::SetPhyPrpty(Name prptyName, double prptyValue)
	{
		switch (prptyName)
		{
		case ElectrostaticPotential:
			electrostaticPotential = prptyValue;
			break;
		case ConductionBandEnergy:
			conductionBandEnergy = prptyValue;
			break;
		case ValenceBandEnergy:
			valenceBandEnergy = prptyValue;
			break;
		case ElectronMass:
			electronMass = prptyValue;
			break;
		case HoleMass:
			holeMass = prptyValue;
			break;
		case ElectronAffinity:
			electronAffinity = prptyValue;
			break;
		case Bandgap:
			bandgap = prptyValue;
			break;
		}
	}

	double PhysProperty::GetPhyPrpty(Name prptyName) const
	{
		double ret;

		switch (prptyName)
		{
		case ElectrostaticPotential:
			ret = electrostaticPotential;
			break;
		case ConductionBandEnergy:
			ret = conductionBandEnergy;
			break;
		case ValenceBandEnergy:
			ret = valenceBandEnergy;
			break;
		case ElectronMass:
			ret = electronMass;
			break;
		case HoleMass:
			ret = holeMass;
			break;
		case ElectronAffinity:
			ret = electronAffinity;
			break;
		case Bandgap:
			ret = bandgap;
			break;
		default:
			SCTM_ASSERT(1==1, 0);
			// use SCTM_ASSERT for non-existed property
		}

		return ret;
	}

	void SetPhysConstant()
	{
		double mp = GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_HoleMass);
		double mn = GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronMass);
		double bandgap = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_Bandgap);
		//the value of reference potential is
		//electron affinity + Eg/2 + 3/4*kT/q*ln(mp/mn)
		ReferencePotential = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronAffinity)
			+ bandgap / 2 - 3/4 * k0 * T0 / q * SctmMath::ln( mp / mn);
	}

}