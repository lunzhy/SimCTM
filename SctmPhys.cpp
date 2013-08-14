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
#include "DomainDetails.h"

namespace SctmPhys
{
	const double &k0 = BoltzmanConstant;
	const double &h = PlanckConstant;
	const double &eps = VacuumDielectricConstant;
	const double &T0 = RoomTemperature;
	const double &q = ElementaryCharge;
	
	double ReferencePotential;

	void PhysProperty::SetPhysPrpty(Name prptyName, double prptyValue)
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

	double PhysProperty::GetPhysPrpty(Name prptyName) const
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
			// use SCTM_ASSERT for non-existed property
			SCTM_ASSERT(true, 1);
		}

		return ret;
	}

	void PhysProperty::FillVertexPhysUsingMatPropty(FDVertex *vertex, PhysProperty::Name vertexPhys, MaterialDB::MatProperty::Name matPrpty)
	{
		//TODO : need to solve the problem of mutual including.
		//the problem is solved but with some unknowns.
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double physValue = 0; //the final physical value related to vertex

		using MaterialDB::GetMatPrpty;
		FDElement *currElem = NULL;
		currElem = vertex->NortheastElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrpty) * currElem->Area : 0;

		currElem = vertex->NorthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrpty) * currElem->Area : 0;

		currElem = vertex->SoutheastElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrpty) * currElem->Area : 0;

		currElem = vertex->SouthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrpty) * currElem->Area : 0;

		SCTM_ASSERT(tot=0, 4);
		physValue = sum / tot;

		vertex->Phys.SetPhysPrpty(vertexPhys, physValue);
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