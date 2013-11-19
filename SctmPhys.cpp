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
#include "FDDomain.h"

namespace SctmPhys
{
	const double &k0 = BoltzmanConstant;
	const double &h = PlanckConstant;
	const double &eps = VacuumDielectricConstant;
	const double &T0 = RoomTemperature;
	const double &q = ElementaryCharge;
	const double &ni = IntrinsicConcentration;
	
	double ReferencePotential;

	PhysProperty::PhysProperty()
	{
		bandgap = 0;
		electrostaticPotential = 0;
		//conductionBandEnergy = 0;
		//valenceBandEnergy = 0;
		electronAffinity = 0;
		e_mass = 0;
		h_mass = 0;
		//netCharge = 0;
		e_density = 0;
		controlArea = 0;
	}

	void PhysProperty::SetPhysPrpty(Name prptyName, double prptyValue)
	{
		switch (prptyName)
		{
		case ElectrostaticPotential:
			electrostaticPotential = prptyValue;
			break;
		case eMass:
			e_mass = prptyValue;
			break;
		case hMass:
			h_mass = prptyValue;
			break;
		case ElectronAffinity:
			electronAffinity = prptyValue;
			break;
		case Bandgap:
			bandgap = prptyValue;
			break;
		case NetCharge:
			netCharge = prptyValue;
			break;
		case eMobility:
			e_mobility = prptyValue;
			break;
		case eDensity:
			e_density = prptyValue;
			break;
		default:
			SCTM_ASSERT(SCTM_ERROR, 10019);
		}
	}

	double PhysProperty::GetPhysPrpty(Name prptyName) const
	{
		static double RefPotential = SctmPhys::ReferencePotential;

		double energy = 0;
		double affinity = 0;
		double pot = 0;
		double ret = 0;

		switch (prptyName)
		{
		case ElectrostaticPotential:
			ret = electrostaticPotential;
			break;
		case ConductionBandEnergy:
			//Ec = -X-q(phi-phiRef)
			pot = this->electrostaticPotential;
			affinity = this->electronAffinity;
			energy = -affinity - (pot-RefPotential);
			ret = energy;
			break;
		case ValenceBandEnergy:
			//Ev = Ev-Eg = -X-q(phi-phiRef)-Eg
			ret =  GetPhysPrpty(ConductionBandEnergy) - bandgap;
			break;
		case eMass:
			ret = e_mass;
			break;
		case hMass:
			ret = h_mass;
			break;
		case ElectronAffinity:
			ret = electronAffinity;
			break;
		case Bandgap:
			ret = bandgap;
			break;
		case NetCharge:
			//there should be four parts
			ret = e_density;
			break;
		case eMobility:
			ret = e_mobility;
			break;
		case eDensity:
			ret = e_density;
			break;
		case DensityControlArea:
			ret = controlArea;
			break;
		default:
			// use SCTM_ASSERT for non-existed property
			SCTM_ASSERT(SCTM_ERROR, 10001);
		}

		return ret;
	}

	void PhysProperty::FillVertexPhysUsingMatPropty(FDVertex *vertex, PhysProperty::Name vertexPhys,
		MaterialDB::MatProperty::Name matPrpty)
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
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->NorthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->SoutheastElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->SouthwestElem;
		tot += ( currElem != NULL ) ? currElem->Area : 0;
		sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		SCTM_ASSERT(tot>=0, 10004);
		physValue = sum / tot;

		vertex->Phys->SetPhysPrpty(vertexPhys, physValue);
	}		

	void PhysProperty::FillVertexPhysUsingMatPropty(FDVertex *vertex, PhysProperty::Name vertexPhys, 
		MaterialDB::MatProperty::Name matPrpty, FDRegion::RegionType rType)
	{
		//TODO : need to solve the problem of mutual including.
		//the problem is solved but with some unknowns.
		double tot = 0; //total area
		double sum = 0; //sum corresponds to the integral
		double physValue = 0; //the final physical value related to vertex

		using MaterialDB::GetMatPrpty;
		FDElement *currElem = NULL;
		currElem = vertex->NortheastElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->NorthwestElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->SoutheastElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		currElem = vertex->SouthwestElem;
		tot += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? currElem->Area : 0;
		sum += (( currElem != NULL ) && (currElem->Region->Type == rType)) ? GetMatPrpty(currElem->Region->Mat, matPrpty) * currElem->Area : 0;

		if (tot == 0)
		{
			physValue = 0;
		}
		else
		{
			physValue = sum / tot;
		}
		//SCTM_ASSERT(tot>=0, 10004);
		vertex->Phys->SetPhysPrpty(vertexPhys, physValue);
	}

	void PhysProperty::CalculateDensityControlArea(FDVertex *vertex)
	{
		double area = 0;
		FDElement *currElem = NULL;

		currElem = vertex->NorthwestElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertex->NortheastElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertex->SouthwestElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		currElem = vertex->SoutheastElem;
		if ( (currElem != NULL) && (currElem->Region->Type == FDRegion::Trapping) )
		{
			area += 0.25 * currElem->Area;
		}
		this->controlArea = area;
	}

	void PhysProperty::UpdateValue(Name prptyName, double val)
	{
		//only several properties are updated
		SetPhysPrpty(prptyName, val);
	}

	void SetPhysConstant()
	{
		double mp = GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_HoleMass);
		double mn = GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronMass);
		double bandgap = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_Bandgap);
		//the value of reference potential is
		//phi(electron affinity) + Eg/2 + 3/4*kT/q*ln(mp/mn)
		SctmPhys::ReferencePotential = MaterialDB::GetMatPrpty(&MaterialDB::Silicon, MaterialDB::MatProperty::Mat_ElectronAffinity)
										+ bandgap / 2 - 3/4 * k0 * T0 / q * SctmMath::ln( mp / mn);
	}

}