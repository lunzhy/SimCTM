/**
* @file Material.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-24   16:31
* @note
* @todo
*/

#include <iostream>
#include "Material.h"
#include "SctmUtils.h"
#include "Normalization.h"

namespace MaterialDB
{
	Material Silicon = Material("Silicon");
	Material SiO2 = Material("SiO2");
	Material Si3N4 = Material("Si3N4");
	
	void SetMaterials()
	{
		using SctmUtils::Normalization;
		Normalization norm = Normalization();

		//Silicon = Material("Silicon");
		Silicon.Bandgap(norm.PushEnergy(1.12));
		Silicon.DielectricConstant(11.9);
		Silicon.ElectronAffinity(norm.PushEnergy(4.05));
		Silicon.ElectronMass(1); // need to be revised
		Silicon.HoleMass(1);
		Silicon.ElectronMobility(norm.PushMobility(1350));

		//SiO2 = Material("SiO2");
		SiO2.Bandgap(norm.PushEnergy(9.4));
		SiO2.DielectricConstant(3.9);
		SiO2.ElectronAffinity(norm.PushEnergy(0.9));
		SiO2.ElectronMass(0.42); // need to be revised
		SiO2.HoleMass(1);
		SiO2.ElectronMobility(0);

		//Si3N4 = Material("Si3N4");
		Si3N4.Bandgap(norm.PushEnergy(5.0));
		Si3N4.DielectricConstant(7.5);
		Si3N4.ElectronAffinity(norm.PushEnergy(1.9));
		Si3N4.ElectronMass(1); // need to be revised
		Si3N4.HoleMass(1);
		Si3N4.ElectronMobility(norm.PushMobility(0.1));
		Si3N4.ElecTrapEnergyFromCB(norm.PushEnergy(1.5));
		Si3N4.ElecTrapXSection(norm.PushArea(1e-15));
	}

	double GetMatPrpty(Material *theMaterial, MatProperty::Name prptyName)
	{
		double ret = 0;
		
		switch (prptyName)
		{
		case MatProperty::Mat_Bandgap:
			ret = theMaterial->Bandgap();
			break;
		case  MatProperty::Mat_ElectronAffinity:
			ret = theMaterial->ElectronAffinity();
			break;
		case  MatProperty::Mat_DielectricConstant:
			ret = theMaterial->DielectricConstant();
			break;
		case  MatProperty::Mat_ElectronMass:
			ret = theMaterial->ElectronMass();
			break;
		case MatProperty::Mat_HoleMass:
			ret = theMaterial->HoleMass();
			break;
		case MatProperty::Mat_ElectronDiffusion:
			ret = theMaterial->ElectronDiffusion();
			break;
		case MatProperty::Mat_HoleDiffusion:
			ret = theMaterial->HoleDiffusion();
			break;
		case MatProperty::Mat_ElectronMobility:
			ret = theMaterial->ElectronMobility();
			break;
		case MatProperty::Mat_HoleMobility:
			ret = theMaterial->HoleMobility();
			break;
		case MatProperty::Mat_ElecTrapEnergyFromCB:
			ret = theMaterial->ElecTrapEnergyFromCB();
			break;
		case MatProperty::Mat_ElecTrapXSection:
			ret = theMaterial->ElecTrapXSection();
			break;
		default:
			// use SCTM_CHECK for non-existed property
			SCTM_ASSERT(SCTM_ERROR, 10002);
		}

		return ret;
	}

}