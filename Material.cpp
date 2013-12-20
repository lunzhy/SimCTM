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
#include "SctmPhys.h"
#include "SctmMath.h"

using SctmUtils::Normalization;

namespace MaterialDB
{
	std::map<Mat::Name, Material*> MaterialMap;
	Material Silicon_material = Material(Mat::Silicon);
	Material SiO2_material = Material(Mat::SiO2);
	Material Si3N4_material = Material(Mat::Si3N4);
	Material Al2O3_material = Material(Mat::Al2O3);
	Material HfO2_material = Material(Mat::HfO2);
	
	void SetMaterials()
	{
		//Silicon
		Silicon_material.Bandgap(1.12);
		Silicon_material.DielectricConstant(11.9);
		Silicon_material.ElectronAffinity(4.05);
		Silicon_material.ElectronMass(1.08); // need to be revised
		Silicon_material.HoleMass(1);
		Silicon_material.ElectronMobility(1350);

		//SiO2
		SiO2_material.Bandgap(9.4);
		SiO2_material.DielectricConstant(3.9);
		SiO2_material.ElectronAffinity(0.9);
		SiO2_material.ElectronMass(0.32); // need to be revised
		SiO2_material.HoleMass(1);
		SiO2_material.ElectronMobility(0);

		//Al2O3
		Al2O3_material.Bandgap(8.8);
		Al2O3_material.DielectricConstant(9.0);
		Al2O3_material.ElectronAffinity(1.25);
		Al2O3_material.ElectronMass(0.2);

		//Si3N4
		Si3N4_material.Bandgap(5.0);
		Si3N4_material.DielectricConstant(7.5);
		Si3N4_material.ElectronAffinity(1.9);
		Si3N4_material.ElectronMass(0.42); // need to be revised
		Si3N4_material.HoleMass(1);
		Si3N4_material.ElectronMobility(0.1);
		Si3N4_material.ElecTrapEnergyFromCB(0.7);
		Si3N4_material.ElecTrapXSection(1e-14);

		//HfO2
		HfO2_material.Bandgap(5.9);
		HfO2_material.DielectricConstant(20.0);
		HfO2_material.ElectronAffinity(2.05);
		HfO2_material.ElectronMass(0.2);
		HfO2_material.HoleMass(1); // need to be revised
		HfO2_material.ElectronMobility(0.01);
		HfO2_material.ElecTrapEnergyFromCB(0.7);
		HfO2_material.ElecTrapXSection(1e-14);

		MaterialMap[Mat::Silicon] = &Silicon_material;
		MaterialMap[Mat::SiO2] = &SiO2_material;
		MaterialMap[Mat::Si3N4] = &Si3N4_material;
		MaterialMap[Mat::Al2O3] = &Al2O3_material;
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

	Material::Material(Mat::Name _name) :name(_name)
	{
		using SctmUtils::SctmGlobalControl;
		temperature = SctmGlobalControl::Get().Temperature;
	}

	double Material::DielectricConstant() const
	{
		return dielectricConstant;
	}

	void Material::DielectricConstant(double val)
	{
		dielectricConstant = val;
	}

	double Material::Bandgap() const
	{
		return bandgap;
	}

	void Material::Bandgap(double val)
	{
		Normalization norm = Normalization(temperature);
		bandgap = norm.PushEnergy(val);
	}

	double Material::ElectronAffinity() const
	{
		return electronAffinity;
	}

	void Material::ElectronAffinity(double val)
	{
		Normalization norm = Normalization(temperature);
		electronAffinity = norm.PushEnergy(val);
	}

	double Material::ElectronMass() const
	{
		return electronMass;
	}

	void Material::ElectronMass(double val)
	{
		electronMass = val;
	}

	double Material::HoleMass() const
	{
		return holeMass;
	}

	void Material::HoleMass(double val)
	{
		holeMass = val;
	}

	double Material::ElectronDiffusion() const
	{
		return electronDiffusion;
	}

	void Material::ElectronDiffusion(double val)
	{
		Normalization norm = Normalization(temperature);
		electronDiffusion = norm.PushDiffusion(val);
	}

	double Material::HoleDiffusion() const
	{
		return holeDiffusion;
	}

	void Material::HoleDiffusion(double val)
	{
		Normalization norm = Normalization(temperature);
		holeDiffusion = norm.PushDiffusion(val);
	}

	double Material::ElectronMobility() const
	{
		return electronMobility;
	}

	void Material::ElectronMobility(double val)
	{
		Normalization norm = Normalization(temperature);
		electronMobility = norm.PushMobility(val);
	}

	double Material::HoleMobility() const
	{
		return holeMobility;
	}

	void Material::HoleMobility(double val)
	{
		Normalization norm = Normalization(temperature);
		holeMobility = norm.PushMobility(val);
	}

	double Material::ElecTrapXSection() const
	{
		return elecTrapXSection;
	}

	void Material::ElecTrapXSection(double val)
	{
		Normalization norm = Normalization(temperature);
		elecTrapXSection = norm.PushArea(val);
	}

	double Material::ElecTrapEnergyFromCB() const
	{
		return elecTrapEnergyFromCB;
	}

	void Material::ElecTrapEnergyFromCB(double val)
	{
		Normalization norm = Normalization(temperature);
		elecTrapEnergyFromCB = norm.PushEnergy(val);
	}


}