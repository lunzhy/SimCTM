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
using SctmUtils::SctmParameterParser;
using SctmUtils::ParamBase;
using SctmUtils::Param;

namespace MaterialDB
{
	std::map<Mat::Name, Material*> MaterialMap;
	Material Silicon_material = Material(Mat::Silicon);
	Material SiO2_material = Material(Mat::SiO2);
	Material Si3N4_material = Material(Mat::Si3N4);
	Material Al2O3_material = Material(Mat::Al2O3);
	Material HfO2_material = Material(Mat::HfO2);
	
	void SetMaterials_Directly()
	{
		//Silicon
		Silicon_material.Bandgap(1.12);
		Silicon_material.DielectricConstant(11.9);
		Silicon_material.ElectronAffinity(4.05);
		Silicon_material.ElectronMass(1.08); // need to be revised
		Silicon_material.HoleMass(0.59);
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
		Si3N4_material.ElecTrapEnergyFromCB(1.2);
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
		MaterialMap[Mat::HfO2] = &HfO2_material;
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


	Mat::Name Mat::Parse(const std::string &matStr)
	{
		if (matStr == "Silicon" || matStr == "Si")
		{
			return Mat::Silicon;
		}
		else if (matStr == "SiO2")
		{
			return Mat::SiO2;
		}
		else if (matStr == "Al2O3")
		{
			return Mat::Al2O3;
		}
		else if (matStr == "Si3N4")
		{
			return Mat::Si3N4;
		}
		else if (matStr == "HfO2")
		{
			return Mat::HfO2;
		}
		return Mat::ErrorMaterial;
	}

	void SetMaterial_FromParFile()
	{
		ParamBase *parBase = NULL;
		
		//Silicon
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_bandgap);
		Silicon_material.Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_dielectricConstant);
		Silicon_material.DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_electronAffinity);
		Silicon_material.ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_eMass);
		Silicon_material.ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_hMass);
		Silicon_material.HoleMass(dynamic_cast<Param<double> *>(parBase)->Value());
		
		//Silicon_material.HoleMass(1);
		
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si_eMobility);
		Silicon_material.ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		//SiO2
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::SiO2_bandgap);
		SiO2_material.Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::SiO2_dielectricConstant);
		SiO2_material.DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::SiO2_electronAffinity);
		SiO2_material.ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::SiO2_eMass);
		SiO2_material.ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value()); // need to be revised
		
		//SiO2_material.HoleMass(1);

		//Al2O3
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Al2O3_bandgap);
		Al2O3_material.Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Al2O3_dielectricConstant);
		Al2O3_material.DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Al2O3_electronAffinity);
		Al2O3_material.ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Al2O3_eMass);
		Al2O3_material.ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());

		//Si3N4
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_bandgap);
		Si3N4_material.Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_dielectricConstant);
		Si3N4_material.DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_electronAffinity);
		Si3N4_material.ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_eMass);
		Si3N4_material.ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value()); // need to be revised

		//Si3N4_material.HoleMass(1);
		
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_eMobility);
		Si3N4_material.ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_eTrapEnergy);
		Si3N4_material.ElecTrapEnergyFromCB(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::Si3N4_eXsection);
		Si3N4_material.ElecTrapXSection(dynamic_cast<Param<double> *>(parBase)->Value());

		//HfO2
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_bandgap);
		HfO2_material.Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_dielectricConstant);
		HfO2_material.DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_electronAffinity);
		HfO2_material.ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_eMass);
		HfO2_material.ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());
		
		//HfO2_material.HoleMass(1); // need to be revised
		
		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_eMobility);
		HfO2_material.ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_eTrapEnergy);
		HfO2_material.ElecTrapEnergyFromCB(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::GetInstance().GetPar(SctmParameterParser::HfO2_eXsection);
		HfO2_material.ElecTrapXSection(dynamic_cast<Param<double> *>(parBase)->Value());


		MaterialMap[Mat::Silicon] = &Silicon_material;
		MaterialMap[Mat::SiO2] = &SiO2_material;
		MaterialMap[Mat::Si3N4] = &Si3N4_material;
		MaterialMap[Mat::Al2O3] = &Al2O3_material;
		MaterialMap[Mat::HfO2] = &HfO2_material;
	}

}