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
	//std::map<Mat::Name, Material*> MaterialMap;
	//Material MaterialMap(Mat::Silicon) = Material(Mat::Silicon);
	//Material SiO2_material = Material(Mat::SiO2);
	//Material Si3N4_material = Material(Mat::Si3N4);
	//Material Al2O3_material = Material(Mat::Al2O3);
	//Material HfO2_material = Material(Mat::HfO2);
	
	Material* MaterialMap(Mat::Name matname)
	{
		static Material *Silicon_material = new Material(Mat::Silicon);
		static Material *SiO2_material = new Material(Mat::SiO2);
		static Material *Si3N4_material = new Material(Mat::Si3N4);
		static Material *Al2O3_material = new Material(Mat::Al2O3);
		static Material *HfO2_material = new Material(Mat::HfO2);

		if (matname == Mat::Silicon)
		{
			return Silicon_material;
		}
		if (matname == Mat::SiO2)
		{
			return SiO2_material;
		}
		if (matname == Mat::Si3N4)
		{
			return Si3N4_material;
		}
		if (matname == Mat::Al2O3)
		{
			return Al2O3_material;
		}
		if (matname == Mat::HfO2)
		{
			return HfO2_material;
		}
		return NULL;
	}

	void SetMaterials_Directly()
	{
		//Silicon
		MaterialMap(Mat::Silicon)->Bandgap(1.12);
		MaterialMap(Mat::Silicon)->Bandgap(1.12);
		MaterialMap(Mat::Silicon)->DielectricConstant(11.9);
		MaterialMap(Mat::Silicon)->ElectronAffinity(4.05);
		MaterialMap(Mat::Silicon)->ElectronMass(1.08); // need to be revised
		MaterialMap(Mat::Silicon)->HoleMass(0.59);
		MaterialMap(Mat::Silicon)->ElectronMobility(1350);

		//SiO2
		MaterialMap(Mat::SiO2)->Bandgap(9.4);
		MaterialMap(Mat::SiO2)->DielectricConstant(3.9);
		MaterialMap(Mat::SiO2)->ElectronAffinity(0.9);
		MaterialMap(Mat::SiO2)->ElectronMass(0.32); // need to be revised
		MaterialMap(Mat::SiO2)->HoleMass(1);
		MaterialMap(Mat::SiO2)->ElectronMobility(0);

		//Al2O3
		MaterialMap(Mat::Al2O3)->Bandgap(8.8);
		MaterialMap(Mat::Al2O3)->DielectricConstant(9.0);
		MaterialMap(Mat::Al2O3)->ElectronAffinity(1.25);
		MaterialMap(Mat::Al2O3)->ElectronMass(0.2);

		//Si3N4
		MaterialMap(Mat::Si3N4)->Bandgap(5.0);
		MaterialMap(Mat::Si3N4)->DielectricConstant(7.5);
		MaterialMap(Mat::Si3N4)->ElectronAffinity(1.9);
		MaterialMap(Mat::Si3N4)->ElectronMass(0.42); // need to be revised
		MaterialMap(Mat::Si3N4)->HoleMass(1);
		MaterialMap(Mat::Si3N4)->ElectronMobility(0.1);
		MaterialMap(Mat::Si3N4)->ElecTrapEnergyFromCB(1.2);
		MaterialMap(Mat::Si3N4)->ElecTrapXSection(1e-14);

		//HfO2
		MaterialMap(Mat::HfO2)->Bandgap(5.9);
		MaterialMap(Mat::HfO2)->DielectricConstant(20.0);
		MaterialMap(Mat::HfO2)->ElectronAffinity(2.05);
		MaterialMap(Mat::HfO2)->ElectronMass(0.2);
		MaterialMap(Mat::HfO2)->HoleMass(1); // need to be revised
		MaterialMap(Mat::HfO2)->ElectronMobility(0.01);
		MaterialMap(Mat::HfO2)->ElecTrapEnergyFromCB(0.7);
		MaterialMap(Mat::HfO2)->ElecTrapXSection(1e-14);


		//MaterialMap[Mat::Silicon] = &Silicon_material;
		//MaterialMap[Mat::SiO2] = &SiO2_material;
		//MaterialMap[Mat::Si3N4] = &Si3N4_material;
		//MaterialMap[Mat::Al2O3] = &Al2O3_material;
		//MaterialMap[Mat::HfO2] = &MaterialMap(Mat::HfO2);
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
		defaultPara:
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
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_bandgap);
		MaterialMap(Mat::Silicon)->Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_dielectricConstant);
		MaterialMap(Mat::Silicon)->DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_electronAffinity);
		MaterialMap(Mat::Silicon)->ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_eMass);
		MaterialMap(Mat::Silicon)->ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_hMass);
		MaterialMap(Mat::Silicon)->HoleMass(dynamic_cast<Param<double> *>(parBase)->Value());
		
		//MaterialMap(Mat::Silicon)->HoleMass(1);
		
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si_eMobility);
		MaterialMap(Mat::Silicon)->ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		//SiO2
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::SiO2_bandgap);
		MaterialMap(Mat::SiO2)->Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::SiO2_dielectricConstant);
		MaterialMap(Mat::SiO2)->DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::SiO2_electronAffinity);
		MaterialMap(Mat::SiO2)->ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::SiO2_eMass);
		MaterialMap(Mat::SiO2)->ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value()); // need to be revised
		
		//MaterialMap(Mat::SiO2)->HoleMass(1);

		//Al2O3
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Al2O3_bandgap);
		MaterialMap(Mat::Al2O3)->Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Al2O3_dielectricConstant);
		MaterialMap(Mat::Al2O3)->DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Al2O3_electronAffinity);
		MaterialMap(Mat::Al2O3)->ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Al2O3_eMass);
		MaterialMap(Mat::Al2O3)->ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());

		//Si3N4
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_bandgap);
		MaterialMap(Mat::Si3N4)->Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_dielectricConstant);
		MaterialMap(Mat::Si3N4)->DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_electronAffinity);
		MaterialMap(Mat::Si3N4)->ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_eMass);
		MaterialMap(Mat::Si3N4)->ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value()); // need to be revised

		//MaterialMap(Mat::Si3N4)->HoleMass(1);
		
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_eMobility);
		MaterialMap(Mat::Si3N4)->ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_eTrapEnergy);
		MaterialMap(Mat::Si3N4)->ElecTrapEnergyFromCB(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::Si3N4_eXsection);
		MaterialMap(Mat::Si3N4)->ElecTrapXSection(dynamic_cast<Param<double> *>(parBase)->Value());

		//HfO2
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_bandgap);
		MaterialMap(Mat::HfO2)->Bandgap(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_dielectricConstant);
		MaterialMap(Mat::HfO2)->DielectricConstant(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_electronAffinity);
		MaterialMap(Mat::HfO2)->ElectronAffinity(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_eMass);
		MaterialMap(Mat::HfO2)->ElectronMass(dynamic_cast<Param<double> *>(parBase)->Value());
		
		//MaterialMap(Mat::HfO2)->HoleMass(1); // need to be revised
		
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_eMobility);
		MaterialMap(Mat::HfO2)->ElectronMobility(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_eTrapEnergy);
		MaterialMap(Mat::HfO2)->ElecTrapEnergyFromCB(dynamic_cast<Param<double> *>(parBase)->Value());

		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::HfO2_eXsection);
		MaterialMap(Mat::HfO2)->ElecTrapXSection(dynamic_cast<Param<double> *>(parBase)->Value());


		//MaterialMap[Mat::Silicon] = &MaterialMap(Mat::Silicon);
		//MaterialMap[Mat::SiO2] = &SiO2_material;
		//MaterialMap[Mat::Si3N4] = &Si3N4_material;
		//MaterialMap[Mat::Al2O3] = &Al2O3_material;
		//MaterialMap[Mat::HfO2] = &HfO2_material;
	}

}