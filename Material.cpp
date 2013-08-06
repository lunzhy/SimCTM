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
namespace MaterialDB
{
	Material Silicon = Material("Silicon");
	Material SiO2 = Material("SiO2");
	Material Si3N4 = Material("Si3N4");
	
	void SetMaterials()
	{
		//Silicon = Material("Silicon");
		Silicon.Bandgap(1.12);
		Silicon.DielectricConstant(11.9);
		Silicon.ElectronAffinity(4.05);
		Silicon.ElectronMass(1); // need to be revised
		Silicon.HoleMass(1);

		//SiO2 = Material("SiO2");
		SiO2.Bandgap(9.4);
		SiO2.DielectricConstant(3.9);
		SiO2.ElectronAffinity(0.9);
		SiO2.ElectronMass(1); // need to be revised
		SiO2.HoleMass(1);

		//Si3N4 = Material("Si3N4");
		Si3N4.Bandgap(5.0);
		Si3N4.DielectricConstant(7.5);
		Si3N4.ElectronAffinity(1.9);
		Si3N4.ElectronMass(1); // need to be revised
		Si3N4.HoleMass(1);
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
		default:;
			// use SCTM_CHECK for non-existed property
		}

		return ret;
	}

}