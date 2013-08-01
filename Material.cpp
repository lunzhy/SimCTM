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

		//SiO2 = Material("SiO2");
		SiO2.Bandgap(9.4);
		SiO2.DielectricConstant(3.9);
		SiO2.ElectronAffinity(0.9);

		//Si3N4 = Material("Si3N4");
		Si3N4.Bandgap(5.0);
		Si3N4.DielectricConstant(7.5);
		Si3N4.ElectronAffinity(1.9);

		std::cout << Si3N4.Bandgap() << std::endl;
	}
}