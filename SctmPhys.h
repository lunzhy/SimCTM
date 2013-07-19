/**
* @file SctmPhys.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-19   10:36
* @note
* @todo
*/

#include <cmath>
namespace SctmPhys
{
	const double ElementaryCharge			= 1.602176487e-19;
	const double PlanckConstant				= 6.62606896e-34;
	const double BoltzmanConstant			= 1.3806504e-23;
	const double ElectronMass				= 9.10938215e-31;
	const double VacuumDielectricConstant	= 8.854187817e-12;
	const double nm_in_cm					= 1e-7;

	const double &h = PlanckConstant;
	const double hbar = h / 2 / M_PI;
	const double &epsilon = VacuumDielectricConstant;
};