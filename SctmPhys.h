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
#pragma once
#include <cmath>
#include "SctmMath.h"
/// @brief This namespace contains all the physics used in the simulation 
///
/// The common physical parameters are defined here. And other classes and structs related to
/// the physics problems are defined here.
namespace SctmPhys
{
	// Physical constant
	const double ElementaryCharge			= 1.602176487e-19; // in [C]
	const double PlanckConstant				= 6.62606896e-34; // in [J.s]
	const double BoltzmanConstant			= 1.3806504e-23; // in [J/K]
	const double ElectronMass				= 9.10938215e-31; // in [kg]
	const double VacuumDielectricConstant	= 8.854187817e-12; // in [F/m]
	
	// unit conversion
	const double cm_in_m					= 1e-2;
	const double nm_in_cm					= 1e-7;
	const double per_sqr_m_in_per_sqr_cm	= 1e-4; // per square meter in per square centimeter 1/m^2 = 1e-4 /cm^2
	
	// simplified physical constant
	const double h = PlanckConstant;
	const double hbar = h / 2 / SctmMath::PI;
	const double eps = VacuumDielectricConstant;

	/// @brief This struct is a data structure to store the physical parameters of the specified vertex
	///
	///
	struct PhysProperty
	{
		double ElectrostaticPotential;
		double ConductionBandEnergy; ///< i.e. conduction band edge
		double ValenceBandEnergy; ///< i.e. valence band edge
		double ElectronMass;
		double HoleMass;
	};
};