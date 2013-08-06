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
	const double ElementaryCharge			= 1.602176487e-19;		// in [C]
	const double PlanckConstant				= 6.62606896e-34;		// in [J.s]
	const double BoltzmanConstant			= 1.3806504e-23;		// in [J/K]
	const double ElectronMass				= 9.10938215e-31;		// in [kg]
	const double VacuumDielectricConstant	= 8.854187817e-12;		// in [F/m]
	const double RoomTemperature			= 300;					// in [K]
	
	// unit conversion
	const double cm_in_m					= 1e-2;
	const double nm_in_cm					= 1e-7;
	const double per_sqr_m_in_per_sqr_cm	= 1e-4; // per square meter in per square centimeter 1/m^2 = 1e-4 /cm^2
	
	// simplified physical constant
	extern const double &k0;
	extern const double &h;
	extern const double &eps;
	extern const double &T0;
	extern const double &q;
	const double hbar = h / 2 / SctmMath::PI;

	
	extern double ReferencePotential;

	void SetPhysConstant();

	/// @brief This class is a data structure to store the physical parameters of the specified vertex
	///
	///
	class PhysProperty
	{
	public:
		/// @brief
		///
		///
		enum Name
		{
			ElectrostaticPotential, ///<
			ConductionBandEnergy, ///<
			ValenceBandEnergy, ///<
			ElectronAffinity, ///<
			ElectronMass, ///<
			HoleMass, ///<
			Bandgap ///<
		};

		void SetPhyPrpty(Name prptyName, double prptyValue);
		double GetPhyPrpty(Name prptyName) const;
	private:
		//TODO : initialize these values when constructing the object. Then we can judge the value when they are used
		double bandgap;
		double electrostaticPotential;
		double conductionBandEnergy; ///< i.e. conduction band edge
		double valenceBandEnergy; ///< i.e. valence band edge
		double electronAffinity;
		double electronMass;
		double holeMass;
	};
}