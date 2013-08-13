/**
* @file SctmPhys.h 
* @brief This file contains the common physical problems occurred in the simulation
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
#include "Material.h"
#include "DomainDetails.h"
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

	
	extern double ReferencePotential; // the reference potential used in calculating the conduction/valence band energy

	/// @brief SetPhysConstant is used to set the constant parameters which are initialized with the material parameters
	/// 
	/// These physical constant cannot be set to const, because their definitions needs to call the other methods.
	/// 
	/// @pre
	/// @return void
	/// @note
	void SetPhysConstant();

	/// @brief This class is a data structure to store the physical parameters of the specified vertex
	///
	///
	class PhysProperty
	{
	public:
		/// @brief this enum is use to store the name of physical properties linked to each vertex
		///
		/// The name of the physical property is used to get and set related value of the vertex
		enum Name
		{
			ElectrostaticPotential, ///< potential
			ConductionBandEnergy, ///< conduction band energy
			ValenceBandEnergy, ///< valence band energy
			ElectronAffinity, ///< electron affinity
			ElectronMass, ///< electron effective mass
			HoleMass, ///< hole effective mass
			Bandgap ///< bandgap
		};

		/// @brief PhysProperty is the construction method for this class
		/// 
		/// All the physical values are initialized to zero in the construction method.
		/// 
		/// @pre
		/// @return 
		/// @note
		PhysProperty()
		{
			bandgap = 0;
			electrostaticPotential = 0;
			conductionBandEnergy = 0;
			valenceBandEnergy = 0;
			electronAffinity = 0;
			electronMass = 0;
			holeMass = 0;
		}
		/// @brief SetPhyPrpty is used to set the value of physical property related to each vertex
		/// 
		/// The name of the specified property is given in enum name. 
		/// 
		/// @param Name prptyName
		/// @param double prptyValue
		/// @pre
		/// @return void
		/// @note
		void SetPhysPrpty(Name prptyName, double prptyValue);
		/// @brief GetPhyPrpty is used to get the value of physical property related to the vertex
		/// 
		/// This is const method and no value should be changed in this method.
		/// 
		/// @param Name prptyName
		/// @pre
		/// @return double
		/// @note
		double GetPhysPrpty(Name prptyName) const;


		void FillVertexPhysUsingMatPropty(PhysProperty::Name vertexPhys, MaterialDB::MatProperty::Name matPrpty);
	private:
		//TODO : initialize these values when constructing the object. Then we can judge the value when they are used
		//the value of these physical properties is normalized value.
		double bandgap; ///< bandgap of the material
		double electrostaticPotential; ///< potential, normalized
		double conductionBandEnergy; ///< i.e. conduction band edge
		double valenceBandEnergy; ///< i.e. valence band edge
		double electronAffinity; ///< electron affinity
		double electronMass; ///< effective electron mass
		double holeMass; ///< effective hole mass
	};
}