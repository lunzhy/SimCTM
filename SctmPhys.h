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
#ifndef _SCTMPHYS_H_
#define _SCTMPHYS_H_

#include <cmath>
#include <map>
#include <vector>
#include "SctmMath.h"
#include "Material.h"
#include "DomainDetails.h"

//in order to use the pointer of FDVertex in this head file
class FDVertex;
class FDDomain;

typedef std::map<MaterialDB::Mat::Name, double> PrptyMap;

/// @brief This namespace contains all the physics used in the simulation 
///
/// The common physical parameters are defined here. And other classes related to
/// the physics problems are defined here.
namespace SctmPhys
{
	// Physical constant
	const double ElementaryCharge			= 1.602176487e-19;		// in [C]
	const double PlanckConstant				= 6.62606896e-34;		// in [J.s]
	const double BoltzmanConstant			= 1.3806504e-23;			// in [J/K]
	const double ElectronMass				= 9.10938215e-31;		// in [kg]
	const double VacuumDielectricConstant	= 8.854187817e-12;		// in [F/m]
	const double RoomTemperature			= 300;					// in [K]
	const double IntrinsicConcentration		= 1e10;					// in [cm^-3]
	
	// unit conversion
	const double cm_in_m					= 1e-2;
	const double nm_in_cm					= 1e-7;
	const double cm_in_nm = 1e7;
	const double per_sqr_m_in_per_sqr_cm	= 1e-4; // per square meter in per square centimeter 1/m^2 = 1e-4 /cm^2
	const double per_cm3_in_per_m3			= 1e6; // 1 cm^-3 = 1e6 m^-3
	
	// simplified physical constant
	extern const double &m0;
	extern const double &k0;
	extern const double &h;
	extern const double &eps;
	extern const double &T0;
	extern const double &q;
	extern const double &ni;
	const double hbar = h / 2 / SctmMath::PI;

	
	extern double ReferencePotential; // the reference potential used in calculating the conduction/valence band energy, stored in normalized value

	/// @brief SetPhysConstant is used to set the constant parameters which are initialized with the material parameters
	/// 
	/// These physical constant cannot be set to const, because their definitions needs to call the other methods.
	/// 
	/// @pre
	/// @return void
	/// @note
	void SetPhysConstant();

	double CalculateFlatbandShift_domain(FDDomain *domain);
	double CalculateFlatbandShift_slice_for1D(FDVertex *channelVert);
	double CalculateFlatbandShift_slice_for2D(FDVertex *channelVert);

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
			Bandgap, ///< bandgap
			NetFreeCharge, ///< the total net charge belongs to the vertex 
			DensityControlArea, ///< density control area of trapping layer
			ElectricField, ///< the magnitude of electric field
			ElectricField_X, ///< electric field in X direction
			ElectricField_Y, ///< electric field in Y direction
			ElectricFieldTrap_Y, ///< electric field restricted in trapping layer in Y direction 
			ElectricFieldTrap_X, ///< electric field restricted in trapping layer in X direction
			ElectricFieldTrap, ///< electric field restricted in trapping layer
			DielectricConstant,
			Temperature, ///< the lattice temperature

			//for electrons
			eDOSMass, ///< electron effective DOS mass
			eMass, ///< electron effective mass
			eMobility, ///< the electron mobility
			eDensity, ///< the electron density
			eCurrentDensity_X, ///< electron current density in X direction
			eCurrentDensity_Y, ///< electron current density in Y direction
			eCurrentDensity, ///< the magnitude of electron current density
			eEffDOS, ///< effective electron density of states, in [cm^-3]
			eThermalVelocity, ///< electron thermal velocity
			//the physics properties below are properties needed by the solver pack
			//CAUTION: currently, the tunneling coefficient stores tunneling-out coefficient of the boundary vertex.
			eTunnelCoeff, ///< the tunneling coefficient of this vertex, in [A*cm]
			eCurrDensMFN_X, ///< the x-direction value of MFN tunneling current density of inner vertex 
			eCurrDensMFN_Y, ///< the y-direction value of MFN tunneling current density of inner vertex
			eSubsCurrDensB2T, ///< electron current density from substrate in calculation of band-to-trap tunneling 

			//for holes
			hMass,
			hDOSMass,
			hMobility,
			hDensity,
			hTunnelCoeff,
			hCurrentDensity_X,
			hCurrentDensity_Y,
			hCurrentDensity,
			hThermalVelocity,
			hEffDOS,
		};

		/// @brief PhysProperty is the construction method for this class
		/// 
		/// All the physical values are initialized to zero in the construction method.
		/// 
		/// @pre
		/// @return 
		/// @note
		PhysProperty(FDVertex *_vert);
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
		double GetPhysPrpty(Name prptyName, MaterialDB::Mat::Name matName = MaterialDB::Mat::ErrorMaterial) const;
		/// @brief FillVertexPhysUsingMatPropty is used to set vertex-based physical value using material-based property
		/// 
		/// Vertex-based physical value is set in consideration of four(or less) adjacent elements of the specific vertex.
		/// value = (weighted sum of the material property with respect of element area) / total area
		/// 
		/// @param FDVertex * vertex the vertex to set
		/// @param PhysProperty::Name vertexPhys
		/// @param MaterialDB::MatProperty::Name matPrpty
		/// @param FDRegion::RegionType rType the designated type of region
		/// @pre
		/// @return void
		/// @note This method is not checked until now. || Used.
		void FillVertexPhysUsingMatPrpty(PhysProperty::Name vertexPhys,
			MaterialDB::MatProperty::Name matPrpty, bool onlyTrapRegion = false);
		void CalculateDensityControlArea();
		void UpdateValue(Name prptyName, double val);
		void SetMultiPrpty(PhysProperty::Name vertPhy, MaterialDB::MatProperty::Name matPrpty);
		bool HasMultiPrpty(PhysProperty::Name prptyName) const;
		std::vector<MaterialDB::Mat::Name> const& GetRelatedMaterialNames();
	private:
		FDVertex *vertSelf; ///< the vertex this physics property belongs to
		//TODO : initialize these values when constructing the object. Then we can judge the value when they are used.
		//the value of these physical properties is normalized value.
		double bandgap; ///< bandgap of the material, normalized, in[eV]
		double electrostaticPotential; ///< potential, normalized
		double electronAffinity; ///< electron affinity
		double densControlArea; ///< density control area, only valid in trapping layer. Only sum up the area in adjacent trapping layers.
		double epsilon; ///< dielectric constant
		double temperature; ///< the lattice temperature of this vertex

		//for electrons
		double e_DOSmass; ///< effective DOS mass
		double e_mass; ///< effective electron mass
		double e_mobility; ///< electron mobility
		double e_density; ///< the electron density
		//the physics properties below are properties needed by the solver pack
		double e_tunnelCoeff; ///< the tunneling coefficient of this vertex, in [A*cm]
		double e_currdensMFN_X; ///< the x-direction value of MFN tunneling current density of inner vertex
		double e_currdensMFN_Y; ///< the y-direction value of MFN tunneling current density of inner vertex
		double e_subsCurrDensB2T; ///< band-to-trap electron current density from substrate in the calculation of band-to-trap tunneling 

		//for holes
		double h_mass;
		double h_DOSmass;
		double h_mobility;
		double h_density;
		double h_tunnelCoeff;

		//the maps below is used to store the properties of vertex that belongs to different materials
		PrptyMap multiElectronAffinity;
		PrptyMap multiBandgap;
		PrptyMap multiDielectricConstant; //not used
		std::vector<MaterialDB::Mat::Name> relatedMatName;

		double getMultiPrptyValue(PhysProperty::Name vertPhy, MaterialDB::Mat::Name matName) const;

	};

	class TrapProperty
	{
	public:
		enum Name
		{
			EpsilonTrapping, ///< dielectric constant for trapping layer
			NetTrappedCharge,
			TrapDensity,

			//for electrons
			eTrapped,
			eTrapDensity,
			eCrossSection,
			eTrapCrossSection,
			eEnergyFromCondBand,
			eOccupation,
			eEmptyTrapDens,
			eCaptureCoeff_J_Model,
			eCaptureCoeff_V_Model,
			eTrappedCapCoeff_V_Model,
			eTrappedCapCoeff_J_Model,
			eEmissionCoeff_BasicSRH, ///< the electron emission rate of basic SRH process
			eCoeff_B2T, ///< electron coefficient in band-to-trap tunneling, in [1/s]
			eFrequency_T2B, ///< electron trap-to-band tunneling out frequency
			eEmissionCoeff_T2B, ///< electron trap-to-band tunneling out coefficient
			eTransCoeffT2B, ///< electron trap-to-band tunneling out transmission coefficient
			eEmissionCoeff_PF, ///< electron emission coefficient of Poole-Frenkel  effect
			//below are parameters to store
			eFrequency_PF, ///< electron emission frequency of Poole-Frenkel effect
			eTrapEnergyDecreasePF, ///< electron trap energy decrease due to Poole-Frenkel effect

			//for holes
			hCrossSection,
			hTrapCrossSection,
			hFrequency_T2B,
			hFrequency_PF,
			hEnergyFromValeBand,
			hTrapped,
			hTrapDensity,
			hOccupation,
			hEmptyTrapDens,
			hCaptureCoeff_J_Model,
			hCaptureCoeff_V_Model,
			hTrappedCapCoeff_V_Model,
			hTrappedCapCoeff_J_Model,
			hEmissionCoeff_BasicSRH,
			hTrapEnergyDecreasePF,
			hEmissionCoeff_PF,
		};
		TrapProperty(FDVertex *_vert);
		/// @brief GetTrapPrpty is used to get the specific trap property value.
		/// 
		///
		/// 
		/// @param TrapProperty::Name trapPrpty
		/// @pre
		/// @return double, the return value is in normalized value.
		/// @note
		double GetTrapPrpty(TrapProperty::Name trapPrpty) const;
		void SetTrapPrpty(TrapProperty::Name trapPrpty, double val);

		void FillTrapPrptyUsingMatPrpty(TrapProperty::Name trapPrpty, MaterialDB::MatProperty::Name matPrpty);

	private:
		FDVertex *vertSelf;
		double epsTrapping;
		double trapDensity; ///< total trap density both for electrons and holes

		double e_trapped;
		double e_trapDensity;
		double e_crossSection;
		double e_trapCrossSection;
		double e_energyFromCondBand; ///< electron trap energy
		double e_frequencyT2B; ///< electron Trap-to-Band tunneling-out frequency
		double e_transCoeffT2B; ///< electron transmission coefficient in trap-to-band tunneling
		double e_frequencyPF; ///< electron emission frequency of Poole - Frenkel effect

		double h_trapped;
		double h_trapDensity;
		double h_crossSection;
		double h_trapCrossSection;
		double h_frequencyT2B;
		double h_frequencyPF;
		double h_energyFromValeBand; ///< hole trap energy
	};
}

#endif