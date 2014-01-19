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
	const double per_sqr_m_in_per_sqr_cm	= 1e-4; // per square meter in per square centimeter 1/m^2 = 1e-4 /cm^2
	
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

	double CalculateFlatbandShift(FDDomain *domain);

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
			eMass, ///< electron effective mass
			hMass, ///< hole effective mass
			Bandgap, ///< bandgap
			NetCharge, ///< the total net charge belongs to the vertex 
			eMobility, ///< the electron mobility
			eDensity, ///< the electron density
			DensityControlArea, ///< density control area of trapping layer
			ElectricField, ///< the magnitude of electric field
			ElectricField_X, ///< electric field in X direction
			ElectricField_Y, ///< electric field in Y direction
			eCurrentDensity_X, ///< electron current density in X direction
			eCurrentDensity_Y, ///< electron current density in Y direction
			eCurrentDensity, ///< the magnitude of electron current density
			DielectricConstant,
			eEffDOS, ///< effective electron density of states, in [cm^-3]
			eThermalVelocity, ///< electron thermal velocity
			Temperature, ///< the lattice temperature

			//the physics properties below are properties needed by the solver pack
			//CAUTION: currently, the tunneling coefficient stores tunneling-out coefficient of the boundary vertex.
			TunnelCoeff, ///< the tunneling coefficient of this vertex, in [A*cm]
			eCurrDensMFN_X, ///< the x-direction value of MFN tunneling current density of inner vertex 
			eCurrDensMFN_Y, ///< the y-direction value of MFN tunneling current density of inner vertex
			eSubsCurrDensB2T, ///< electron current density from substrate in calculation of band-to-trap tunneling 
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
			MaterialDB::MatProperty::Name matPrpty);
		void FillVertexPhysUsingMatPrpty(PhysProperty::Name vertexPhys,
			MaterialDB::MatProperty::Name matPrpty, FDRegion::TypeName rType);
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
		//double conductionBandEnergy; ///< i.e. conduction band edge
		//double valenceBandEnergy; ///< i.e. valence band edge
		double electronAffinity; ///< electron affinity
		double e_mass; ///< effective electron mass
		double h_mass; ///< effective hole mass
		//double netCharge; ///< total net charge belongs to the vertex
		double e_mobility; ///< electron mobility
		double e_density; ///< the electron density
		double controlArea; ///< density control area, only valid in trapping layer. Only sum up the area in adjacent trapping layers.
		double epsilon; ///< dielectric constant
		double temperature; ///< the lattice temperature of this vertex

		//the physics properties below are properties needed by the solver pack
		double tunnelCoeff; ///< the tunneling coefficient of this vertex, in [A*cm]
		double e_currdensMFN_X; ///< the x-direction value of MFN tunneling current density of inner vertex
		double e_currdensMFN_Y; ///< the y-direction value of MFN tunneling current density of inner vertex
		double e_subsCurrDensB2T; ///< band-to-trap electron current density from substrate in the calculation of band-to-trap tunneling 

		//the maps below is used to store the properties of vertex that belongs to different materials
		PrptyMap multiElectronAffinity;
		PrptyMap multiBandgap;
		std::vector<MaterialDB::Mat::Name> relatedMatName;

		double getMultiPrptyValue(PhysProperty::Name vertPhy, MaterialDB::Mat::Name matName) const;

	};

	class TrapProperty
	{
	public:
		enum Name
		{
			eTrapped,
			eTrapDensity,
			eCrossSection,
			EnergyFromCondBand,
			eOccupation,
			eEmptyTrapDens,
			NetCharge,
			eCaptureCoeff_J_Model,
			eCaptureCoeff_V_Model,
			eEmissionCoeff_BasicSRH, ///< the electron emission rate of basic SRH process
			eCoeff_B2T, ///< electron coefficient in band-to-trap tunneling, in [1/s]
			eFrequencyT2B, ///< electron trap-to-band tunneling out frequency
			eEmissionCoeff_T2B, ///< electron trap-to-band tunneling out coefficient
			eTransCoeffT2B, ///< electron trap-to-band tunneling out transmission coefficient

		};
		TrapProperty(FDVertex *_vert);
		double GetTrapPrpty(TrapProperty::Name trapPrpty) const;
		void SetTrapPrpty(TrapProperty::Name trapPrpty, double val);

		void FillTrapPrptyUsingMatPrpty(TrapProperty::Name trapPrpty, MaterialDB::MatProperty::Name matPrpty);

	private:
		FDVertex *vertSelf;
		double e_trapped;
		double e_trapDensity;
		double e_crossSection;
		double energyFromCondBand;
		double e_frequencyT2B; ///< electron Trap-to-Band tunneling-out frequency
		double e_transCoeffT2B; ///< electron transmission coefficient in trap-to-band tunneling
	};
}

#endif