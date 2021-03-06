/**
* @file Material.h
* @brief this file is used to define the structure of class related to materials used in the simulation
*
*
*
* @author
* @version 
* @date 2013-7-2   10:14
* @note
* @todo
*/
#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <map>
#include <string>

/// @brief This namespace contains the material parameters.
///
/// This is used as a material database.
namespace MaterialDB
{
	class Mat
	{
	public:
		enum Name
		{
			ErrorMaterial,
			Silicon,
			SiO2,
			Si3N4,
			HfO2,
			Al2O3,
		};
		static Mat::Name Parse(const std::string &matStr);
	};
	/// @brief This class just provides a container to restrict the enum of material property name.
	class MatProperty
	{
	public:
		enum Name
		{
			Mat_DielectricConstant, ///< dielectric constant of the material
			Mat_HighFrqDielConst, ///< high-frequency dielectric constant
			Mat_Bandgap, ///< bandgap of the material
			Mat_ElectronAffinity, ///< electron affinity of the material
			Mat_ElectronMass, ///< electron effective mass of the material
			Mat_HoleMass, ///< hole effective mass of the material
			Mat_ElecDOSMass, ///< electron effective DOS mass
			Mat_HoleDOSMass, ///< hole effective DOS mass
			Mat_ElectronDiffusion, ///< electron diffusion coefficient of the material
			Mat_HoleDiffusion, ///< hole diffusion coefficient of the material
			Mat_ElectronMobility, ///< electron mobility of the material
			Mat_HoleMobility, ///< hole mobility of the material
			Mat_ElecXSection, ///< electron trap cross section
			Mat_HoleXSection, ///< hole trap cross section
			Mat_ElecTrappedXSection, ///< cross section for trapped electrons
			Mat_HoleTrappedXSection, ///< cross section for trapped holes
			Mat_ElecTrapEnergyFromCB, ///< electron trap energy from conduction band
			Mat_HoleTrapEnergyFromVB, ///< hole trap energy from valence band
			Mat_ElecFrequencyT2B, ///< electron Trap-to-Band tunneling out frequency
			Mat_HoleFrequencyT2B, ///< hole Trap-to-Band tunneling out frequency
			Mat_ElecFrequencyPF, ///< electron emission frequency in Poole-Frenkel effect
			Mat_HoleFrequencyPF, ///< hole emission frequency in Poole-Frenkel effect
		};
	};

	/// @brief This class is a data structure to store the parameters for different materials used in the simulation.
	///
	/// The data is stored in normalized value
	class Material
	{
	private:
		double temperature;
		Mat::Name name; ///< material name
		double dielectricConstant; ///< dielectric constant
		double highFrqDC; /// high-frequency dielectric constant
		double bandgap; ///< bandgap, in [eV]
		double electronAffinity; ///< electron affinity energy, in [eV]
		double elecMass; ///< electron effective tunneling mass, in [m0]
		double holeMass; ///< hole effective tunneling mass, in [m0]
		double elecDOSMass; ///< electron effective DOS mass, in [m0]
		double holeDOSMass; ///< hole effect DOS mass, in [m0]
		double electronDiffusion; ///< electron diffusion coefficient, in [D0]
		double holeDiffusion; ///< hole diffusion coefficient, in [D0]
		double electronMobility; ///< electron mobility, in [cm^2/V/s]
		double holeMobility; ///< hole mobility, in [cm^2/V/s]
		double elecXSection; ///< electron trap cross section, in [cm^2]
		double holeXSection; ///< hole trap cross section, in [cm^2]
		double elecTrappedXSection; ///< cross section for trapped electrons
		double holeTrappedXSection; ///< cross section for trapped holes
		double elecTrapEnergyFromCB; ///< electron trap energy from conduction band, in eV
		double holeTrapEnergyFromVB; ///< hole trap energy from valence band
		double elecFrequencyT2B; ///< electron Trap-to-Band tunneling out frequency
		double holeFrequencyT2B; ///< hole Trap-to-Band tunneling out frequency
		double elecFrequencyPF; ///< electron emission frequency in Poole-Frenkel effect
		double holeFrequencyPF; ///< hole emission frequency in Poole-Frenkel effect
		
	public:
		/// @brief Material is the construction method of this class
		///  
		/// The object of this class is constructed with given material name.
		/// 
		/// @param string _name
		/// @pre
		/// @return 
		/// @note
		Material(Mat::Name _name);
		Mat::Name MatName() { return name; }

		//The methods below are used to encapsulate the private members of this class.
		double		DielectricConstant() const;
		void		DielectricConstant(double val);
		double		HighFrqDielConst() const;
		void		HighFrqDielConst(double val);
		double		Bandgap() const;
		void		Bandgap(double val);
		double		ElectronAffinity() const;
		void		ElectronAffinity(double val);
		
		double		ElecDOSMass() const;
		void		ElecDOSMass(double val);
		double		HoleDOSMass() const;
		void		HoleDOSMass(double val);
		
		double		ElectronMass() const;
		void		ElectronMass(double val);
		double		HoleMass() const;
		void		HoleMass(double val);
		
		double		ElectronDiffusion() const;
		void		ElectronDiffusion(double val);
		double		HoleDiffusion() const;
		void		HoleDiffusion(double val);
		
		double		ElectronMobility() const;
		void		ElectronMobility(double val);
		double		HoleMobility() const;
		void		HoleMobility(double val);
		
		double		ElecXSection() const;
		void		ElecXSection(double val);
		double		HoleXSection() const;
		void		HoleXSection(double val);

		double		ElecTrappedXSection() const;
		void		ElecTrappedXSection(double val);
		double		HoleTrappedXSection() const;
		void		HoleTrappedXSection(double val);
		
		double		ElecTrapEnergyFromCB() const;
		void		ElecTrapEnergyFromCB(double val);
		double		HoleTrapEnergyFromVB() const;
		void		HoleTrapEnergyFromVB(double val);

		double		ElecFrequencyT2B() const;
		void		ElecFrequencyT2B(double val);
		double		HoleFrequencyT2B() const;
		void		HoleFrequencyT2B(double val);

		double		ElecFrequencyPF() const;
		void		ElecFrequencyPF(double val);
		double		HoleFrequencyPF() const;
		void		HoleFrequencyPF(double val);
		
	};

	/// @brief GetMatPrpty is called to get the value of material property with given material
	/// 
	/// This method is used with the pointer to the material that is examined. And the property name is given
	/// in the name within the enum of material property name
	/// 
	/// @param Material * theMaterial
	/// @param MatProperty::Name prptyName
	/// @pre
	/// @return double
	/// @note
	double GetMatPrpty(Material *theMaterial, MatProperty::Name prptyName);
	/// @brief SetMatPrpty is called to set the value of material property with given material
	/// 
	/// This method is barely used. Because the material properties are set within this namespace using other methods.
	/// The initialization of material is not used after the preparation of the simulation.
	/// 
	/// @param Material * theMaterial
	/// @param MatProperty::Name prptyName
	/// @pre
	/// @return void
	/// @note
	void SetMatPrpty(Material *theMaterial, MatProperty::Name prptyName);
	/// @brief SetMaterials_Directly is used to set the values of material properties directly for debugging.
	/// 
	///
	/// 
	/// @pre
	/// @return void
	/// @note
	void SetMaterials_Directly();

	void SetMaterial_FromParFile();

	//TODO: consider the method that accounts for material specification
	//extern std::map<Mat::Name, Material*> MaterialMap;
	Material* GetMaterial(Mat::Name matname);
}

#endif