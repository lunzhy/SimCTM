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

#include <string>
using std::string;

/// @brief This namespace contains the material parameters.
///
/// This is used as a material database.
namespace MaterialDB
{
	/// @brief This class is a data structure to store the parameters for different materials used in the simulation.
	class Material
	{
	private:
		string		name; ///< material name
		double		dielectricConstant; ///< dielectric constant
		double		bandgap; ///< bandgap, in [eV]
		double		electronAffinity; ///< electron affinity energy, in [eV]
		double		electronMass; ///< electron effective mass, in [m0]
		double		holeMass; ///< hole effective mass, in [m0]
		double		electronDOS;
		double		holeDOS;
	public:
		/// @brief Material is the construction method of this class
		///  
		/// The object of this class is constructed with given material name.
		/// 
		/// @param string _name
		/// @pre
		/// @return 
		/// @note
		Material(string _name):name(_name) {};

		//The methods below are used to encapsulate the private members of this class.
		double		DielectricConstant() const			{ return dielectricConstant;	}
		void		DielectricConstant(double val)		{ dielectricConstant = val;		}
		double		Bandgap() const						{ return bandgap;				}
		void		Bandgap(double val)					{ bandgap = val;				}
		double		ElectronAffinity() const			{ return electronAffinity;		}
		void		ElectronAffinity(double val)		{ electronAffinity = val;		}
		double		ElectronMass() const				{ return electronMass;			}
		void		ElectronMass(double val)			{ electronMass = val;			}
		double		HoleMass() const					{ return holeMass;				}
		void		HoleMass(double val)				{ holeMass = val;				}
		double		ElectronDOS() const					{ return electronDOS;			}
		void		ElectronDOS(double val)				{ electronDOS = val;			}
		double		HoleDOS() const						{ return holeDOS;				}
		void		HoleDOS(double val)					{ holeDOS = val;				}
	};


	/// @brief This class just provides a container to restrict the enum of material property name.
	class MatProperty
	{
	public:
		enum Name
		{
			Mat_DielectricConstant, ///< dielectric constant of the material
			Mat_Bandgap, ///< bandgap of the material
			Mat_ElectronAffinity, ///< electron affinity of the material
			Mat_ElectronMass, ///< electron effective mass of the material
			Mat_HoleMass, ///< hole effective mass of the material
			Mat_ElectronDOS, ///< electron density of states of the material
			Mat_HoleDOS ///< hoe density of states of the material
		};
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
	/// @brief SetMaterials is used to set the values of material properties.
	/// 
	///
	/// 
	/// @pre
	/// @return void
	/// @note
	void SetMaterials();

	//TODO: conceive the method that accounts for material specification
	extern Material Silicon;
	extern Material SiO2;
	extern Material Si3N4;
}

#endif