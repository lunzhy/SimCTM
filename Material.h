/**
* @file Material.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-2   10:14
* @note
* @todo
*/
#pragma once
#include <string>


/// @brief This namespace contains the material parameters.
///
/// This is used as a material database.
using std::string;
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
		double		electronMass;
		double		holeMass;
		double		electronDOS;
		double		holeDOS;
	public:
		Material(string _name):name(_name) {};

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

	//TODO: conceive the method that accounts for material specification
	extern Material Silicon;
	extern Material SiO2;
	extern Material Si3N4;
	void SetMaterials();
}