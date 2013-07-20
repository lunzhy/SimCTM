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

#include <string>


/// @brief This namespace contains the material parameters.
///
/// This is used as a material database.
namespace MaterialDB
{
	/// @brief This struct is a data structure to store the parameters for different materials used in the simulation.
	struct Material
	{
		std::string name; ///< material name
		double epsilon; ///dielectric constant
	};

	Material Silicon = {"Silicon", 
						11.9};
}