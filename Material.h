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
namespace MaterialDB
{
	struct Material
	{
		std::string name;
		double epsilon;
	};

	Material Silicon = {"Silicon", 
						11.9};
}