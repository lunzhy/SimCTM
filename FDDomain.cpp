/**
* @file FDDomain.cpp
* @brief This file contains the definition/implementation of the member methods.
*
* 
*
* @author
* @version 
* @date 2013-7-4   19:40
* @note
* @todo
*/

#include "FDDomain.h"
#include "Material.h"

FDElement * FDDomain::getElement(unsigned int id)
{
	//need to judge if the id is appropriate
	return elements.at(id);
	//return elements[id];
}

FDVertex * FDDomain::getVertex(unsigned int id)
{
	return vertices.at(id);
	//return vertices[id];
}

FDRegion * FDDomain::getRegion(unsigned int id)
{
	return regions.at(id);
	//return regions[id];
}

FDContact * FDDomain::getContact(unsigned int id)
{
	return contacts.at(id);
}
