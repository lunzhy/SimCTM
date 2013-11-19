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
#include "DomainDetails.h"
#include "SctmUtils.h"
#include "SctmPhys.h"
using SctmPhys::PhysProperty;
using namespace SctmUtils;

FDElement * FDDomain::GetElement(unsigned int id)
{
	//need to judge if the id is appropriate
	return elements.at(id);
	//return elements[id];
}

FDVertex * FDDomain::GetVertex(unsigned int id)
{
	return vertices.at(id);
	//return vertices[id];
}

FDRegion * FDDomain::GetRegion(unsigned int id)
{
	return regions.at(id);
	//return regions[id];
}

FDContact * FDDomain::GetContact(unsigned int id)
{
	return contacts.at(id);
}

std::vector<FDVertex *> & FDDomain::GetVertices()
{
	return this->vertices;
}

bool FDDomain::isValidElem(FDElement *elem)
{
	return elem != NULL;
}

bool FDDomain::isNotTrappingElem(FDElement *elem)
{
	if (elem == NULL)
		return true;
	else
	{
		return elem->Region->Type != FDRegion::Trapping;
	}
}

void FDDomain::setBoundaryCondition()
{
	FDVertex *currVertex;
	for (std::size_t iVer = 0; iVer != vertices.size(); ++iVer)
	{
		currVertex = GetVertex(iVer);
		setVertBC_Potential(currVertex);
		setVertBC_eDensity(currVertex);
	}
}

void FDDomain::setVertBC_Potential(FDVertex *vert)
{
	bool isValid_NW = isValidElem(vert->NorthwestElem);
	bool isValid_NE = isValidElem(vert->NortheastElem);
	bool isValid_SE = isValidElem(vert->SoutheastElem);
	bool isValid_SW = isValidElem(vert->SouthwestElem);

	//firstly, to decide if the vertex is at a contact.
	if ( vert->IsAtContact() )
	{
		//the gate name is in accordance with the name specified in setting domain details
		if (vert->Contact->ContactName == "Gate")
		{
			//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
			vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Dirichlet, 0);
			return;
		}
		else if (vert->Contact->ContactName == "Channel")
		{
			//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
			vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Dirichlet, 0);
			return;
		}
	}

	//for vertex that is not correlated to a contact. The value of 0 indicated an artificial boundary there.
	//Northwest corner
	if ( !isValid_NW && !isValid_NE &&
		!isValid_SW && isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-vert->EastLength, vert->SouthLength));
		return;
	}

	//Northeast corner
	if ( !isValid_NW && !isValid_NE && 
		isValid_SW && !isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(vert->WestLength, vert->SouthLength));
		return;
	}

	//Southeast corner
	if ( isValid_NW && !isValid_NE && 
		!isValid_SW && !isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(vert->WestLength, -vert->NorthLength));
		return;
	}

	//Southwest corner
	if ( !isValid_NW && isValid_NE && 
		!isValid_SW && !isValid_SE)
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-vert->EastLength, -vert->NorthLength));
		return;
	}

	//North side
	if ( !isValid_NW && !isValid_NE && 
		isValid_SW && isValid_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(0, 1));
		return;
	}

	//East side
	if ( isValid_NW && !isValid_NE && 
		isValid_SW && !isValid_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(1, 0));
		return;
	}

	//South side
	if ( isValid_NW && isValid_NE && 
		!isValid_SW && !isValid_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(0, -1));
		return;
	}

	//West side
	if ( !isValid_NW && isValid_NE && 
		!isValid_SW && isValid_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-1, 0));
		return;
	}
}

void FDDomain::setVertBC_eDensity(FDVertex *vert)
{
	bool notTrapping_NW = isNotTrappingElem(vert->NorthwestElem);
	bool notTrapping_NE = isNotTrappingElem(vert->NortheastElem);
	bool notTrapping_SE = isNotTrappingElem(vert->SoutheastElem);
	bool notTrapping_SW = isNotTrappingElem(vert->SouthwestElem);

	bool valid_NW = isValidElem(vert->NorthwestElem);
	bool valid_NE = isValidElem(vert->NortheastElem);
	bool valid_SE = isValidElem(vert->SoutheastElem);
	bool valid_SW = isValidElem(vert->SouthwestElem);

	//When dealing with the normal direction of the boundary condition in terms eDensity problem, actually, the vector value
	//of the boundary condition is not the real normal vector value of the trapping region. For, example, to a corner vertex
	//of the trapping region, the normal vector is determined considering if its adjacent element is valid.

	//Northwest corner
	if ( notTrapping_NW && notTrapping_NE &&
		notTrapping_SW && !notTrapping_SE )
	{
		if (              valid_NE &&
			!valid_SW )
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, 1));
			return;
		}
		if (              !valid_NE &&
			valid_SW )
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-vert->EastLength, vert->SouthLength));
		return;
	}

	//Northeast corner
	if ( notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && notTrapping_SE )
	{
		if ( valid_NW &&
			!valid_SE)
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0 ,1));
			return;
		}
		if ( !valid_NW &&
			valid_SE)
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(vert->WestLength, vert->SouthLength));
		return;
	}

	//Southeast corner
	if ( !notTrapping_NW && notTrapping_NE && 
		notTrapping_SW && notTrapping_SE )
	{
		if (			!valid_NE &&
			valid_SW)
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
			return;
		}
		if (			valid_NE &&
			!valid_SW)
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(vert->WestLength, -vert->NorthLength));
		return;
	}

	//Southwest corner
	if ( notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && notTrapping_SE)
	{
		if ( valid_NW &&
			!valid_SE )
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
			return;
		}
		if ( !valid_NW &&
			valid_SE )
		{
			vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-vert->EastLength, -vert->NorthLength));
		return;
	}

	//North side
	if ( notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && !notTrapping_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, 1));
		return;
	}

	//East side
	if ( !notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && notTrapping_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
		return;
	}

	//South side
	if ( !notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && notTrapping_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
		return;
	}

	//West side
	if ( notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && !notTrapping_SE)
	{
		vert->BndCond.SetBndCond(FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
		return;
	}

}

void FDDomain::BuildDomain()
{
	UtilsTimer.Set();
	//Initialize the vectors in FDDomain
	vertices.clear();
	vertsTrapping.clear();
	elements.clear();
	regions.clear();
	contacts.clear();

	//build the data and mesh structure of simulated region, this is a pure virtual method
	buildStructure();
	//set the physics value related to vertex
	setVertexPhysics();
	//set the boundary condition, the specific value is not considered in this class.
	setBoundaryCondition();

	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());
}

void FDDomain::setVertexPhysics()
{
	FDVertex * currVertex = NULL;
	FDElement * currElem = NULL;
	double tot = 0; //total area
	double sum = 0; //sum corresponds to integral value
	double physValue = 0;

	using namespace MaterialDB;
	using namespace SctmPhys;
	using std::vector;

	vector<MatProperty::Name> matPrptys;
	vector<PhysProperty::Name> verPrptys; //vertex-based physical property
	matPrptys.push_back(MatProperty::Mat_ElectronAffinity); verPrptys.push_back(PhysProperty::ElectronAffinity);
	matPrptys.push_back(MatProperty::Mat_ElectronMass); verPrptys.push_back(PhysProperty::eMass);
	matPrptys.push_back(MatProperty::Mat_Bandgap); verPrptys.push_back(PhysProperty::Bandgap);
	matPrptys.push_back(MatProperty::Mat_ElectronMobility); verPrptys.push_back(PhysProperty::eMobility);

	//iteration over the vertices
	for (std::size_t iVer = 0; iVer != this->vertices.size(); ++iVer)
	{
		currVertex = GetVertex(iVer);
		//iteration over the physical properties to be set from material property
		for (std::size_t iPrpty = 0; iPrpty != matPrptys.size(); ++iPrpty)
		{
			//filling vertex physics using material property
			//The method for filling vertex-based physical value using material-based value is ready
			//electron mobility is only valid in the trapping region
			if (matPrptys.at(iPrpty) == MatProperty::Mat_ElectronMobility)
			{
				currVertex->Phys->FillVertexPhysUsingMatPropty(currVertex, verPrptys.at(iPrpty), matPrptys.at(iPrpty), FDRegion::Trapping);
			}
			else
			{
				currVertex->Phys->FillVertexPhysUsingMatPropty(currVertex, verPrptys.at(iPrpty), matPrptys.at(iPrpty));
			}
			currVertex->Phys->CalculateDensityControlArea(currVertex);
			/*
			tot = 0; sum = 0;
			currElem = currVertex->SouthwestElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->SoutheastElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->NortheastElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->NorthwestElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->Mat, matPrptys.at(iPrpty)) * currElem->Area : 0;

			physValue = sum / tot;
			currVertex->Phys->SetPhysPrpty(verPrptys.at(iPrpty), physValue);
			*/
		}
	}
}
