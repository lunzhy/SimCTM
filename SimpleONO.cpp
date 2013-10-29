/**
* @file DomainTest.cpp
* @brief
*
* For rectangle grids
*
* @author
* @version 
* @date 2013-7-3   10:00
* @note
* @todo
*/

#include "SimpleONO.h"
#include "Material.h"
#include "SctmUtils.h"
#include <iostream>
#include "SctmPhys.h"
using namespace SctmUtils;

void SimpleONO::BuildDomain()
{
	//Initialize the vectors
	vertices.clear();
	elements.clear();
	regions.clear();

	//set the data structure of simulated region
	prepareStructure(); //mainly set the structure and mesh parameters
	setDomainDetails(); //fill in the regions, vertices and elements
	setAdjacency(); //set the adjacency of vertices and elements

	//fill in the physical parameters and values
	setVertexPhysics();
	setBoundaryCondition(true);
	//stuffPotential();
	//refreshBandEnergy();
	//printStructure();
}

void SimpleONO::prepareStructure()
{
	setParameters();
}

void SimpleONO::setParameters()
{
	double nm_in_cm = SctmPhys::nm_in_cm;

	////////////////////////////////////////////////////////////////////////////
	//modify here to change the structures
	//current parameters is set according to the structure prepared in Sentaurus
	double xLength_in_nm = 10;
	double yLengthTunnel_in_nm = 4;
	double yLengthTrap_in_nm = 10;
	double yLengthBlock_in_nm = 9;
	int xGridNumber = 50; //the grid number, not vertex number
	int yGridNumberTunnel = 5;
	int yGridNumberTrap = 50;
	int yGridNumberBlock = 5;
	////////////////////////////////////////////////////////////////////////////
	//set geometric class members
	//here, the length of the parameter is conversed to [cm]
	xLength = xLength_in_nm * nm_in_cm;
	xCntVertex = xGridNumber + 1;
	yLengthTunnel = yLengthTunnel_in_nm * nm_in_cm;
	yCntVertexTunnel = yGridNumberTunnel + 1;
	yLengthTrap = yLengthTrap_in_nm * nm_in_cm;
	yCntVertexTrap = yGridNumberTrap + 1;
	yLengthBlock = yLengthBlock_in_nm * nm_in_cm;
	yCntVertexBlock = yGridNumberBlock + 1;
	yCntTotalVertex = yCntVertexTunnel + yCntVertexTrap + yCntVertexBlock - 1 - 1;

	xGrid = xLength / ( xCntVertex - 1 );
	yGridTunnel = yLengthTunnel / ( yCntVertexTunnel - 1 );
	yGridTrap = yLengthTrap / ( yCntVertexTrap - 1 );
	yGridBlock = yLengthBlock / ( yCntVertexBlock - 1 );
}

void SimpleONO::printStructure()
{
	for (std::vector<FDVertex *>::size_type ix = 0; ix != this->vertices.size(); ++ix)
	{
		std::cout << "id=" << vertices.at(ix)->GetID() 
			<< '\t' << (vertices.at(ix)->WestVertex == NULL ? -1 : vertices.at(ix)->WestVertex->GetID())
			<< '\t' << (vertices.at(ix)->SouthVertex == NULL ? -1 : vertices.at(ix)->SouthVertex->GetID())
			<< '\t' << (vertices.at(ix)->EastVertex == NULL ? -1 : vertices.at(ix)->EastVertex->GetID())
			<< '\t' << (vertices.at(ix)->NorthVertex == NULL ? -1 : vertices.at(ix)->NorthVertex->GetID())
			<< '\t' << (vertices.at(ix)->SouthwestElem == NULL ? -1 : vertices.at(ix)->SouthwestElem->GetInternalID())
			<< '\t' << (vertices.at(ix)->SoutheastElem == NULL ? -1 : vertices.at(ix)->SoutheastElem->GetInternalID())
			<< '\t' << (vertices.at(ix)->NortheastElem == NULL ? -1 : vertices.at(ix)->NortheastElem->GetInternalID())
			<< '\t' << (vertices.at(ix)->NorthwestElem == NULL ? -1 : vertices.at(ix)->NorthwestElem->GetInternalID())
			<< '\t' << vertices.at(ix)->WestLength
			<< '\t' << vertices.at(ix)->SouthLength
			<< '\t' << vertices.at(ix)->EastLength
			<< '\t' << vertices.at(ix)->NorthLength
			<< std::endl;
	}

	for (std::vector<FDElement *>::size_type ix = 0; ix != this->elements.size(); ++ix)
	{
		std::cout << "id=" << elements.at(ix)->GetInternalID() << '\t' << elements.at(ix)->Region->Type
			<< '\t' << elements.at(ix)->SouthwestVertex->GetID()
			<< '\t' << elements.at(ix)->SoutheastVertex->GetID()
			<< '\t' << elements.at(ix)->NortheastVertex->GetID()
			<< '\t' << elements.at(ix)->NorthwestVertex->GetID()
			<< '\t' << elements.at(ix)->WestLength
			<< '\t' << elements.at(ix)->SouthLength
			<< '\t' << elements.at(ix)->EastLength
			<< '\t' << elements.at(ix)->NorthLength
			<< std::endl;
	}
}

void SimpleONO::setDomainDetails()
{
	int cntVertex = 0;
	int cntElement = 0;
	int cntRegion = 0;
	int cntContact = 0;

	double currCoordX = 0.0;
	double currCoordY = 0.0;

	double nm_in_cm = SctmPhys::nm_in_cm;

	////////////////////////////////////////////////////////////////////
	//set vertices
	////////////////////////////////////////////////////////////////////
	Normalization theNorm = Normalization();

	double normCoordX = 0.0;
	double normCoordY = 0.0;
	//firstly, scan the x direction and secondly increment the y coordinate
	//this consequence is in accordance with the FDDomainHelper
	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		//currCoordY = iy * yGridTunnel;
		//currCoordY = theNorm.PushLength(currCoordY);
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			//currCoordX = ix * xGrid;
			//the coordinates are normalized before pushing into vertex
			normCoordX = theNorm.PushLength(currCoordX); //the value of currCoordX is already in cm
			normCoordY = theNorm.PushLength(currCoordY); //the value of currCoordY is already in cm
			vertices.push_back(new FDVertex(cntVertex, normCoordX, normCoordY));
			cntVertex++;
			currCoordX += xNextGridLength(ix);//non-normalized value, in cm
		}
		currCoordX = 0;
		currCoordY += yNextGridLength(iy);//non-normalized value, in cm
	}

	////////////////////////////////////////////////////////////////////
	//set contacts
	////////////////////////////////////////////////////////////////////
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntTotalVertex);
	//the third parameter of FDContact construction is of no use, because the voltage is set in the following process.
	contacts.push_back(new FDContact(cntContact, "Gate", 0));
	cntContact++;
	contacts.push_back(new FDContact(cntContact, "Channel", 0));//here channel is an imagined contact
	FDContact *currContact = NULL;
	FDVertex *currVertex = NULL;
	int id = 0;

	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			if ( iy == yCntTotalVertex-1 ) { currContact = contacts[0]; }//gate contact
			else if ( iy == 0 ) { currContact = contacts[1]; }//imagined channel contact
			else { currContact = NULL; }

			id = vertexHelper.IdAt(ix, iy);
			currVertex = GetVertex(id);

			if ( currContact != NULL )
			{
				currVertex->SetContact(currContact);
				currContact->AddVertex(currVertex);
			}
		}
	}

	/////////////////////////////////////////////////////////////////////
	//set regions
	////////////////////////////////////////////////////////////////////
	regions.push_back(new FDRegion(cntRegion, FDRegion::Tunneling));
	regions.back()->Mat = &MaterialDB::SiO2;
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::Trapping));
	regions.back()->Mat = &MaterialDB::Si3N4;
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::Blocking));
	regions.back()->Mat = &MaterialDB::SiO2;
	cntRegion++;

	/////////////////////////////////////////////////////////////////////
	//set elements with specified region
	////////////////////////////////////////////////////////////////////
	FDRegion *currRegion = NULL;
	FDVertex *swVertex = NULL;
	FDVertex *seVertex = NULL;
	FDVertex *nwVertex = NULL;
	FDVertex *neVertex = NULL;

	//the number of element always one less than the corresponding vertex number.
	for (int iy = 0; iy != yCntTotalVertex-1; ++iy)
	{
		for (int ix = 0; ix != xCntVertex-1; ++ix)
		{
			swVertex = GetVertex(vertexHelper.IdAt(ix, iy));
			seVertex = GetVertex(vertexHelper.IdAt(ix+1, iy));
			neVertex = GetVertex(vertexHelper.IdAt(ix+1, iy+1));
			nwVertex = GetVertex(vertexHelper.IdAt(ix, iy+1));
			elements.push_back(new FDElement(cntElement, swVertex, seVertex, neVertex, nwVertex));

			//set inner member of element and region
			currRegion = thisRegion(iy);//get current region through the use of a helper class
			GetElement(cntElement)->SetRegion(currRegion);
			currRegion->AddElement(GetElement(cntElement));
			cntElement++;
		}
	}
}

void SimpleONO::setAdjacency()
{
	//set vertex properties
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntTotalVertex);
	FDDomainHelper elementHelper = FDDomainHelper(xCntVertex-1, yCntTotalVertex-1); //element number = vertex number - 1
	int id = 0;
	FDVertex *currVertex = NULL;
	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			id = vertexHelper.IdAt(ix, iy);
			currVertex = GetVertex(id);
			////////////////////////////////
			//set adjacent vertex and length
			//if current vertex doesn't have a west/east/south/north edge, the length is set 0.
			if ( ix-1 >= 0 )
			{
				currVertex->WestVertex = GetVertex(vertexHelper.IdAt(ix-1, iy)); 
				currVertex->WestLength = FDVertex::Distance(currVertex, currVertex->WestVertex);
			}
			else
			{
				currVertex->WestVertex = NULL; 
				currVertex->WestLength = 0;
			}
			if ( ix+1 <= vertexHelper.GetMaxX() )
			{
				currVertex->EastVertex = GetVertex(vertexHelper.IdAt(ix+1, iy));
				currVertex->EastLength = FDVertex::Distance(currVertex, currVertex->EastVertex);
			}
			else									
			{
				currVertex->EastVertex = NULL;
				currVertex->EastLength = 0;
			}
			if ( iy-1 >= 0 )
			{
				currVertex->SouthVertex = GetVertex(vertexHelper.IdAt(ix, iy-1));
				currVertex->SouthLength = FDVertex::Distance(currVertex, currVertex->SouthVertex);
			}
			else
			{
				currVertex->SouthVertex = NULL;
				currVertex->SouthLength = 0;
			}
			if ( iy+1 <= vertexHelper.GetMaxY() )
			{
				currVertex->NorthVertex = GetVertex(vertexHelper.IdAt(ix, iy+1));
				currVertex->NorthLength = FDVertex::Distance(currVertex, currVertex->NorthVertex);
			}
			else
			{
				currVertex->NorthVertex = NULL;
				currVertex->NorthLength = 0;
			}									
			//////////////////////
			//set adjacent element
			if ( (ix-1 >= 0) && (iy-1 >=0) )						{ currVertex->SouthwestElem = GetElement(elementHelper.IdAt(ix-1, iy-1)); }
			else													{ currVertex->SouthwestElem = NULL; }
			if ( (ix < vertexHelper.GetMaxX()) && (iy-1 >= 0) )		{ currVertex->SoutheastElem = GetElement(elementHelper.IdAt(ix, iy-1)); }
			else													{ currVertex->SoutheastElem = NULL; }
			if ( (ix < vertexHelper.GetMaxX()) 
				&& (iy < vertexHelper.GetMaxY()) )					{ currVertex->NortheastElem = GetElement(elementHelper.IdAt(ix, iy)); }
			else													{ currVertex->NortheastElem = NULL;}
			if ( (ix-1 >= 0) && (iy < vertexHelper.GetMaxY()) )		{ currVertex->NorthwestElem = GetElement(elementHelper.IdAt(ix-1, iy)); }
			else													{ currVertex->NorthwestElem = NULL; }
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//set element adjacent properties
	//this is actually already done in the construction of every element.
	FDElement *currElem = NULL;
	for (int iy = 0; iy != yCntVertexTunnel-1; ++iy)
	{
		for (int ix = 0; ix != xCntVertex - 1; ++ix)
		{
			id = elementHelper.IdAt(ix, iy);
			currElem = GetElement(id);

			currElem->SouthwestVertex = GetVertex(vertexHelper.IdAt(ix ,iy));
			currElem->SoutheastVertex = GetVertex(vertexHelper.IdAt(ix+1, iy));
			currElem->NortheastVertex = GetVertex(vertexHelper.IdAt(ix+1, iy+1));
			currElem->NorthwestVertex = GetVertex(vertexHelper.IdAt(ix, iy+1));
		}
	}
}

double SimpleONO::yNextGridLength(int vertexY)
{
	if ( vertexY < yCntVertexTunnel - 1)
		return yGridTunnel;
	else
		vertexY -= yCntVertexTunnel - 1;

	if ( vertexY < yCntVertexTrap - 1)
		return yGridTrap;
	else
		vertexY -= yCntVertexTrap - 1;

	if ( vertexY < yCntVertexBlock - 1)
		return yGridBlock;
	else
		return 0;//0 means that current vertex is the last vertex in this direction
}

double SimpleONO::xNextGridLength(int vertexX)
{
	return xGrid;
}

FDRegion * SimpleONO::thisRegion(int elemY)
{
	if ( elemY < yCntVertexTunnel - 1 ) //transfer the vertex count into element count
		return GetRegion(0);
	else
		elemY -= yCntVertexTunnel - 1;

	if ( elemY < yCntVertexTrap - 1 )
		return GetRegion(1);
	else
		elemY -= yCntVertexTrap - 1;

	if ( elemY < yCntVertexBlock - 1)
		return GetRegion(2);
	else
		return NULL;
}

void SimpleONO::stuffPotential()
{
	//TODO: the voltage of contact should be regarded differently                
	//this value is obtained from Sentaurus result for the current condition
	double channelPotential = 0.6345;
	double elecFieldTunnel = 9.54e6; // in [V/cm]
	double elecFieldTrap = 4.96e6; // in [V/cm]
	double elecFieldBlock = 9.55e6; // in [V/cm]

	Normalization theNorm = Normalization();
	//elecFieldTunnel = theNorm.PushElecField(elecFieldTunnel);
	//elecFieldTrap = theNorm.PushElecField(elecFieldTrap);
	//elecFieldBlock = theNorm.PushElecField(elecFieldBlock);
	//channelPotential = theNorm.PushPotential(channelPotential);

	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntTotalVertex);
	int id = 0;
	FDVertex * currVertex = NULL;
	double potential = channelPotential;
	double nextElecField = 0;
	double iyForElecField = 0;
	double normPotential = 0;

	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			id = vertexHelper.IdAt(ix, iy);
			currVertex = GetVertex(id);
			// the calculated potential is in [V], so normalization is needed here when stuffing.
			normPotential = theNorm.PushPotential(potential);
			currVertex->Phys.SetPhysPrpty(PhysProperty::ElectrostaticPotential, normPotential);
		}
		
		//for next electric field
		iyForElecField = iy;
		if ( iyForElecField < yCntVertexTunnel - 1 )
			nextElecField = elecFieldTunnel;
		else
		{
			iyForElecField -= yCntVertexTunnel - 1;
			if ( iyForElecField < yCntVertexTrap - 1 )
				nextElecField = elecFieldTrap;
			else
			{
				iyForElecField -= yCntVertexTrap - 1;
				if ( iyForElecField < yCntVertexBlock - 1 )
					nextElecField = elecFieldBlock;
				else
					nextElecField = 0;
			}
		}

		double i = yNextGridLength(iy); // the return value of this method is in [cm]
		potential += yNextGridLength(iy) * nextElecField; // [V] = [cm] * [V/cm]
	}
}

void SimpleONO::setVertexPhysics()
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
			currVertex->Phys.FillVertexPhysUsingMatPropty(currVertex, verPrptys.at(iPrpty), matPrptys.at(iPrpty));
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
			currVertex->Phys.SetPhysPrpty(verPrptys.at(iPrpty), physValue);
			*/
		}
	}
}

void SimpleONO::refreshBandEnergy()
{
	double RefPotential = SctmPhys::ReferencePotential;
	double q = SctmPhys::q;

	double energy = 0;
	double affinity = 0;
	double bandgap = 0;
	double potential = 0;

	FDVertex *currVertex = NULL;
	for (size_t iver = 0; iver != vertices.size(); ++iver)
	{
		currVertex = GetVertex(iver);

		potential = currVertex->Phys.GetPhysPrpty(PhysProperty::ElectrostaticPotential);
		affinity = currVertex->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity);
		bandgap = currVertex->Phys.GetPhysPrpty(PhysProperty::Bandgap);

		//Ec = -X-q(phi-phiRef)
		energy = -affinity - q*(potential-RefPotential);
		currVertex->Phys.SetPhysPrpty(PhysProperty::ConductionBandEnergy, energy);
		
		//Ev = -X-q(phi-phiRef)-Eg
		energy = energy - bandgap;
		currVertex->Phys.SetPhysPrpty(PhysProperty::ValenceBandEnergy, energy);
	}
}

void SimpleONO::setBoundaryCondition()
{
	/////////////////////////////////////////////////////////////////////////////
	//set physical class members
	//in [V]
	//TODO: the gate potential should be obtained with gate voltage and work function.
	//Currently, the gate voltage is not considered in the structure.
	/////////////////////////////////////////////////////////////////////////////
	this->gatePotential = 16.526;
	this->channelPotential = 0.634;
	double potentialValue = 0;

	//normalization is needed here because all the values related to the domain details (i.e. the stored value) are normalized
	Normalization theNorm = Normalization();

	FDVertex *currVertex;
	for (std::size_t iVer = 0; iVer != vertices.size(); ++iVer)
	{
		currVertex = GetVertex(iVer);
		
		//firstly, to decide if the vertex is at contact.
		if (currVertex->IsAtContact())
		{
			//the gate name is in accordance with the name specified in setting domain details
			if (currVertex->Contact->ContactName == "Gate")
			{
				potentialValue = theNorm.PushPotential(this->gatePotential);
				//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
				currVertex->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Dirichlet, potentialValue);
			}
			else if (currVertex->Contact->ContactName == "Channel")
			{
				potentialValue = theNorm.PushPotential(this->channelPotential);
				//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
				currVertex->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Dirichlet, potentialValue);
			}
		}
		else //the vertex is not related to a contact
		{
			//check the existence of adjacent element to check if the vertex is boundary vertex
			//the following is used to set the potential boundary conditions
			if ( (currVertex->NorthwestElem == NULL) || (currVertex->NortheastElem == NULL) ||
				 (currVertex->SouthwestElem == NULL) || (currVertex->SoutheastElem == NULL)) 
			{
				//both two boundary condition values can be neglected in artificial boundary conditions
				currVertex->BndCond.SetBndCond(FDBoundary::Potential, FDBoundary::BC_Artificial);
				//currVertex->BndCond.SetBndCond(FDBoundary::eCurrentDensity, FDBoundary::BC_Artificial);
			}

			//the following is used to set the current boundary conditions in simpleONO structures.
			//TODO: for complicated structures, this method is not checked.
			if ( 
				(( currVertex->NortheastElem == NULL ? false : currVertex->NortheastElem->Region->Type == FDRegion::Blocking )
				|| ( currVertex->NorthwestElem == NULL ? false : currVertex->NorthwestElem->Region->Type == FDRegion::Blocking ))
				&& 
				(( currVertex->SoutheastElem == NULL ? false : currVertex->SoutheastElem->Region->Type == FDRegion::Trapping )
				|| ( currVertex->SouthwestElem == NULL ? false : currVertex->SouthwestElem->Region->Type == FDRegion::Trapping) )
			   )
			{
				currVertex->BndCond.SetBndCond(FDBoundary::eCurrentDensity, FDBoundary::BC_Dirichlet);
			}

			if ( 
				(( currVertex->NortheastElem == NULL ? false : currVertex->NortheastElem->Region->Type == FDRegion::Trapping )
				|| ( currVertex->NorthwestElem == NULL ? false : currVertex->NorthwestElem->Region->Type == FDRegion::Trapping ))
				&& 
				(( currVertex->SoutheastElem == NULL ? false : currVertex->SoutheastElem->Region->Type == FDRegion::Tunneling )
				|| ( currVertex->SouthwestElem == NULL ? false : currVertex->SouthwestElem->Region->Type == FDRegion::Tunneling) )
				)
			{
				currVertex->BndCond.SetBndCond(FDBoundary::eCurrentDensity, FDBoundary::BC_Dirichlet);
			}
		}
	}
}

void SimpleONO::setBoundaryCondition(bool fake)
{
	/////////////////////////////////////////////////////////////////////////////
	//set physical class members
	//in [V]
	//TODO: the gate potential should be obtained with gate voltage and work function.
	//Currently, the gate voltage is not considered in the structure.
	/////////////////////////////////////////////////////////////////////////////
	this->gatePotential = 16.526;
	this->channelPotential = 0.634;

	//normalization is needed here because all the values related to the domain details (i.e. the stored value) are normalized
	//Normalization theNorm = Normalization();

	FDVertex *currVertex;
	for (std::size_t iVer = 0; iVer != vertices.size(); ++iVer)
	{
		currVertex = GetVertex(iVer);
		setVertBC_Potential(currVertex);
		setVertBC_eDensity(currVertex);
	}
}

bool SimpleONO::isValidElem(FDElement *elem)
{
	return elem != NULL;
}

bool SimpleONO::isNotElemOf(FDRegion::RegionType rType, FDElement *elem)
{
	if (elem == NULL)
		return true;
	else
	{
		return elem->Region->Type != rType;
	}
}

void SimpleONO::setVertBC_eDensity(FDVertex *vert)
{
	bool notTrapping_NW = isNotElemOf(FDRegion::Trapping, vert->NorthwestElem);
	bool notTrapping_NE = isNotElemOf(FDRegion::Trapping, vert->NortheastElem);
	bool notTrapping_SE = isNotElemOf(FDRegion::Trapping, vert->SoutheastElem);
	bool notTrapping_SW = isNotElemOf(FDRegion::Trapping, vert->SouthwestElem);

	bool valid_NW = isValidElem(vert->NorthwestElem);
	bool valid_NE = isValidElem(vert->NortheastElem);
	bool valid_SE = isValidElem(vert->SoutheastElem);
	bool valid_SW = isValidElem(vert->SouthwestElem);
	
	//Northwest corner
	if ( notTrapping_NW && notTrapping_NE &&
		notTrapping_SW && !notTrapping_SE )
	{
		if (              valid_NE &&
	         !valid_SW )
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, 1));
			return;
		}
		if (              !valid_NE &&
			 valid_SW )
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-vert->EastLength, vert->SouthLength));
		return;
	}

	//Northeast corner
	if ( notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && notTrapping_SE )
	{
		if ( valid_NW &&
							!valid_SE)
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0 ,1));
			return;
		}
		if ( !valid_NW &&
							valid_SE)
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(vert->WestLength, vert->SouthLength));
		return;
	}

	//Southeast corner
	if ( !notTrapping_NW && notTrapping_NE && 
		notTrapping_SW && notTrapping_SE )
	{
		if (			!valid_NE &&
			valid_SW)
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
			return;
		}
		if (			valid_NE &&
			!valid_SW)
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(vert->WestLength, -vert->NorthLength));
		return;
	}

	//Southwest corner
	if ( notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && notTrapping_SE)
	{
		if ( valid_NW &&
						!valid_SE )
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
			return;
		}
		if ( !valid_NW &&
						valid_SE )
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
			return;
		}
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-vert->EastLength, -vert->NorthLength));
		return;
	}

	//North side
	if ( notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && !notTrapping_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, 1));
		return;
	}

	//East side
	if ( !notTrapping_NW && notTrapping_NE && 
		!notTrapping_SW && notTrapping_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
		return;
	}

	//South side
	if ( !notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && notTrapping_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(0, -1));
		return;
	}

	//West side
	if ( notTrapping_NW && !notTrapping_NE && 
		notTrapping_SW && !notTrapping_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(-1, 0));
		return;
	}

}

void SimpleONO::setVertBC_Potential(FDVertex *vert)
{
	//normalization is needed here because all the values related to the domain details (i.e. the stored value) are normalized
	static Normalization theNorm = Normalization();
	static double potentialValue = 0;

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
			potentialValue = theNorm.PushPotential(this->gatePotential);
			//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
			vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Dirichlet, potentialValue);
			return;
		}
		else if (vert->Contact->ContactName == "Channel")
		{
			potentialValue = theNorm.PushPotential(this->channelPotential);
			//the second value has the default value of 0 in setting BC_Dirichlet boundary condition.
			vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Dirichlet, potentialValue);
			return;
		}
	}

	//Northwest corner
	if ( !isValid_NW && !isValid_NE &&
		!isValid_SW && isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-vert->EastLength, vert->SouthLength));
		return;
	}

	//Northeast corner
	if ( !isValid_NW && !isValid_NE && 
		isValid_SW && !isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(vert->WestLength, vert->SouthLength));
		return;
	}

	//Southeast corner
	if ( isValid_NW && !isValid_NE && 
		!isValid_SW && !isValid_SE )
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(vert->WestLength, -vert->NorthLength));
		return;
	}

	//Southwest corner
	if ( !isValid_NW && isValid_NE && 
		!isValid_SW && !isValid_SE)
	{
		//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-vert->EastLength, -vert->NorthLength));
		return;
	}

	//North side
	if ( !isValid_NW && !isValid_NE && 
		isValid_SW && isValid_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(0, 1));
		return;
	}

	//East side
	if ( isValid_NW && !isValid_NE && 
		isValid_SW && !isValid_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(1, 0));
		return;
	}

	//South side
	if ( isValid_NW && isValid_NE && 
		!isValid_SW && !isValid_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(0, -1));
		return;
	}

	//West side
	if ( !isValid_NW && isValid_NE && 
		!isValid_SW && isValid_SE)
	{
		vert->BndCond.SetBndCond(true, FDBoundary::Potential, FDBoundary::BC_Neumann, 0, VectorValue(-1, 0));
		return;
	}

}
