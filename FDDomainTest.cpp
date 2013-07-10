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

#include "FDDomainTest.h"
#include "Material.h"
#include "Normalization.h"
#include <iostream>

void FDDomainTest::BuildDomain()
{
	//Initialize the vectors
	vertices.clear();
	elements.clear();
	regions.clear();

	prepareStructures(); //mainly set the structure and mesh parameters
	setDomainDetails(); //fill in the regions, vertices and elements
	setAdjacency(); //set the adjacency of vertices and elements

	printStructure();
}

void FDDomainTest::prepareStructures()
{
	setParameters();
}

void FDDomainTest::setParameters()
{
	double nm_in_cm = GeneralMath::nm_in_cm;

	////////////////////////////////////////////////////////////////////////////
	//modify here to change the structures
	double xLength_in_nm = 10;
	double yLengthTunnel_in_nm = 10;
	double yLengthTrap_in_nm = 10;
	double yLengthBlock_in_nm = 10;
	int xGridNumber = 5;
	int yGridNumberTunnel = 5;
	int yGridNumberTrap = 5;
	int yGridNumberBlock = 5;
	////////////////////////////////////////////////////////////////////////////

	xLength = xLength_in_nm * nm_in_cm;
	xCntVertex = xGridNumber + 1;
	yLengthTunnel = yLengthTunnel_in_nm * nm_in_cm;
	yCntVertexTunnel = yGridNumberTunnel + 1;
	yLengthTrap = yLengthTrap_in_nm * nm_in_cm;
	yCntVertexTrap = yGridNumberTrap + 1;
	yLengthBlock = yLengthBlock_in_nm * nm_in_cm;
	yCntVertexBlock = yGridNumberBlock + 1;

	xGrid = xLength / ( xCntVertex - 1 );
	yGridTunnel = yLengthTunnel / ( yCntVertexTunnel - 1 );
	yGridTrap = yLengthTrap / ( yCntVertexTrap - 1 );
	yGridBlock = yLengthBlock / ( yCntVertexBlock - 1 );
}

void FDDomainTest::printStructure()
{
	for (std::vector<FDVertex *>::size_type ix = 0; ix != this->vertices.size(); ++ix)
	{
		std::cout << "id=" << vertices.at(ix)->GetInternalID() 
			<< '\t' << (vertices.at(ix)->WestVertex == NULL ? -1 : vertices.at(ix)->WestVertex->GetInternalID())
			<< '\t' << (vertices.at(ix)->SouthVertex == NULL ? -1 : vertices.at(ix)->SouthVertex->GetInternalID())
			<< '\t' << (vertices.at(ix)->EastVertex == NULL ? -1 : vertices.at(ix)->EastVertex->GetInternalID())
			<< '\t' << (vertices.at(ix)->NorthVertex == NULL ? -1 : vertices.at(ix)->NorthVertex->GetInternalID())
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
		std::cout << "id=" << elements.at(ix)->GetInternalID()
			<< '\t' << elements.at(ix)->SouthwestVertex->GetInternalID()
			<< '\t' << elements.at(ix)->SoutheastVertex->GetInternalID()
			<< '\t' << elements.at(ix)->NortheastVertex->GetInternalID()
			<< '\t' << elements.at(ix)->NorthwestVertex->GetInternalID()
			<< '\t' << elements.at(ix)->WestLength
			<< '\t' << elements.at(ix)->SouthLength
			<< '\t' << elements.at(ix)->EastLength
			<< '\t' << elements.at(ix)->NorthLength
			<< std::endl;
	}
}

void FDDomainTest::setDomainDetails()
{
	int cntVertex = 0;
	int cntElement = 0;
	int cntRegion = 0;

	double currCoordX = 0.0;
	double currCoordY = 0.0;

	Utility::Normalization theNorm = Utility::Normalization();
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntVertexTunnel);
	FDDomainHelper elementHelper = FDDomainHelper(xCntVertex-1, yCntVertexTunnel-1);

	//set vertices
	int yCntTotalVertex = yCntVertexTunnel + yCntVertexTrap + yCntVertexBlock - 1 - 1;
	double normCoordX = 0.0;
	double normCoordY = 0.0;
	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		//currCoordY = iy * yGridTunnel;
		//currCoordY = theNorm.PushLength(currCoordY);
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			//currCoordX = ix * xGrid;
			normCoordX = theNorm.PushLength(currCoordX);
			normCoordY = theNorm.PushLength(currCoordY);
			vertices.push_back(new FDVertex(cntVertex, normCoordX, normCoordY));
			cntVertex++;
			currCoordX += xNextGridLength(ix);
		}
		currCoordX = 0;
		currCoordY += yNextGridLength(iy);
	}

	//set regions
	regions.push_back(new FDRegion(cntRegion, FDRegion::TunnelingOxide));
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::TrappingLayer));
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::BlockingOxide));
	cntRegion++;
	//set elements with specified region
	FDRegion *currRegion = NULL;
	FDVertex *swVertex = NULL;
	FDVertex *seVertex = NULL;
	FDVertex *nwVertex = NULL;
	FDVertex *neVertex = NULL;

	for (int iy = 0; iy != yCntVertexTunnel-1; ++iy)
	{
		for (int ix = 0; ix != xCntVertex-1; ++ix)
		{
			swVertex = getVertex(vertexHelper.IdAt(ix, iy));
			seVertex = getVertex(vertexHelper.IdAt(ix+1, iy));
			neVertex = getVertex(vertexHelper.IdAt(ix+1, iy+1));
			nwVertex = getVertex(vertexHelper.IdAt(ix, iy+1));
			elements.push_back(new FDElement(cntElement, swVertex, seVertex, neVertex, nwVertex));
			//set inner member of element and region
			currRegion = thisRegion(iy);
			getElement(cntElement)->SetRegion(currRegion);
			currRegion->AddElement(getElement(cntElement));
			cntElement++;
		}
	}
}

void FDDomainTest::setAdjacency()
{
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntVertexTunnel);
	FDDomainHelper elementHelper = FDDomainHelper(xCntVertex-1, yCntVertexTunnel-1);
	//set vertex properties
	int id = 0;
	FDVertex *currVertex = NULL;
	for (int iy = 0; iy != yCntVertexTunnel; ++iy)
	{
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			id = vertexHelper.IdAt(ix, iy);
			currVertex = getVertex(id);
			/////////////////////
			//set adjacent vertex and length
			if ( ix-1 >= 0 )
			{
				currVertex->WestVertex = getVertex(vertexHelper.IdAt(ix-1, iy)); 
				currVertex->WestLength = FDVertex::Distance(currVertex, currVertex->WestVertex);
			}
			else
			{
				currVertex->WestVertex = NULL; 
				currVertex->WestLength = 0;
			}
			if ( ix+1 <= vertexHelper.GetMaxX() )
			{
				currVertex->EastVertex = getVertex(vertexHelper.IdAt(ix+1, iy));
				currVertex->EastLength = FDVertex::Distance(currVertex, currVertex->EastVertex);
			}
			else									
			{
				currVertex->EastVertex = NULL;
				currVertex->EastLength = 0;
			}
			if ( iy-1 >= 0 )
			{
				currVertex->SouthVertex = getVertex(vertexHelper.IdAt(ix, iy-1));
				currVertex->SouthLength = FDVertex::Distance(currVertex, currVertex->SouthVertex);
			}
			else
			{
				currVertex->SouthVertex = NULL;
				currVertex->SouthLength = 0;
			}
			if ( iy+1 <= vertexHelper.GetMaxY() )
			{
				currVertex->NorthVertex = getVertex(vertexHelper.IdAt(ix, iy+1));
				currVertex->NorthLength = FDVertex::Distance(currVertex, currVertex->NorthVertex);
			}
			else
			{
				currVertex->NorthVertex = NULL;
				currVertex->NorthLength = 0;
			}									
			/////////////////////
			//set adjacent element
			if ( (ix-1 >= 0) && (iy-1 >=0) )						{ currVertex->SouthwestElem = getElement(elementHelper.IdAt(ix-1, iy-1)); }
			else													{ currVertex->SouthwestElem = NULL; }
			if ( (ix < vertexHelper.GetMaxX()) && (iy-1 >= 0) )		{ currVertex->SoutheastElem = getElement(elementHelper.IdAt(ix, iy-1)); }
			else													{ currVertex->SoutheastElem = NULL; }
			if ( (ix < vertexHelper.GetMaxX()) 
				&& (iy < vertexHelper.GetMaxY()) )					{ currVertex->NortheastElem = getElement(elementHelper.IdAt(ix, iy)); }
			else													{ currVertex->NortheastElem = NULL;}
			if ( (ix-1 >= 0) && (iy < vertexHelper.GetMaxY()) )		{ currVertex->NorthwestElem = getElement(elementHelper.IdAt(ix-1, iy)); }
			else													{ currVertex->NorthwestElem = NULL; }
		}
	}

	//set element adjacent properties
	FDElement *currElem = NULL;
	for (int iy = 0; iy != yCntVertexTunnel-1; ++iy)
	{
		for (int ix = 0; ix != xCntVertex - 1; ++ix)
		{
			id = elementHelper.IdAt(ix, iy);
			currElem = getElement(id);

			currElem->SouthwestVertex = getVertex(vertexHelper.IdAt(ix ,iy));
			currElem->SoutheastVertex = getVertex(vertexHelper.IdAt(ix+1, iy));
			currElem->NortheastVertex = getVertex(vertexHelper.IdAt(ix+1, iy+1));
			currElem->NorthwestVertex = getVertex(vertexHelper.IdAt(ix, iy+1));
		}
	}
}

double FDDomainTest::yNextGridLength(int iy)
{
	if ( iy < yCntVertexTunnel - 1)
		return yGridTunnel;
	else
		iy -= yCntVertexTunnel - 1;

	if ( iy < yCntVertexTrap - 1)
		return yGridTrap;
	else
		iy -= yCntVertexTrap - 1;

	if ( iy < yCntVertexBlock - 1)
		return yGridBlock;
	else
		return 0;//0 means that current vertex is the last vertex in this direction
}

double FDDomainTest::xNextGridLength(int ix)
{
	return xGrid;
}

FDRegion * FDDomainTest::thisRegion(int elemY)
{
	if ( elemY < yCntVertexTunnel - 1 )
		return getRegion(0);
	else
		elemY -= yCntVertexTunnel - 1;

	if ( elemY < yCntVertexTrap - 1 )
		return getRegion(1);
	else
		elemY -= yCntVertexTrap - 1;

	if ( elemY < yCntVertexBlock - 1)
		return getRegion(2);
	else
		return NULL;
}