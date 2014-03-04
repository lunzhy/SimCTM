/**
* @file TripleCells.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2014-3-4   11:38
* @note
* @todo
*/

#include "TripleCells.h"
#include "FDDomain.h"
#include "DomainDetails.h"
#include "Normalization.h"
#include "SctmUtils.h"
#include "SctmPhys.h"
#include "Material.h"

using namespace SctmUtils;

TripleCells::TripleCells()
{
	this->temperature = SctmGlobalControl::Get().Temperature;
}

void TripleCells::setDomainDetails()
{
	double nm_in_cm = SctmPhys::nm_in_cm;

	int indexVertex = 0;
	int indexElement = 0;
	int indexRegion = 0;
	int indexContact = 0;

	Normalization norm = Normalization(this->temperature);
	
	double coordX_in_nm = 0; // real value (without normalization), in [nm]
	double coordY_in_nm = 0;
	double normCoordX = 0; // normalized value, in [cm]
	double normCoordY = 0;
	//////////////////////////////////////////////////////////////////
	//set vertices
	//////////////////////////////////////////////////////////////////
	int xIndexMax = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	int yIndexMax = gridThickTunnel + gridThickTrap + gridThickBlock + gridThickIso;
	for (size_t iy = 0; iy <= yIndexMax; ++iy)
	{
		for (size_t ix = 0; ix <= xIndexMax; ++ix)
		{
			//the iteration sequence guarantees the use of FDDomainHelper
			coordX_in_nm = getCoordX(ix, iy);
			coordY_in_nm = getCoordY(ix, iy);
			if (coordX_in_nm < 0 ) //so coordY_in_nm will be less than o
			{
				indexVertex++; 
				//count this vertex as well, so FDDomainHelper can be used to position the vertex by idX and idY.
				//so the indexes of the vertices are not coherent.
				continue;
			}
			normCoordX = norm.PushLength(coordX_in_nm * nm_in_cm);
			normCoordY = norm.PushLength(coordY_in_nm * nm_in_cm);
			vertices.push_back(new FDVertex(indexVertex, normCoordX, normCoordY));
			indexVertex++;
		}
	}

	///////////////////////////////////////////////////////////////////
	//set contacts
	///////////////////////////////////////////////////////////////////
	double voltage = voltageGate1;
	double workfun = workfunctionGate1;
	contacts.push_back(new FDContact(indexContact++, "Gate1", norm.PushPotential(voltage), norm.PushPotential(workfun)));
	voltage = voltageGate2;
	workfun = workfunctionGate2;
	contacts.push_back(new FDContact(indexContact++, "Gate2", norm.PushPotential(voltage), norm.PushPotential(workfun)));
	voltage = voltageGate3;
	workfun = workfunctionGate3;
	contacts.push_back(new FDContact(indexContact++, "Gate3", norm.PushPotential(voltage), norm.PushPotential(workfun)));

	FDDomainHelper vertexHelper = FDDomainHelper(xIndexMax + 1, yIndexMax + 1);
	FDContact *currContact = NULL;
	FDVertex *currVertex = NULL;
	int vertID = 0;
	int indexX = 0;
	int indexY = 0;

	/////gate1
	currContact = this->GetContact("Gate1");
	//gate1, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (size_t ix = 0; ix <= gridWidthGate1; ++ix)
	{
		indexX = ix;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate1, y direction
	indexX = gridWidthGate1;
	for (size_t iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate2
	currContact = this->GetContact("Gate2");
	//gate2, y direction, left side
	indexX = gridWidthGate1 + gridWidthIso2;
	for (size_t iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (size_t ix = 0; ix <= gridWidthGate2; ++ix)
	{
		indexX = gridWidthGate1 + gridWidthIso2 + ix;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, y direction, right side
	indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2;
	for (size_t iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate 3
	currContact = this->GetContact("Gate3");
	//gate3, y direction
	indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso2;
	for (size_t iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate 3, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (size_t ix = 0; ix <= gridWidthGate3; ++ix)
	{
		indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + ix;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////channel
	currContact = this->GetContact("Channel");
	indexY = 0;
	for (size_t ix = 0; ix <= xIndexMax; ++ix)
	{
		indexX = ix;
		vertID = vertexHelper.IdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	///////////////////////////////////////////////////////////////////
	//set regions and elements
	///////////////////////////////////////////////////////////////////
	using MaterialDB::GetMaterial;
	using MaterialDB::Mat;
	Mat::Name currMatName;//ErrorMaterial
	currMatName = Mat::Parse(materialTunnel);
	regions.push_back(new FDRegion(indexRegion++, "Tunnel", GetMaterial(currMatName)));
	currMatName = Mat::Parse(materialTrap);
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate1", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso3", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate3", GetMaterial(currMatName)));
	currMatName = Mat::Parse(materialBlock);
	regions.push_back(new FDRegion(indexRegion++, "Block", GetMaterial(currMatName)));
	currMatName = Mat::Parse(materialIso);
	regions.push_back(new FDRegion(indexRegion++, "Iso2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Iso3", GetMaterial(currMatName)));

	/////set element of each region
	FDRegion *currRegion = NULL;
	int gridTotalWidth = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	int ixBegin = 0;
	int ixEnd = 0; // exclude the end
	int iyBegin = 0;
	int iyEnd = 0;

	//Tunnel
	currRegion = this->GetRegion("Tunnel");
	ixBegin = 0; ixEnd = gridTotalWidth;
	iyBegin = 0; iyEnd = gridThickTunnel;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate1
	currRegion = this->GetRegion("Trap.Gate1");
	ixBegin = 0; ixEnd = gridWidthGate1;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso2
	currRegion = this->GetRegion("Trap.Iso2");
	ixBegin = gridWidthGate1; ixEnd = gridWidthGate1 + gridWidthIso2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate2
	currRegion = this->GetRegion("Trap.Gate2");
	ixBegin = gridWidthGate1 + gridWidthIso2; ixEnd = ixBegin + gridWidthGate2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso3
	currRegion = this->GetRegion("Trap.Iso3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate3
	currRegion = this->GetRegion("Trap.Gate3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3; ixEnd = ixBegin + gridWidthGate3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Block
	currRegion = this->GetRegion("Block");
	ixBegin = 0; ixEnd = gridTotalWidth;
	iyBegin = gridThickTunnel + gridThickTrap; iyEnd = iyBegin + gridThickBlock;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso2
	currRegion = this->GetRegion("Iso2");
	ixBegin = gridWidthGate1; ixEnd = ixBegin + gridWidthIso2;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso3
	currRegion = this->GetRegion("Iso3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement++, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);


}

double TripleCells::getCoordX(int idX, int idY)
{
	double currdX = 0; // the coordinate in x-direction, in [nm], without normalization 

	static double lengthPerGridGate1 = widthGate1 / gridWidthGate1;
	static double lengthPerGridIso2 = widthIso2 / gridWidthIso2;
	static double lengthPerGridGate2 = widthGate2 / gridWidthGate2;
	static double lengthPerGridIso3 = widthIso3 / gridWidthIso3;
	static double lengthPerGridGate3 = widthGate3 / gridWidthGate3;

	int yBndIndexMain = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation

	if (idY > yBndIndexMain) //invalid vertex
	{
		if (((idX >= gridWidthGate1) && (idX <= gridWidthGate1 + gridWidthIso2)) ||
			((idX >= gridWidthGate1 + gridWidthIso2 + gridWidthGate2) && (idX <= gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3)))
		{
			//continue the process
		}
		else
		{
			return -1;// negative values means that no vertex exists at this point.
		}
	}
	// up to here, invalid point is wiped out.
	if (idX <= gridWidthGate1)
	{
		currdX = idX * lengthPerGridGate1;
		return currdX;
	}
	idX -= gridWidthGate1;
	currdX = widthGate1;
	if (idX <= gridWidthIso2)
	{
		currdX += idX * lengthPerGridIso2;
		return currdX;
	}
	idX -= gridWidthIso2;
	currdX += widthIso2;
	if (idX <= gridWidthGate2)
	{
		currdX += idX * lengthPerGridGate2;
		return currdX;
	}
	idX -= gridWidthGate2;
	currdX += widthGate2;
	if (idX <= gridWidthIso3)
	{
		currdX += idX * lengthPerGridIso3;
		return currdX;
	}
	idX -= gridWidthIso3;
	currdX += widthIso3;
	// idX <= gridWidthGate3 is always true, if the checking runs to here.
	currdX += idX * lengthPerGridGate3;
	return currdX;
}

double TripleCells::getCoordY(int idX, int idY)
{
	double currdY = 0; // the coordinate in y-direction, in [nm], without normalization
	static double lengthPerGridTunnel = thickTunnel / gridThickTunnel;
	static double lengthPerGridTrap = thickTrap / gridThickTrap;
	static double lengthPerGridBlock = thickBlock / gridThickBlock;
	static double lengthPerGridIso = thickIso / gridThickIso;

	int yBndIndexMain = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation

	if (idY > yBndIndexMain) //invalid vertex
	{
		if (((idX >= gridWidthGate1) && (idX <= gridWidthGate1 + gridWidthIso2)) ||
			((idX >= gridWidthGate1 + gridWidthIso2 + gridWidthGate2) && (idX <= gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3)))
		{
			//continue the process
		}
		else
		{
			return -1;// negative values means that no vertex exists at this point.
		}
	}
	// up to here, invalid point is wiped out.
	if (idY <= gridThickTunnel)
	{
		currdY = idY * lengthPerGridTunnel;
		return currdY;
	}
	idY -= gridThickTunnel;
	currdY = thickTunnel;
	if (idY <= gridThickTrap)
	{
		currdY += idY * lengthPerGridTrap;
		return currdY;
	}
	idY -= gridThickTrap;
	currdY += thickTrap;
	if (idY <= gridThickBlock)
	{
		currdY += idY * lengthPerGridBlock;
		return currdY;
	}
	idY -= gridThickBlock;
	currdY += thickBlock;
	//up to here, idY <= gridThickIso is always true
	currdY += idY * lengthPerGridIso;
	return currdY;
}

void TripleCells::setSingleElement(int idElem, FDRegion *region, int xbegin, int xend, int ybegin, int yend)
{
	static int gridTotalX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	static int gridTotalY = gridThickTunnel + gridThickTrap + gridThickBlock + gridThickIso;
	static FDDomainHelper vertexHelper = FDDomainHelper(gridTotalX + 1, gridTotalY + 1);

	FDElement *currElem = NULL;
	FDVertex *swVertex = NULL;
	FDVertex *seVertex = NULL;
	FDVertex *nwVertex = NULL;
	FDVertex *neVertex = NULL;

	for (size_t iyElem = ybegin; iyElem != yend; ++iyElem)
	{
		for (size_t ixElem = xbegin; ixElem != xend; ++ixElem)
		{
			swVertex = GetVertex(vertexHelper.IdAt(ixElem, iyElem));
			seVertex = GetVertex(vertexHelper.IdAt(ixElem + 1, iyElem));
			neVertex = GetVertex(vertexHelper.IdAt(ixElem + 1, iyElem + 1));
			nwVertex = GetVertex(vertexHelper.IdAt(ixElem, iyElem + 1));
			currElem = new FDElement(idElem, swVertex, seVertex, neVertex, nwVertex);
			this->elements.push_back(currElem);
			currElem->SetVertexAdjacent();
			currElem->SetRegion(region);
			region->AddElement(currElem);
		}
	}
}

