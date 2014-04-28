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

void TripleCells::buildStructure()
{
	setParametersFromParamParser();
	setDomainDetails();
	setAdjacency();
}

void TripleCells::setParametersFromParamParser()
{
	ParamBase *parBase = NULL;
	
	//gate 1 voltage
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_voltage);
	voltageGate1 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 1 work function
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_workfunction);
	workfunctionGate1 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 2 voltage
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_voltage);
	voltageGate2 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 2 work function
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_workfunction);
	workfunctionGate2 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 3 voltage
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_voltage);
	voltageGate3 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 3 work function
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_workfunction);
	workfunctionGate3 = dynamic_cast<Param<double> *>(parBase)->Value();

	//isolation layer material name
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso_material);
	matNameStrIso = dynamic_cast<Param<string> *>(parBase)->Value();

	//tunnel layer material name
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_tunnel_material);
	matNameStrTunnel = dynamic_cast<Param<string> *>(parBase)->Value();

	//trap layer material name
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_trap_material);
	matNameStrTrap = dynamic_cast<Param<string> *>(parBase)->Value();

	//block layer material name
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_block_material);
	matNameStrBlock = dynamic_cast<Param<string> *>(parBase)->Value();

	//gate 1 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_width);
	widthGate1 = dynamic_cast<Param<double> *>(parBase)->Value();

	//isolation 2 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso2_width);
	widthIso2 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 2 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_width);
	widthGate2 = dynamic_cast<Param<double> *>(parBase)->Value();

	//isolation 3 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso3_width);
	widthIso3 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 3 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_width);
	widthGate3 = dynamic_cast<Param<double> *>(parBase)->Value();

	//gate 1 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_width_grid);
	gridWidthGate1 = dynamic_cast<Param<int> *>(parBase)->Value();

	//isolation 2 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso2_width_grid);
	gridWidthIso2 = dynamic_cast<Param<int> *>(parBase)->Value();

	//gate 2 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_width_grid);
	gridWidthGate2 = dynamic_cast<Param<int> *>(parBase)->Value();

	//isolation 3 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso3_width_grid);
	gridWidthIso3 = dynamic_cast<Param<int> *>(parBase)->Value();

	//gate 3 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_width_grid);
	gridWidthGate3 = dynamic_cast<Param<int> *>(parBase)->Value();

	//isolation layer thick
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso_thick);
	thickIso = dynamic_cast<Param<double> *>(parBase)->Value();

	//block layer thick
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_block_thick);
	thickBlock = dynamic_cast<Param<double> *>(parBase)->Value();

	//trap layer thick
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_trap_thick);
	thickTrap = dynamic_cast<Param<double> *>(parBase)->Value();

	//tunnel layer thick
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_tunnel_thick);
	thickTunnel = dynamic_cast<Param<double> *>(parBase)->Value();

	//grid number of isolation in thick direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso_thick_grid);
	gridThickIso = dynamic_cast<Param<int> *>(parBase)->Value();

	//grid number of block layer in thick direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_block_thick_grid);
	gridThickBlock = dynamic_cast<Param<int> *>(parBase)->Value();

	//grid number of trap layer in thick direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_trap_thick_grid);
	gridThickTrap = dynamic_cast<Param<int> *>(parBase)->Value();

	//grid number of tunnel layer in thick direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_tunnel_thick_grid);
	gridThickTunnel = dynamic_cast<Param<int> *>(parBase)->Value();

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
	//the iteration sequence determined the method to get id width specific (idx, idy)
	for (int iy = 0; iy <= yIndexMax; ++iy)
	{
		for (int ix = 0; ix <= xIndexMax; ++ix)
		{
			//the iteration sequence guarantees the use of FDDomainHelper
			coordX_in_nm = getVertCoordX(ix, iy);
			coordY_in_nm = getVertCoordY(ix, iy);
			if (coordX_in_nm < 0 ) //so coordY_in_nm will be less than o
			{
				//indexVertex++; 
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
	contacts.push_back(new FDContact(indexContact++, "Channel", 0, 0));//explicitly set voltage and workfunction to be 0 

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
	for (int ix = 0; ix <= gridWidthGate1; ++ix)
	{
		indexX = ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate1, y direction
	indexX = gridWidthGate1;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate2
	currContact = this->GetContact("Gate2");
	//gate2, y direction, left side
	indexX = gridWidthGate1 + gridWidthIso2;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (int ix = 0; ix <= gridWidthGate2; ++ix)
	{
		indexX = gridWidthGate1 + gridWidthIso2 + ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, y direction, right side
	indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate 3
	currContact = this->GetContact("Gate3");
	//gate3, y direction
	indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate 3, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (int ix = 0; ix <= gridWidthGate3; ++ix)
	{
		indexX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////channel
	currContact = this->GetContact("Channel");
	indexY = 0;
	for (int ix = 0; ix <= xIndexMax; ++ix)
	{
		indexX = ix;
		vertID = getVertIdAt(indexX, indexY);
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
	currMatName = Mat::Parse(matNameStrTunnel);
	regions.push_back(new FDRegion(indexRegion++, "Tunnel", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrTrap);
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate1", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso3", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate3", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrBlock);
	regions.push_back(new FDRegion(indexRegion++, "Block", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrIso);
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
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate1
	currRegion = this->GetRegion("Trap.Gate1");
	ixBegin = 0; ixEnd = gridWidthGate1;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso2
	currRegion = this->GetRegion("Trap.Iso2");
	ixBegin = gridWidthGate1; ixEnd = gridWidthGate1 + gridWidthIso2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate2
	currRegion = this->GetRegion("Trap.Gate2");
	ixBegin = gridWidthGate1 + gridWidthIso2; ixEnd = ixBegin + gridWidthGate2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso3
	currRegion = this->GetRegion("Trap.Iso3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate3
	currRegion = this->GetRegion("Trap.Gate3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3; ixEnd = ixBegin + gridWidthGate3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Block
	currRegion = this->GetRegion("Block");
	ixBegin = 0; ixEnd = gridTotalWidth;
	iyBegin = gridThickTunnel + gridThickTrap; iyEnd = iyBegin + gridThickBlock;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso2
	currRegion = this->GetRegion("Iso2");
	ixBegin = gridWidthGate1; ixEnd = ixBegin + gridWidthIso2;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso3
	currRegion = this->GetRegion("Iso3");
	ixBegin = gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);


}

double TripleCells::getVertCoordX(int idX, int idY)
{
	double currdX = 0; // the coordinate in x-direction, in [nm], without normalization 

	static double lengthPerGridGate1 = widthGate1 / gridWidthGate1;
	static double lengthPerGridIso2 = widthIso2 / gridWidthIso2;
	static double lengthPerGridGate2 = widthGate2 / gridWidthGate2;
	static double lengthPerGridIso3 = widthIso3 / gridWidthIso3;
	static double lengthPerGridGate3 = widthGate3 / gridWidthGate3;

	if (!isValidVertex(idX, idY))
	{
		return -1;
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

double TripleCells::getVertCoordY(int idX, int idY)
{
	double currdY = 0; // the coordinate in y-direction, in [nm], without normalization
	static double lengthPerGridTunnel = thickTunnel / gridThickTunnel;
	static double lengthPerGridTrap = thickTrap / gridThickTrap;
	static double lengthPerGridBlock = thickBlock / gridThickBlock;
	static double lengthPerGridIso = thickIso / gridThickIso;

	if (!isValidVertex(idX, idY))
	{
		return -1;
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

void TripleCells::setSingleElement(int &idElem, FDRegion *region, int xbegin, int xend, int ybegin, int yend)
{
	FDElement *currElem = NULL;
	FDVertex *swVertex = NULL;
	FDVertex *seVertex = NULL;
	FDVertex *nwVertex = NULL;
	FDVertex *neVertex = NULL;

	for (int iyElem = ybegin; iyElem != yend; ++iyElem)
	{
		for (int ixElem = xbegin; ixElem != xend; ++ixElem)
		{
			swVertex = GetVertex(getVertIdAt(ixElem, iyElem));
			seVertex = GetVertex(getVertIdAt(ixElem + 1, iyElem));
			neVertex = GetVertex(getVertIdAt(ixElem + 1, iyElem + 1));
			nwVertex = GetVertex(getVertIdAt(ixElem, iyElem + 1));
			currElem = new FDElement(idElem++, swVertex, seVertex, neVertex, nwVertex);
			this->elements.push_back(currElem);
			currElem->SetVertexAdjacent();
			currElem->SetRegion(region);
			region->AddElement(currElem);
		}
	}
}

bool TripleCells::isValidVertex(int idX, int idY)
{
	bool ret = true;
	static int gridTotalX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	static int gridMainY = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation
	static int gridTotalY = gridMainY + gridThickIso;
	if (idX < 0 || idX > gridTotalX)
	{
		ret = false;
	}
	if (idY < 0 || idY > gridTotalY)
	{
		ret = false;
	}
	if (idY > gridMainY)
	{
		//check gate1 region
		if (idX >= 0 && idX < gridWidthGate1)
		{
			ret = false;
		}
		//check gate2 region
		if (idX > gridWidthGate1 + gridWidthIso2 && idX < gridWidthGate1 + gridWidthIso2 + gridWidthGate2)
		{
			ret = false;
		}
		//check gate3 region
		if (idX > gridTotalX - gridWidthGate3 && idX <= gridTotalX)
		{
			ret = false;
		}
	}
	return ret;

}

void TripleCells::setAdjacency()
{
	int gridTotalX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	int gridTotalY = gridThickTunnel + gridThickTrap + gridThickBlock + gridThickIso;
	
	FDVertex *currVert = NULL;
	FDVertex *adjacentVert = NULL;

	int vertID = 0;
	int adjacentID = 0;
	for (int iy = 0; iy <= gridTotalY; ++iy)
	{
		for (int ix = 0; ix <= gridTotalX; ++ix)
		{
			vertID = getVertIdAt(ix, iy);
			if (vertID < 0)//(ix, iy) does not corresponds to a valid vertex
			{
				continue;
			}
			currVert = this->GetVertex(vertID);
			//west vertex
			adjacentID = getVertIdAt(ix - 1, iy);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->WestVertex = adjacentVert;
				currVert->WestLength = FDVertex::Distance(currVert, currVert->WestVertex);
			}
			//north vertex
			adjacentID = getVertIdAt(ix, iy + 1);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->NorthVertex = adjacentVert;
				currVert->NorthLength = FDVertex::Distance(currVert, currVert->NorthVertex);
			}
			//east vertex
			adjacentID = getVertIdAt(ix + 1, iy);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->EastVertex = adjacentVert;
				currVert->EastLength = FDVertex::Distance(currVert, currVert->EastVertex);
			}
			//south vertex
			adjacentID = getVertIdAt(ix, iy - 1);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->SouthVertex = adjacentVert;
				currVert->SouthLength = FDVertex::Distance(currVert, currVert->SouthVertex);
			}
		}
	}
}

void TripleCells::setTrapDistribution()
{
	//TODO: Setting the distribution of trap density is temporarily considered here
	//This is only for uniform distribution
	using SctmPhys::TrapProperty;

	Normalization norm = Normalization(this->temperature);
	double unifromTrapDens = SctmGlobalControl::Get().UniformTrapDens;
	double eTrapDens = norm.PushDensity(unifromTrapDens); // in [cm^-3]

	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != ddVerts.size(); ++iVert)
	{
		currVert = ddVerts.at(iVert);
		currVert->Trap->SetTrapPrpty(TrapProperty::eTrapDensity, eTrapDens);
	}
}

int TripleCells::getVertIdAt(int idX, int idY)
{
	int vertID = 0;
	static int gridTotalX = gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	static int gridMainY = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation
	static int gridTotalY = gridMainY + gridThickIso;

	//out of boundary
	if (idX < 0 || idX > gridTotalX)
	{
		return -1;
	}
	if (idY < 0 || idY > gridTotalY)
	{
		return -1;
	}

	if (idY <= gridMainY)// the vertex is in the main part
	{
		vertID = (gridTotalX + 1)*idY + idX;
	}
	else// the vertex is in the isolation layer 
	{
		vertID = (gridTotalX + 1)*(gridMainY + 1);//count the main part
		vertID += (gridWidthIso2 + 1 + gridWidthIso3 + 1)*(idY - gridMainY - 1);//count the isolation vertex
		if (idX < gridWidthGate1)
		{
			return -1;//-1 means invalid id
		}
		else if (idX <= gridWidthGate1 + gridWidthIso2)//in isolation 1
		{
			vertID += idX - gridWidthGate1;
		}
		else if (idX < gridWidthGate1 + gridWidthIso2 + gridWidthGate2)
		{
			return -1;
		}
		else if (idX <= gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3)//in isolation 2
		{
			vertID += gridWidthIso2 + 1;
			vertID += idX - (gridWidthGate1 + gridWidthIso2 + gridWidthGate2);
		}
		else
		{
			return -1;
		}
	}
	return vertID;
}

void TripleCells::postProcessOfDomain()
{
	clearTrapExceptMainCell();
	//writing channel points is done in Pytaurus
	//writeChannelPoints();
}

void TripleCells::writeChannelPoints()
{
	FDContact *channelCont = NULL;
	channelCont = this->GetContact("Channel");

	SctmData::Get().WriteVertexInfo(channelCont->GetContactVerts());
}

void TripleCells::clearTrapExceptMainCell()
{
	string mainCell = "Trap.Gate2";
	bool relatedToMain = false;
	FDVertex *vert = NULL;
	FDElement *northWestElem = NULL;
	FDElement *northEastElem = NULL;
	FDElement *southEastElem = NULL;
	FDElement *southWestElem = NULL;

	for (size_t iVert = 0; iVert != this->ddVerts.size(); ++iVert)
	{
		relatedToMain = false;
		vert = ddVerts.at(iVert);
		northWestElem = vert->NorthwestElem;
		if (northWestElem != NULL && northWestElem->Region->RegName == mainCell)
		{
			relatedToMain = true;
		}
		northEastElem = vert->NortheastElem;
		if (northEastElem != NULL && northEastElem->Region->RegName == mainCell)
		{
			relatedToMain = true;
		}
		southEastElem = vert->SoutheastElem;
		if (southEastElem != NULL && southEastElem->Region->RegName == mainCell)
		{
			relatedToMain = true;
		}
		southWestElem = vert->SouthwestElem;
		if (southWestElem != NULL && southWestElem->Region->RegName == mainCell)
		{
			relatedToMain = true;
		}
		if (!relatedToMain)
		{
			vert->Trap->SetTrapPrpty(TrapProperty::eTrapped, 0);
		}
	}
}

bool TripleCells::IsEndOfEffectiveCapacitor(FDVertex *vert)
{
	FDElement *elem = NULL;

	elem = vert->NorthwestElem;
	if (elem != NULL && (elem->Region->RegName == "Iso2" || elem->Region->RegName == "Iso3"))
	{
		return true;
	}
	elem = vert->NortheastElem;
	if (elem != NULL && (elem->Region->RegName == "Iso2" || elem->Region->RegName == "Iso3"))
	{
		return true;
	}
	elem = vert->SouthwestElem;
	if (elem != NULL && (elem->Region->RegName == "Iso2" || elem->Region->RegName == "Iso3"))
	{
		return true;
	}
	elem = vert->SoutheastElem;
	if (elem != NULL && (elem->Region->RegName == "Iso2" || elem->Region->RegName == "Iso3"))
	{
		return true;
	}
	return false;
}



void TripleCellsFull::setAdditionalParameters()
{
	setParametersFromParamParser();
	ParamBase *parBase = NULL;

	//isolation 1 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso1_width);
	widthIso1 = dynamic_cast<Param<double> *>(parBase)->Value();

	//isolation 1 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso1_width_grid);
	gridWidthIso1 = dynamic_cast<Param<int> *>(parBase)->Value();

	//isolation 4 width
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso4_width);
	widthIso4 = dynamic_cast<Param<double> *>(parBase)->Value();

	//isolation 4 grid number in width direction
	parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_iso4_width_grid);
	gridWidthIso4 = dynamic_cast<Param<int> *>(parBase)->Value();
}

void TripleCellsFull::buildStructure()
{
	setAdditionalParameters();
	setDomainDetails();
	setAdjacency();
}

void TripleCellsFull::setDomainDetails()
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
	int xIndexMax = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + 
		gridWidthIso3 + gridWidthGate3 + gridWidthIso4;
	int yIndexMax = gridThickTunnel + gridThickTrap + gridThickBlock + gridThickIso;
	//the iteration sequence determined the method to get id width specific (idx, idy)
	for (int iy = 0; iy <= yIndexMax; ++iy)
	{
		for (int ix = 0; ix <= xIndexMax; ++ix)
		{
			//the iteration sequence guarantees the use of FDDomainHelper
			coordX_in_nm = getVertCoordX(ix, iy);
			coordY_in_nm = getVertCoordY(ix, iy);
			if (coordX_in_nm < 0) //so coordY_in_nm will be less than o
			{
				//indexVertex++; 
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
	contacts.push_back(new FDContact(indexContact++, "Channel", 0, 0));//explicitly set voltage and workfunction to be 0 

	FDDomainHelper vertexHelper = FDDomainHelper(xIndexMax + 1, yIndexMax + 1);
	FDContact *currContact = NULL;
	FDVertex *currVertex = NULL;
	int vertID = 0;
	int indexX = 0;
	int indexY = 0;

	/////gate1
	currContact = this->GetContact("Gate1");
	//gate1, y direction, left side
	indexX = gridWidthIso1;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate1, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (int ix = 0; ix <= gridWidthGate1; ++ix)
	{
		indexX = gridWidthIso1 + ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate1, y direction, right side
	indexX = gridWidthIso1 + gridWidthGate1;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate2
	currContact = this->GetContact("Gate2");
	//gate2, y direction, left side
	indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (int ix = 0; ix <= gridWidthGate2; ++ix)
	{
		indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate2, y direction, right side
	indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////gate 3
	currContact = this->GetContact("Gate3");
	//gate3, y direction, left side
	indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate 3, x direction
	indexY = gridThickTunnel + gridThickTrap + gridThickBlock;
	for (int ix = 0; ix <= gridWidthGate3; ++ix)
	{
		indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + ix;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}
	//gate3, y direction, right side
	indexX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3;
	for (int iy = 1; iy <= gridThickIso; ++iy)
	{
		indexY = gridThickTunnel + gridThickTrap + gridThickBlock + iy;
		vertID = getVertIdAt(indexX, indexY);
		currVertex = this->GetVertex(vertID);
		currVertex->SetContact(currContact);
		currContact->AddVertex(currVertex);
	}

	/////channel
	currContact = this->GetContact("Channel");
	indexY = 0;
	for (int ix = 0; ix <= xIndexMax; ++ix)
	{
		indexX = ix;
		vertID = getVertIdAt(indexX, indexY);
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
	currMatName = Mat::Parse(matNameStrTunnel);
	regions.push_back(new FDRegion(indexRegion++, "Tunnel", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrTrap);
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso1", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate1", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso3", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Gate3", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Trap.Iso4", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrBlock);
	regions.push_back(new FDRegion(indexRegion++, "Block", GetMaterial(currMatName)));
	currMatName = Mat::Parse(matNameStrIso);
	regions.push_back(new FDRegion(indexRegion++, "Iso1", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Iso2", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Iso3", GetMaterial(currMatName)));
	regions.push_back(new FDRegion(indexRegion++, "Iso4", GetMaterial(currMatName)));

	/////set element of each region
	FDRegion *currRegion = NULL;
	int gridTotalWidth = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + 
						gridWidthIso3 + gridWidthGate3 + gridWidthIso4;
	int ixBegin = 0;
	int ixEnd = 0; // exclude the end
	int iyBegin = 0;
	int iyEnd = 0;

	//Tunnel
	currRegion = this->GetRegion("Tunnel");
	ixBegin = 0; ixEnd = gridTotalWidth;
	iyBegin = 0; iyEnd = gridThickTunnel;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso1
	currRegion = this->GetRegion("Trap.Iso1");
	ixBegin = 0; ixEnd = gridWidthIso1;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate1
	currRegion = this->GetRegion("Trap.Gate1");
	ixBegin = gridWidthIso1; ixEnd = ixBegin + gridWidthGate1;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso2
	currRegion = this->GetRegion("Trap.Iso2");
	ixBegin = gridWidthIso1 + gridWidthGate1; ixEnd = ixBegin + gridWidthIso2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate2
	currRegion = this->GetRegion("Trap.Gate2");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2; ixEnd = ixBegin + gridWidthGate2;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso3
	currRegion = this->GetRegion("Trap.Iso3");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Gate3
	currRegion = this->GetRegion("Trap.Gate3");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3; ixEnd = ixBegin + gridWidthGate3;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Trap.Iso4
	currRegion = this->GetRegion("Trap.Iso4");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3; ixEnd = ixBegin + gridWidthIso4;
	iyBegin = gridThickTunnel; iyEnd = iyBegin + gridThickTrap;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Block
	currRegion = this->GetRegion("Block");
	ixBegin = 0; ixEnd = gridTotalWidth;
	iyBegin = gridThickTunnel + gridThickTrap; iyEnd = iyBegin + gridThickBlock;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso1
	currRegion = this->GetRegion("Iso1");
	ixBegin = 0; ixEnd = ixBegin + gridWidthIso1;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso2
	currRegion = this->GetRegion("Iso2");
	ixBegin = gridWidthIso1 + gridWidthGate1; ixEnd = ixBegin + gridWidthIso2;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso3
	currRegion = this->GetRegion("Iso3");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2; ixEnd = ixBegin + gridWidthIso3;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);

	//Iso4
	currRegion = this->GetRegion("Iso4");
	ixBegin = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3; ixEnd = ixBegin + gridWidthIso4;
	iyBegin = gridThickTunnel + gridThickTrap + gridThickBlock; iyEnd = iyBegin + gridThickIso;
	setSingleElement(indexElement, currRegion, ixBegin, ixEnd, iyBegin, iyEnd);


}

double TripleCellsFull::getVertCoordX(int idX, int idY)
{
	double currdX = 0; // the coordinate in x-direction, in [nm], without normalization 

	static double lengthPerGridIso1 = widthIso1 / gridWidthIso1;
	static double lengthPerGridGate1 = widthGate1 / gridWidthGate1;
	static double lengthPerGridIso2 = widthIso2 / gridWidthIso2;
	static double lengthPerGridGate2 = widthGate2 / gridWidthGate2;
	static double lengthPerGridIso3 = widthIso3 / gridWidthIso3;
	static double lengthPerGridGate3 = widthGate3 / gridWidthGate3;
	static double lengthPerGridIso4 = widthIso4 / gridWidthIso4;

	if (!isValidVertex(idX, idY))
	{
		return -1;
	}

	// up to here, invalid point is wiped out.
	if (idX <= gridWidthIso1)
	{
		currdX = idX * lengthPerGridIso1;
		return currdX;
	}
	idX -= gridWidthIso1;
	currdX = widthIso1;
	if (idX <= gridWidthGate1)
	{
		currdX += idX * lengthPerGridGate1;
		return currdX;
	}
	idX -= gridWidthGate1;
	currdX += widthGate1;
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
	if (idX <= gridWidthGate3)
	{
		currdX += idX * lengthPerGridGate3;
		return currdX;
	}
	idX -= gridWidthGate3;
	currdX += widthGate3;

	// idX <= gridWidthIso4 is always true, if the checking runs to here.
	currdX += idX * lengthPerGridIso4;
	return currdX;
}

bool TripleCellsFull::isValidVertex(int idX, int idY)
{
	bool ret = true;
	static int gridTotalX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + 
		gridWidthIso3 + gridWidthGate3 + gridWidthIso4;
	static int gridMainY = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation
	static int gridTotalY = gridMainY + gridThickIso;
	if (idX < 0 || idX > gridTotalX)
	{
		ret = false;
	}
	if (idY < 0 || idY > gridTotalY)
	{
		ret = false;
	}
	if (idY > gridMainY)
	{
		//check gate1 region
		if (idX > gridWidthIso1 && idX < gridWidthIso1 + gridWidthGate1)
		{
			ret = false;
		}
		//check gate2 region
		if (idX > gridWidthIso1 + gridWidthGate1 + gridWidthIso2 && idX < gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2)
		{
			ret = false;
		}
		//check gate3 region
		if (idX > gridTotalX - gridWidthGate3 - gridWidthIso4 && idX < gridTotalX - gridWidthIso4)
		{
			ret = false;
		}
	}
	return ret;

}

void TripleCellsFull::setAdjacency()
{
	int gridTotalX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + 
					gridWidthIso3 + gridWidthGate3 + gridWidthIso4;
	int gridTotalY = gridThickTunnel + gridThickTrap + gridThickBlock + gridThickIso;

	FDVertex *currVert = NULL;
	FDVertex *adjacentVert = NULL;

	int vertID = 0;
	int adjacentID = 0;
	for (int iy = 0; iy <= gridTotalY; ++iy)
	{
		for (int ix = 0; ix <= gridTotalX; ++ix)
		{
			vertID = getVertIdAt(ix, iy);
			if (vertID < 0)//(ix, iy) does not corresponds to a valid vertex
			{
				continue;
			}
			currVert = this->GetVertex(vertID);
			//west vertex
			adjacentID = getVertIdAt(ix - 1, iy);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->WestVertex = adjacentVert;
				currVert->WestLength = FDVertex::Distance(currVert, currVert->WestVertex);
			}
			//north vertex
			adjacentID = getVertIdAt(ix, iy + 1);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->NorthVertex = adjacentVert;
				currVert->NorthLength = FDVertex::Distance(currVert, currVert->NorthVertex);
			}
			//east vertex
			adjacentID = getVertIdAt(ix + 1, iy);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->EastVertex = adjacentVert;
				currVert->EastLength = FDVertex::Distance(currVert, currVert->EastVertex);
			}
			//south vertex
			adjacentID = getVertIdAt(ix, iy - 1);
			if (adjacentID >= 0)
			{
				adjacentVert = this->GetVertex(adjacentID);
				currVert->SouthVertex = adjacentVert;
				currVert->SouthLength = FDVertex::Distance(currVert, currVert->SouthVertex);
			}
		}
	}
}

int TripleCellsFull::getVertIdAt(int idX, int idY)
{
	int vertID = 0;
	static int gridTotalX = gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + 
							gridWidthIso3 + gridWidthGate3 + gridWidthIso4;
	static int gridMainY = gridThickTunnel + gridThickTrap + gridThickBlock; // the top boundary index in y direction for the layers except isolation
	static int gridTotalY = gridMainY + gridThickIso;

	//out of boundary
	if (idX < 0 || idX > gridTotalX)
	{
		return -1;
	}
	if (idY < 0 || idY > gridTotalY)
	{
		return -1;
	}

	if (idY <= gridMainY)// the vertex is in the main part
	{
		vertID = (gridTotalX + 1)*idY + idX;
	}
	else// the vertex is in the isolation layer 
	{
		vertID = (gridTotalX + 1)*(gridMainY + 1);//count the main part
		vertID += (gridWidthIso1 + 1 + gridWidthIso2 + 1 + gridWidthIso3 + 1 + gridWidthIso4 + 1)*(idY - gridMainY - 1);//count the isolation vertex
		if (idX <= gridWidthIso1)//in isolation 1
		{
			vertID += idX;
		}
		else if (idX < gridWidthIso1 + gridWidthGate1)//in gate1
		{
			return -1;//-1 means invalid id
		}
		else if (idX <= gridWidthIso1 + gridWidthGate1 + gridWidthIso2)//in isolation2
		{	
			vertID += gridWidthIso1 + 1;
			vertID += idX - (gridWidthIso1 + gridWidthGate1);
		}
		else if (idX < gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2)//in gate2
		{
			return -1;
		}
		else if (idX <= gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3)//in isolation3
		{
			vertID += gridWidthIso1 + 1 + gridWidthIso2 + 1;
			vertID += idX - (gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2);
		}
		else if (idX < gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3)//in gate3
		{
			return -1;
		}
		else
		{
			vertID += gridWidthIso1 + 1 + gridWidthIso2 + 1 + gridWidthIso3 + 1;
			vertID += idX - (gridWidthIso1 + gridWidthGate1 + gridWidthIso2 + gridWidthGate2 + gridWidthIso3 + gridWidthGate3);
		}
	}
	return vertID;
}

TripleCellsFull::TripleCellsFull()
{

}

void TripleCellsFull::setSingleElement(int &idElem, FDRegion *region, int xbegin, int xend, int ybegin, int yend)
{
	FDElement *currElem = NULL;
	FDVertex *swVertex = NULL;
	FDVertex *seVertex = NULL;
	FDVertex *nwVertex = NULL;
	FDVertex *neVertex = NULL;

	for (int iyElem = ybegin; iyElem != yend; ++iyElem)
	{
		for (int ixElem = xbegin; ixElem != xend; ++ixElem)
		{
			swVertex = GetVertex(getVertIdAt(ixElem, iyElem));
			seVertex = GetVertex(getVertIdAt(ixElem + 1, iyElem));
			neVertex = GetVertex(getVertIdAt(ixElem + 1, iyElem + 1));
			nwVertex = GetVertex(getVertIdAt(ixElem, iyElem + 1));
			currElem = new FDElement(idElem++, swVertex, seVertex, neVertex, nwVertex);
			this->elements.push_back(currElem);
			currElem->SetVertexAdjacent();
			currElem->SetRegion(region);
			region->AddElement(currElem);
		}
	}
}

double TripleCellsFull::getVertCoordY(int idX, int idY)
{
	double currdY = 0; // the coordinate in y-direction, in [nm], without normalization
	static double lengthPerGridTunnel = thickTunnel / gridThickTunnel;
	static double lengthPerGridTrap = thickTrap / gridThickTrap;
	static double lengthPerGridBlock = thickBlock / gridThickBlock;
	static double lengthPerGridIso = thickIso / gridThickIso;

	if (!isValidVertex(idX, idY))
	{
		return -1;
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

bool TripleCellsFull::IsEndOfEffectiveCapacitor(FDVertex *vert)
{
	FDElement *elem = NULL;

	elem = vert->NorthwestElem;
	if (elem != NULL && (elem->Region->RegName == "Iso1" || elem->Region->RegName == "Iso2"
							|| elem->Region->RegName == "Iso3" || elem->Region->RegName == "Iso4"))
	{
		return true;
	}
	elem = vert->NortheastElem;
	if (elem != NULL && (elem->Region->RegName == "Iso1" || elem->Region->RegName == "Iso2"
							|| elem->Region->RegName == "Iso3" || elem->Region->RegName == "Iso4"))
	{
		return true;
	}
	elem = vert->SouthwestElem;
	if (elem != NULL && (elem->Region->RegName == "Iso1" || elem->Region->RegName == "Iso2"
							|| elem->Region->RegName == "Iso3" || elem->Region->RegName == "Iso4"))
	{
		return true;
	}
	elem = vert->SoutheastElem;
	if (elem != NULL && (elem->Region->RegName == "Iso1" || elem->Region->RegName == "Iso2"
							|| elem->Region->RegName == "Iso3" || elem->Region->RegName == "Iso4"))
	{
		return true;
	}
	return false;
}
