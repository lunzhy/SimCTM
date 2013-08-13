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
	stuffPotential();
	refreshBandEnergy();
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
	double xLength_in_nm = 90;
	double yLengthTunnel_in_nm = 4;
	double yLengthTrap_in_nm = 7;
	double yLengthBlock_in_nm = 9;
	int xGridNumber = 10;
	int yGridNumberTunnel = 15;
	int yGridNumberTrap = 15;
	int yGridNumberBlock = 15;
	////////////////////////////////////////////////////////////////////////////
	//here, the length of the parameter is conversed to [cm]
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

	yCntTotalVertex = yCntVertexTunnel + yCntVertexTrap + yCntVertexBlock - 1 - 1;
}

void SimpleONO::printStructure()
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
		std::cout << "id=" << elements.at(ix)->GetInternalID() << '\t' << elements.at(ix)->Region->Type
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

void SimpleONO::setDomainDetails()
{
	int cntVertex = 0;
	int cntElement = 0;
	int cntRegion = 0;

	double currCoordX = 0.0;
	double currCoordY = 0.0;

	double nm_in_cm = SctmPhys::nm_in_cm;

	////////////////////////////////////////////////////////////////////
	//set vertices
	Normalization theNorm = Normalization();
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntTotalVertex);
	FDDomainHelper elementHelper = FDDomainHelper(xCntVertex-1, yCntTotalVertex-1); // element number = vertex number - 1
	double normCoordX = 0.0;
	double normCoordY = 0.0;
	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		//currCoordY = iy * yGridTunnel;
		//currCoordY = theNorm.PushLength(currCoordY);
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			//currCoordX = ix * xGrid;
			//the coordinates are normalized before pushing into vertex
			normCoordX = theNorm.PushLength(currCoordX * nm_in_cm);
			normCoordY = theNorm.PushLength(currCoordY * nm_in_cm);
			vertices.push_back(new FDVertex(cntVertex, normCoordX, normCoordY));
			cntVertex++;
			currCoordX += xNextGridLength(ix);//non-normalized value, in cm
		}
		currCoordX = 0;
		currCoordY += yNextGridLength(iy);//non-normalized value, in cm
	}

	/////////////////////////////////////////////////////////////////////
	//set regions
	regions.push_back(new FDRegion(cntRegion, FDRegion::Tunneling));
	regions.back()->RegionMaterial = &MaterialDB::SiO2;
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::Trapping));
	regions.back()->RegionMaterial = &MaterialDB::Si3N4;
	cntRegion++;
	regions.push_back(new FDRegion(cntRegion, FDRegion::Blocking));
	regions.back()->RegionMaterial = &MaterialDB::SiO2;
	cntRegion++;

	/////////////////////////////////////////////////////////////////////
	//set elements with specified region
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

void SimpleONO::setAdjacency()
{
	//set vertex properties
	FDDomainHelper vertexHelper = FDDomainHelper(xCntVertex, yCntTotalVertex);
	FDDomainHelper elementHelper = FDDomainHelper(xCntVertex-1, yCntTotalVertex-1);
	int id = 0;
	FDVertex *currVertex = NULL;
	for (int iy = 0; iy != yCntTotalVertex; ++iy)
	{
		for (int ix = 0; ix != xCntVertex; ++ix)
		{
			id = vertexHelper.IdAt(ix, iy);
			currVertex = getVertex(id);
			////////////////////////////////
			//set adjacent vertex and length
			//if current vertex doesn't have a west/east/south/north edge, the length is set 0.
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
			//////////////////////
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

	/////////////////////////////////////////////////////////////////////////////////////////////////
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

void SimpleONO::stuffPotential()
{
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
			currVertex = getVertex(id);
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
	double tot = 0; // total area
	double sum = 0; // sum corresponds to integral value
	double physValue = 0;

	using namespace MaterialDB;
	using namespace SctmPhys;
	using std::vector;

	vector<MatProperty::Name> matPrptys;
	vector<PhysProperty::Name> vertexPhyPrpty;
	matPrptys.push_back(MatProperty::Mat_ElectronAffinity); vertexPhyPrpty.push_back(PhysProperty::ElectronAffinity);
	matPrptys.push_back(MatProperty::Mat_ElectronMass); vertexPhyPrpty.push_back(PhysProperty::ElectronMass);
	matPrptys.push_back(MatProperty::Mat_Bandgap); vertexPhyPrpty.push_back(PhysProperty::Bandgap);

	//iteration over the vertices
	for (std::size_t iVer = 0; iVer != this->vertices.size(); ++iVer)
	{
		currVertex = getVertex(iVer);
		//iteration over the physical properties to be set from material property
		for (std::size_t iPrpty = 0; iPrpty != matPrptys.size(); ++iPrpty)
		{
			tot = 0; sum = 0;
			currElem = currVertex->SouthwestElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->SoutheastElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->NortheastElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrptys.at(iPrpty)) * currElem->Area : 0;
			
			currElem = currVertex->NorthwestElem;
			tot += ( currElem != NULL ) ? currElem->Area : 0;
			sum += ( currElem != NULL ) ? GetMatPrpty(currElem->Region->RegionMaterial, matPrptys.at(iPrpty)) * currElem->Area : 0;

			physValue = sum / tot;
			currVertex->Phys.SetPhysPrpty(vertexPhyPrpty.at(iPrpty), physValue);
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
		currVertex = getVertex(iver);

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
