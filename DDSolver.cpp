/**
* @file DDSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-29   19:14
* @note
* @todo
*/

#include "DDSolver.h"
#include "FDDomain.h"
#include "Material.h"
#include "SctmPhys.h"
using namespace SctmPhys;

DriftDiffusionSolver::DriftDiffusionSolver(FDDomain *domain): totalVertices(domain->GetVertices())
{
	getDDVertices(domain);
	prepareSolver();
}

void DriftDiffusionSolver::prepareSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");
	
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->eDensity.resize(vertSize);

	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	buildVertexMap();
	buildCoefficientMatrix();
	buildRhsVector();
	//the structure does not change, so the final coefficient matrix after refreshing does not change either
	refreshCoefficientMatrix();
}

void DriftDiffusionSolver::buildVertexMap()
{
	std::pair<MapForVertex::iterator, bool> insertPairVertex;
	std::pair<MapForPrpty::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->vertMap.insert(MapForVertex::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(MapForPrpty::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(MapForPrpty::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastDensityMap.insert(MapForPrpty::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eDensity)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &current)
{

}

void DriftDiffusionSolver::buildCoefficientMatrix()
{
	int matrixSize = this->vertices.size();
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));

	FDVertex *currVert = NULL;
	int vertID = 0;
	int indexEquation = 0;
	int indexCoefficient = 0;
	
	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	double KT_div_q = SctmPhys::k0 * this->temperature / SctmPhys::q;

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = vertMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);

		coeff_center = 0;
		// in the boundary cases, the length of one or two direction may be zero 
		deltaX = (currVert->EastLength + currVert->WestLength) / 2;
		deltaY = (currVert->NorthLength + currVert->SouthLength) / 2;

		//the initial filling method can be used for boundary vertices and non-boundary vertices
		if ( currVert->EastVertex != NULL )
		{
			mobility = (mobilityMap[currVert->EastVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = vertMap[currVert->EastVertex->GetID()];
			coeff_adjacent = mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				+ KT_div_q );
			coeff_center += mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				- KT_div_q );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if ( currVert->WestVertex != NULL )
		{
			mobility = (mobilityMap[currVert->WestVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = vertMap[currVert->WestVertex->GetID()];
			coeff_adjacent = mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 
				+ KT_div_q );
			coeff_center += mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 
				- KT_div_q );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if ( currVert->NorthVertex != NULL )
		{
			mobility = (mobilityMap[currVert->NorthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = vertMap[currVert->NorthVertex->GetID()];
			coeff_adjacent = mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				+ KT_div_q );
			coeff_center += mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				- KT_div_q );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if ( currVert->SouthVertex != NULL )
		{
			mobility = (mobilityMap[currVert->SouthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = vertMap[currVert->SouthVertex->GetID()];
			coeff_adjacent = mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				+ KT_div_q );
			coeff_center += mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2
				- KT_div_q );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}


		coeff_center += -1 / timeStep; // from pn/pt, p=partial differential
 		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetInternalID]
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
	}
}

void DriftDiffusionSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	double rhsVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		rhsVal = -lastDensityMap[currVert->GetID()] / timeStep;
		this->rhsVector.at(iVert) = rhsVal;
	}
}

void DriftDiffusionSolver::refreshRhs()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double q = SctmPhys::q;
	double bcVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		if (!currVert->BndCond.Valid(FDBoundary::eCurrentDensity))
			continue;

		//the vertex is at current density boundary condition
		equationID = vertMap[currVert->GetID()]; //equationIndex = iVert
		if (currVert->BndCond.GetBCType(FDBoundary::eCurrentDensity) != FDBoundary::BC_Dirichlet)
			continue;

		//the vertex is at BC_Dirichlet boundary condition
		//only solve DD equation in trapping region
		if ( currVert->WestVertex == NULL
			|| (currVert->WestVertex->NorthwestElem == NULL ? false : currVert->WestVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			|| (currVert->WestVertex->SouthwestElem == NULL ? false : currVert->WestVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += 1 / q * bcVal / (currVert->EastLength / 2);
		}

		if ( currVert->EastVertex == NULL
			|| (currVert->EastVertex->NortheastElem == NULL ? false : currVert->EastVertex->NortheastElem->Region->Type != FDRegion::Trapping)
			|| (currVert->EastVertex->SoutheastElem == NULL ? false : currVert->EastVertex->SoutheastElem->Region->Type != FDRegion::Trapping) )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += -1 / q * bcVal / (currVert->WestLength / 2);
		}

		if ( currVert->SouthVertex == NULL
			|| (currVert->SouthVertex->SoutheastElem == NULL ? false : currVert->SouthVertex->SoutheastElem->Region->Type != FDRegion::Trapping)
			|| (currVert->SouthVertex->SouthwestElem == NULL ? false : currVert->SouthVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += 1 / q * bcVal / (currVert->NorthLength / 2);
		}

		if ( currVert->NorthVertex == NULL
			|| (currVert->NorthVertex->NorthwestElem == NULL ? false : currVert->NorthVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			|| (currVert->NorthVertex->NortheastElem == NULL ? false : currVert->NorthVertex->NortheastElem->Region->Type != FDRegion::Trapping) )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += -1 / q * bcVal / (currVert->SouthLength / 2);
		}
	}
}

void DriftDiffusionSolver::refreshCoefficientMatrix()
{
	//The boundary conditions are always BC_Dirichlet, so there is no need refreshing the matrix.
}

void DriftDiffusionSolver::getDDVertices(FDDomain *domain)
{
	FDVertex *currVert = NULL;

	for (size_t iVert = 0; iVert != this->totalVertices.size(); ++iVert)
	{
		currVert = this->totalVertices.at(iVert);
		if (
			((currVert->NorthwestElem != NULL) 
			&& (currVert->NorthwestElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->NortheastElem != NULL) 
			&& (currVert->NortheastElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->SoutheastElem != NULL) 
			&& (currVert->SoutheastElem->Region->Type == FDRegion::Trapping))
			
			||((currVert->SouthwestElem != NULL) 
			&& (currVert->SouthwestElem->Region->Type == FDRegion::Trapping))
		   )
		{
			this->vertices.push_back(currVert);
		}
	}
}

DDTest::DDTest(FDDomain *_domain) : DriftDiffusionSolver(_domain)
{

}
