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
using namespace SctmUtils;

DriftDiffusionSolver::DriftDiffusionSolver(FDDomain *domain): totalVertices(domain->GetVertices())
{
	getDDVertices(domain);
	//for testing the DD solver, the preparation is done in the DDTest class.
	//prepareSolver();
}

void DriftDiffusionSolver::SolveDD()
{
	
}

void DriftDiffusionSolver::initializeSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");
	
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);

	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	buildVertexMap();
	setTimeStep();
	buildCoefficientMatrix();
	buildRhsVector();
	//the structure does not change, so the final coefficient matrix after refreshing does not change either
	refreshCoefficientMatrix();
}

void DriftDiffusionSolver::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
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
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys.GetPhysPrpty(PhysProperty::eDensity)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &in_current)
{
}

void DriftDiffusionSolver::buildCoefficientMatrix()
{
	int matrixSize = this->vertices.size();
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));

	FDVertex *currVert = NULL;
	int indexEquation = 0;
	int indexCoefficient = 0;
	
	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	double KT_div_q = SctmPhys::k0 * this->temperature / SctmPhys::q;
	//we should consider the normalized current equation
	//for electron, J = -u[ n * p_phi/p_x - p_n/p_x ]

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);

		coeff_center = 0;
		// in the boundary cases, the length of one or two direction may be zero 
		deltaX = (currVert->EastLength + currVert->WestLength) / 2;
		deltaY = (currVert->NorthLength + currVert->SouthLength) / 2;
		SCTM_ASSERT(deltaX!=0 && deltaY!=0, 10015);

		//the initial filling method can be used for boundary vertices and non-boundary vertices
		if (   (currVert->EastVertex != NULL)
			 &&( (currVert->NortheastElem == NULL ? false : currVert->NortheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SoutheastElem == NULL ? false : currVert->SoutheastElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->EastVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->EastVertex->GetID()];
			coeff_adjacent = - mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->EastLength / deltaX *
				( (potentialMap[currVert->EastVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->WestVertex != NULL)
		     &&( (currVert->NorthwestElem == NULL ? false : currVert->NorthwestElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SouthwestElem == NULL ? false : currVert->SouthwestElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->WestVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->WestVertex->GetID()];
			coeff_adjacent = - mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 - 1 );
			coeff_center += - mobility / currVert->WestLength / deltaX *
				( (potentialMap[currVert->WestVertex->GetID()] - potentialMap[currVert->GetID()]) /2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->NorthVertex != NULL)
			 &&( (currVert->NortheastElem == NULL ? false : currVert->NortheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->NorthwestElem == NULL ? false : currVert->NorthwestElem->Region->Type == FDRegion::Trapping) )
		   )
		{
			mobility = (mobilityMap[currVert->NorthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->NorthVertex->GetID()];
			coeff_adjacent = - mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->NorthLength / deltaY *
				( (potentialMap[currVert->NorthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		if (   (currVert->SouthVertex != NULL)
			 &&( (currVert->SoutheastElem == NULL ? false : currVert->SoutheastElem->Region->Type == FDRegion::Trapping)
			   ||(currVert->SouthwestElem == NULL ? false : currVert->SouthwestElem->Region->Type == FDRegion::Trapping))
		   )
		{
			mobility = (mobilityMap[currVert->SouthVertex->GetID()] + mobilityMap[currVert->GetID()]) / 2;
			indexCoefficient = equationMap[currVert->SouthVertex->GetID()];
			coeff_adjacent = - mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / currVert->SouthLength / deltaY *
				( (potentialMap[currVert->SouthVertex->GetID()] - potentialMap[currVert->GetID()]) / 2 + 1 );

			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
		}

		//coeff_center += -1 / timeStep; // from p_n/p_t, p=partial differential

		indexCoefficient = equationMap[currVert->GetID()];
 		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetID]
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
	}
}

void DriftDiffusionSolver::buildCoefficientMatrix(bool newMethod)
{
	int indexEquation = 0;
	int matrixSize = this->vertices.size();
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));

	FDVertex *currVert = NULL;

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);
		//boundary vertex
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			setCoefficientBCVertex(currVert);
			continue;
		}
		//inner vertex
		setCoefficientInnerVertex(currVert);
	}
}

void DriftDiffusionSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	double rhsVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			switch (currVert->BndCond.GetBCType(FDBoundary::eDensity))
			{
			case FDBoundary::BC_Dirichlet:
				rhsVal = currVert->BndCond.GetBCValue(FDBoundary::eDensity);
				break;
			case FDBoundary::BC_Cauchy:
				//for Cauchy boundary condition
				rhsVal = currVert->BndCond.GetBCValue(FDBoundary::eDensity);
				break;	
			}
		}
		else
		{
			//the inner vertex
			rhsVal = -lastElecDensMap[currVert->GetID()] / timeStep;
		}
		this->rhsVector.at(iVert) = rhsVal;
	}
}

void DriftDiffusionSolver::refreshRhsWithBC()
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
		equationID = equationMap[currVert->GetID()]; //equationIndex = iVert
		SCTM_ASSERT(currVert->BndCond.GetBCType(FDBoundary::eCurrentDensity) == FDBoundary::BC_Dirichlet, 10013);
		//if (currVert->BndCond.GetBCType(FDBoundary::eCurrentDensity) != FDBoundary::BC_Dirichlet)
		//	continue;

		//the vertex is at BC_Dirichlet boundary condition
		//only solve DD equation in trapping region
		if ( currVert->WestVertex == NULL
			||(  (currVert->WestVertex->NorthwestElem == NULL ? true : currVert->WestVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->WestVertex->SouthwestElem == NULL ? true : currVert->WestVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += bcVal / (currVert->EastLength / 2);
		}

		if ( currVert->EastVertex == NULL
			||(  (currVert->EastVertex->NortheastElem == NULL ? true : currVert->EastVertex->NortheastElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->EastVertex->SoutheastElem == NULL ? true : currVert->EastVertex->SoutheastElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueWestEast(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += - bcVal / (currVert->WestLength / 2);
		}

		if ( currVert->SouthVertex == NULL
			||(  (currVert->SouthVertex->SoutheastElem == NULL ? true : currVert->SouthVertex->SoutheastElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->SouthVertex->SouthwestElem == NULL ? true : currVert->SouthVertex->SouthwestElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += bcVal / (currVert->NorthLength / 2);
		}

		if ( currVert->NorthVertex == NULL
			||(  (currVert->NorthVertex->NorthwestElem == NULL ? true : currVert->NorthVertex->NorthwestElem->Region->Type != FDRegion::Trapping)
			   &&(currVert->NorthVertex->NortheastElem == NULL ? true : currVert->NorthVertex->NortheastElem->Region->Type != FDRegion::Trapping) )
		   )
		{
			bcVal = currVert->BndCond.GetBCValueSouthNorth(FDBoundary::eCurrentDensity);
			rhsVector.at(equationID) += - bcVal / (currVert->SouthLength / 2);
		}
	}
}

void DriftDiffusionSolver::refreshCoefficientMatrix()
{
	//The boundary conditions are always BC_Cauchy, so there is no need to refresh the matrix for this reason.
	//However, due to that the time step is likely to be different in each simulation step, the time dependent item
	//in the coefficient matrix has to be add to the corresponding item.
	FDVertex *currVert = NULL;
	int indexEquation = 0;
	int indexCoefficient = 0;

	double coeffToAdd = 0; // coefficient for center vertex
	coeffToAdd = -1 / timeStep; // from p_n/p_t, p=partial differential
	
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);

		//there is no need refreshing coefficient related to boundary condition for this reason
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
			continue;

		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);

		indexCoefficient = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetID]
		//=iVert (at present)
		
		matrixSolver.RefreshMatrixValue(indexEquation, indexCoefficient, coeffToAdd, SctmSparseMatrixSolver::Add);
	}
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

void DriftDiffusionSolver::setTimeStep()
{
	timeStep = UtilsTimeStep.NextTimeStep();
}

void DriftDiffusionSolver::fillBackElecDens()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double edens = 0;
	int VertID = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		VertID = currVert->GetID();
		equationID = equationMap[VertID];
		edens = this->elecDensity.at(equationID);
		currVert->Phys.SetPhysPrpty(PhysProperty::eDensity, edens);
		//It is also essential to refresh property map
		lastElecDensMap[equationID] = edens;
	}
}

void DriftDiffusionSolver::setCoefficientBCVertex(FDVertex *vert)
{
	//for Cauchy boundary condition
	int indexEquation = 0;
	int indexCoeff = 0;

	double coeff_center = 0;
	double coeff_east = 0;
	double coeff_west = 0;
	double coeff_south = 0;
	double coeff_north = 0;

	double val_alpha = 0;
	double val_beta = 0;

	double mobility = mobilityMap[vert->GetID()]; 
	//the vector of normal direction is (alpha, beta)
	double norm_alpha = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).X();
	double norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).Y();

	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
	case FDBoundary::BC_Dirichlet:
		coeff_center = 1;
		indexEquation = equationMap[vert->GetID()];
		indexCoeff = equationMap[vert->GetID()];
		matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_center;
		break;
	case FDBoundary::BC_Cauchy:
		//for center vertex
		if (norm_alpha > 0)
		{
			coeff_center = -norm_alpha * ( potentialMap[vert->GetID()] - potentialMap[vert->WestVertex->GetID()] ) / vert->WestLength;
			coeff_center += norm_alpha * 1 / vert->WestLength;
			coeff_west = norm_alpha * (-1) / vert->WestLength;
		}
		else if (norm_alpha < 0)
		{
			coeff_center = -norm_alpha * ( potentialMap[vert->GetID()] - potentialMap[vert->EastVertex->GetID()] ) / (-vert->EastLength);
			coeff_center += norm_alpha * 1 / (-vert->EastLength);
			coeff_east = norm_alpha * (-1) / (-vert->EastLength);
		}
		else
			coeff_center = 0;

		if (norm_beta > 0)
		{
			coeff_center += -norm_beta * ( potentialMap[vert->GetID()] - potentialMap[vert->SouthVertex->GetID()] ) / vert->SouthLength;
			coeff_center += norm_beta * 1 / vert->SouthLength;
			coeff_south = norm_beta * (-1) / vert->SouthLength;
		}
		else if (norm_beta < 0)
		{
			coeff_center += -norm_beta * ( potentialMap[vert->GetID()] -  potentialMap[vert->NorthVertex->GetID()] ) / (-vert->NorthLength);
			coeff_center += norm_beta * 1 / (-vert->NorthLength);
			coeff_north = norm_beta * (-1) / (-vert->NorthLength);
		}
		else
			coeff_center += 0;

		coeff_center *= mobility;
		coeff_east *= mobility;
		coeff_west *= mobility;
		coeff_north *= mobility;
		coeff_south *= mobility;

		indexEquation = equationMap[vert->GetID()];
		indexCoeff = equationMap[vert->GetID()];
		matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_center;
		if (coeff_east != 0)
		{
			indexCoeff = equationMap[vert->EastVertex->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_east;
		}
		if (coeff_west != 0)
		{
			indexCoeff = equationMap[vert->WestVertex->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_west;
		}
		if (coeff_north != 0)
		{
			indexCoeff = equationMap[vert->NorthVertex->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_north;
		}
		if (coeff_south != 0)
		{
			indexCoeff = equationMap[vert->SouthVertex->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_south;
		}
		break;
	}

}

void DriftDiffusionSolver::setCoefficientInnerVertex(FDVertex *vert)
{
	int indexEquation = 0;
	int indexCoefficient = 0;

	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	//we should consider the normalized current equation
	//for electron, J = -u[ n * p_phi/p_x - p_n/p_x ]

	indexEquation = equationMap[vert->GetID()];
	
	coeff_center = 0;
	// in the boundary cases, the length of one or two direction may be zero 
	deltaX = (vert->EastLength + vert->WestLength) / 2;
	deltaY = (vert->NorthLength + vert->SouthLength) / 2;
	SCTM_ASSERT(deltaX!=0 && deltaY!=0, 10015);

	//the initial filling method can be used for boundary vertices and non-boundary vertices
	//related to east vertex
	mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->EastVertex->GetID()];
	coeff_adjacent = - mobility / vert->EastLength / deltaX *
		( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->EastLength / deltaX *
		( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;


	//related to west vertex
	mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->WestVertex->GetID()];
	coeff_adjacent = - mobility / vert->WestLength / deltaX *
		( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 - 1 );
	coeff_center += - mobility / vert->WestLength / deltaX *
		( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 + 1 );

	matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;


	//related to north vertex
	mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->NorthVertex->GetID()];
	coeff_adjacent = - mobility / vert->NorthLength / deltaY *
		( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->NorthLength / deltaY *
		( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;

	//related to south vertex
	mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->SouthVertex->GetID()];
	coeff_adjacent = - mobility / vert->SouthLength / deltaY *
		( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->SouthLength / deltaY *
		( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;

	//coeff_center += -1 / timeStep; // from p_n/p_t, p=partial differential

	//related to current center vertex
	indexCoefficient = equationMap[vert->GetID()];
	SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
	//indexCoefficent = indexEquation = vertMap[currVert->GetID]
	matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
}

DDTest::DDTest(FDDomain *_domain) : DriftDiffusionSolver(_domain)
{
	initializeSolver();
}

void DDTest::initializeSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");

	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);

	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	buildVertexMap();
	
	buildCoefficientMatrix(true);
}

void DDTest::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//TODO: this is only a temporary method
	//the following is used to set the mobility to uniform value, because the original setting method leads to incorrect calculation
	//of the mobility at the trapping layer interface
	Normalization norm = Normalization();
	double mobility = MaterialDB::GetMatPrpty(&MaterialDB::Si3N4, MaterialDB::MatProperty::Mat_ElectronMobility);

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, mobility));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, 0)); // suppose no electronic field exists
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, 0)); // suppose the initial carrier density is 0
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DDTest::setBndCondCurrent()
{
	FDVertex *currVert = NULL;
	int equaID = 0;
	double bcVal_in = 1e-3; //the magnitude of current density vector
	double bcVal_out = 0;
	// current density at boundary for test, in [A/cm^2]
	// note that the current direction is the reversed direction of electron flow

	Normalization norm = Normalization();
	bcVal_in = norm.PushCurrDens(bcVal_in);

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = isNotElemOf(FDRegion::Trapping, currVert->NorthwestElem);
		bool notTrapping_NE = isNotElemOf(FDRegion::Trapping, currVert->NortheastElem);
		bool notTrapping_SE = isNotElemOf(FDRegion::Trapping, currVert->SoutheastElem);
		bool notTrapping_SW = isNotElemOf(FDRegion::Trapping, currVert->SouthwestElem);

		bool valid_NW = isValidElem(currVert->NorthwestElem);
		bool valid_NE = isValidElem(currVert->NortheastElem);
		bool valid_SE = isValidElem(currVert->SoutheastElem);
		bool valid_SW = isValidElem(currVert->SouthwestElem);

		bool inWest = true;
		bool inEast = true;
		bool inSouth = true;
		bool inNorth = true;
		bool inNorthWest = true;
		bool inNorthEast = true;
		bool inSouthEast = true;
		bool inSouthWest = true;

		//Southeast corner
		if (( !notTrapping_NW && notTrapping_NE && 
			notTrapping_SW && notTrapping_SE ) && inSouthEast)
		{
			if (			!valid_NE &&
				valid_SW)
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if (			valid_NE &&
				!valid_SW)
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Southwest corner
		if (( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE) && inSouthWest)
		{
			if ( valid_NW &&
				!valid_SE )
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if ( !valid_NW &&
				valid_SE )
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Northwest corner
		if (( notTrapping_NW && notTrapping_NE &&
			notTrapping_SW && !notTrapping_SE ) && inNorthWest)
		{
			if (              valid_NE &&
				!valid_SW )
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if (              !valid_NE &&
				valid_SW )
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Northeast corner
		if (( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE ) && inNorthEast)
		{
			if ( valid_NW &&
				!valid_SE)
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if ( !valid_NW &&
				valid_SE)
			{
				currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//North side
		if (( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && !notTrapping_SE ) && inNorth)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//South side
		if (( !notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE ) && inSouth)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//East side
		if (( !notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE ) && inEast)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//West side
		if (( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && !notTrapping_SE) && inWest)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, bcVal_in);
			continue;
		}

	}
}

void DDTest::SolveDD()
{
	prepareSolver();

	this->matrixSolver.SolveMatrix(rhsVector, this->elecDensity);
	
	//UtilsDebug.PrintVector(this->rhsVector, "right hand side vector");
	UtilsDebug.PrintVector(this->elecDensity, "electron density");
	
	fillBackElecDens();

	SctmFileOperator write = SctmFileOperator("E:\\PhD Study\\SimCTM\\SctmTest\\DDTest\\eDensity.txt", SctmFileOperator::Write);
	write.WriteDDResult(this->vertices, "electron density");
}

void DDTest::prepareSolver()
{
	setTimeStep();
	setBndCondCurrent();

	//UtilsDebug.PrintSparseMatrixRow(matrixSolver.matrix, 0);
	UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);
	refreshCoefficientMatrix();
	//UtilsDebug.PrintSparseMatrixRow(matrixSolver.matrix, 0);
	//UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);

	//buildRhsVector and refreshRhsWithBC are called together, because for each simulation step, the initial building of Rhs is
	//different due to the difference in last time electron density
	buildRhsVector();
	UtilsDebug.PrintVector(this->rhsVector, "right hand side vector");
	//refreshRhsWithBC();
}

bool DDTest::isValidElem(FDElement *elem)
{
	return elem != NULL;
}

bool DDTest::isNotElemOf(FDRegion::RegionType rType, FDElement *elem)
{
	if (elem == NULL)
		return true;
	else
	{
		return elem->Region->Type != rType;
	}
}

void DDTest::setBndCondDensity()
{
	FDVertex *currVert = NULL;
	int equaID = 0;
	double bcVal_in = 1e12; //electron density in [cm^-3]
	double bcVal_out = 1e11;
	// current density at boundary for test, in [A/cm^2]
	// note that the current direction is the reversed direction of electron flow

	Normalization norm = Normalization();
	bcVal_in = norm.PushDensity(bcVal_in);
	bcVal_out = norm.PushDensity(bcVal_out);

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = isNotElemOf(FDRegion::Trapping, currVert->NorthwestElem);
		bool notTrapping_NE = isNotElemOf(FDRegion::Trapping, currVert->NortheastElem);
		bool notTrapping_SE = isNotElemOf(FDRegion::Trapping, currVert->SoutheastElem);
		bool notTrapping_SW = isNotElemOf(FDRegion::Trapping, currVert->SouthwestElem);

		bool valid_NW = isValidElem(currVert->NorthwestElem);
		bool valid_NE = isValidElem(currVert->NortheastElem);
		bool valid_SE = isValidElem(currVert->SoutheastElem);
		bool valid_SW = isValidElem(currVert->SouthwestElem);

		//Southeast corner
		if ( !notTrapping_NW && notTrapping_NE && 
			notTrapping_SW && notTrapping_SE )
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_in);
			continue;
		}

		//Southwest corner
		if ( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_in);
			continue;
		}

		//South side
		if ( !notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_in);
			continue;
		}

		//Northwest corner
		if ( notTrapping_NW && notTrapping_NE &&
			notTrapping_SW && !notTrapping_SE )
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_out);
			continue;
		}

		//Northeast corner
		if ( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE )
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_out);
			continue;
		}

		//North side
		if ( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && !notTrapping_SE)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, bcVal_out);
			continue;
		}

		//for testing
		/*
		//East side
		if ( !notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE)
		{
			vert->BndCond.SetBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Cauchy, 0, VectorValue(1, 0));
			return;
		}
		
		//West side
		if ( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && !notTrapping_SE)
		{
			currVert->BndCond.RefreshBndCondValue(true, FDBoundary::eDensity, bcVal_in);
			return;
		}
		*/
	}
}
