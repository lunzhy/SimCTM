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
#include "SimpleONO.h"
#include "Normalization.h"
#include "SctmUtils.h"
using namespace SctmPhys;
using namespace SctmUtils;

DriftDiffusionSolver::DriftDiffusionSolver(FDDomain *domain): totalVertices(domain->GetVertices())
{
	this->bcMethod = UsingCurrentDensity;
	this->useCrankNicolsonMethod = true;
	this->lastTimeStep = 0;
	getDDVertices(domain);
	initializeSolver();
}

void DriftDiffusionSolver::SolveDD()
{
	
}

void DriftDiffusionSolver::initializeSolver()
{
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);
	this->temperature = 300; //temporarily used here TODO: modify to accord with the whole simulation
	
	buildVertexMap(); //call the method in DDSolver (Base class), because at this time the derived class is not constructed.
	buildCoefficientMatrix();
}

void DriftDiffusionSolver::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//Previously another method of building vertex map is used in DDTest due to the incorrect calculation of
	//mobility of boundary vertex. This is now solved by introducing a new method of filling vertex physical values.
	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		//equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eDensity)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		equationID += 1;
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &in_current)
{
}

void DriftDiffusionSolver::buildCoefficientMatrix_Obsolete()
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

void DriftDiffusionSolver::buildCoefficientMatrix()
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
			switch (this->bcMethod)
			{
			case DirectDiscretization:
				setCoefficientBCVertex_DirectDiscretization(currVert);
				break;
			case UsingCurrentDensity:
				setCoefficientBCVertex_UsingCurrent(currVert);
				break;
			}
			continue;
		}
		//inner vertex
		setCoefficientInnerVertex(currVert);
	}
}

void DriftDiffusionSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	int equationID = 0;

	double rhsVal = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();
		equationID = equationMap[vertID];

		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			//for boundary vertex
			switch (this->bcMethod)
			{
			case DirectDiscretization:
				//the following function will judge the type of the boundary condition.
				rhsVal = getRhsBCVertex_DirectDiscretiztion(currVert);
				break;
			case UsingCurrentDensity:
				rhsVal = getRhsBCVertex_UsingCurrent(currVert);
				break;
			}
		}
		else
		{
			//for inner vertex
			rhsVal = getRhsInnerVertex(currVert);
		}
		this->rhsVector.at(equationID) = rhsVal;
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
	if ( lastTimeStep == 0 )
		coeffToAdd = -1 / timeStep; // from p_n/p_t, p=partial differential
	else
		coeffToAdd = -1 / timeStep - ( -1 / lastTimeStep );
	
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);

		//there is no need refreshing coefficient related to boundary condition for this reason
		//this is true for two method of differentiating boundary equation
		//if (currVert->IsAtBoundary(FDBoundary::eDensity))
		//	continue;

		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 10012);

		indexCoefficient = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetID]
		//=iVert (at present)
		
		matrixSolver.RefreshMatrixValue(indexEquation, indexCoefficient, coeffToAdd, SctmSparseMatrixSolver::Add);
	}
	lastTimeStep = timeStep; // to revert the matrix coefficient, this position is important
}

void DriftDiffusionSolver::getDDVertices(FDDomain *domain)
{
	FDVertex *currVert = NULL;

	bool notTrapping_NW = false;
	bool notTrapping_NE = false;
	bool notTrapping_SE = false;
	bool notTrapping_SW = false;

	for (size_t iVert = 0; iVert != this->totalVertices.size(); ++iVert)
	{
		currVert = this->totalVertices.at(iVert);

		bool notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NorthwestElem);
		bool notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NortheastElem);
		bool notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SoutheastElem);
		bool notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SouthwestElem);

		if (!(notTrapping_NE && notTrapping_NW && notTrapping_SE && notTrapping_SW))
		{
			this->vertices.push_back(currVert);
		}
	}
}

void DriftDiffusionSolver::setTimeStep()
{
 	timeStep = UtilsTimeStep.TimeStep();
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
		currVert->Phys->SetPhysPrpty(PhysProperty::eDensity, edens);
		//It is also essential to refresh property map
		lastElecDensMap[VertID] = edens;
	}
}

void DriftDiffusionSolver::setCoefficientBCVertex_DirectDiscretization(FDVertex *vert)
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

	//"coeff / 2" for Crank-Nicolson discretization method
	//related to east vertex
	mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->EastVertex->GetID()];
	coeff_adjacent = - mobility / vert->EastLength / deltaX *
		( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->EastLength / deltaX *
		( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	if (useCrankNicolsonMethod)
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
	}
	else
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
	}

	//related to west vertex
	mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->WestVertex->GetID()];
	coeff_adjacent = - mobility / vert->WestLength / deltaX *
		( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 - 1 );
	coeff_center += - mobility / vert->WestLength / deltaX *
		( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 + 1 );

	if (useCrankNicolsonMethod)
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
	}
	else
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
	}

	//related to north vertex
	mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->NorthVertex->GetID()];
	coeff_adjacent = - mobility / vert->NorthLength / deltaY *
		( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->NorthLength / deltaY *
		( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	if (useCrankNicolsonMethod)
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
	}
	else
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
	}

	//related to south vertex
	mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
	indexCoefficient = equationMap[vert->SouthVertex->GetID()];
	coeff_adjacent = - mobility / vert->SouthLength / deltaY *
		( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
	coeff_center += - mobility / vert->SouthLength / deltaY *
		( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

	if (useCrankNicolsonMethod)
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
	}
	else
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
	}

	//time-related modification of coefficient is done in refreshing coefficient
	//coeff_center += -1 / timeStep; // from p_n/p_t, p=partial differential

	//related to current center vertex
	indexCoefficient = equationMap[vert->GetID()];
	SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
	//indexCoefficent = indexEquation = vertMap[currVert->GetID]
	if (useCrankNicolsonMethod)
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center / 2;
	}
	else
	{
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
	}
	
}

void DriftDiffusionSolver::prepareSolver()
{
	setTimeStep();
	refreshBoundary();

	//UtilsDebug.PrintSparseMatrixRow(matrixSolver.matrix, 0);
	//UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);
	refreshCoefficientMatrix();
	//UtilsDebug.PrintSparseMatrixRow(matrixSolver.matrix, 0);
	UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);

	//buildRhsVector and refreshRhsWithBC are called together, because for each simulation step, the initial building of Rhs is
	//different due to the difference in last time electron density
	buildRhsVector();
	UtilsDebug.PrintVector(this->rhsVector, "right hand side vector");
	//refreshRhsWithBC();
}

void DriftDiffusionSolver::refreshBoundary()
{
	//current nothing is done here.
}

void DriftDiffusionSolver::setCoefficientBCVertex_UsingCurrent(FDVertex *vert)
{
	int indexEquation = 0;
	int indexCoefficient = 0;

	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex, could be west, east, south, north
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	bool notTrapping_NW = false;
	bool notTrapping_NE = false;
	bool notTrapping_SE = false;
	bool notTrapping_SW = false;

	//Notice. can not set the row of matrix according to the boundary normal direction, this is because the boundary condition
	//at the corner has been specially treated, which can not be applied to setting the matrix.
	//double norm_alpha =  vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormX();
	//double norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormY();

	//we should consider the normalized current equation
	//for electron, J = -u[ n * p_phi/p_x - p_n/p_x ]

	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
	case FDBoundary::BC_Dirichlet:
		coeff_center = 1;
		indexEquation = equationMap[vert->GetID()];
		indexCoefficient = equationMap[vert->GetID()];
		matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
		break;

	case FDBoundary::BC_Cauchy:
		indexEquation = equationMap[vert->GetID()];

		notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NorthwestElem);
		notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NortheastElem);
		notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SoutheastElem);
		notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SouthwestElem);

		getDeltaXYAtVertex(vert, deltaX, deltaY);

		//West vertex
		if (!(notTrapping_NW && notTrapping_SW))
		{
			mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			indexCoefficient = equationMap[vert->WestVertex->GetID()];
			coeff_adjacent = - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

			if (useCrankNicolsonMethod)
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
			}
			else
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
			}
			
		}

		//East vertex
		if (!(notTrapping_NE && notTrapping_SE))
		{
			mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			indexCoefficient = equationMap[vert->EastVertex->GetID()];
			coeff_adjacent = - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

			if (useCrankNicolsonMethod)
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
			}
			else
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
			}
		}

		//South vertex
		if (!(notTrapping_SE && notTrapping_SW))
		{
			mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			indexCoefficient = equationMap[vert->SouthVertex->GetID()];
			coeff_adjacent = - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

			if (useCrankNicolsonMethod)
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
			}
			else
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
			}
		}

		//North vertex
		if (!(notTrapping_NE && notTrapping_NW))
		{
			mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			indexCoefficient = equationMap[vert->NorthVertex->GetID()];
			coeff_adjacent = - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
			coeff_center += - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );

			if (useCrankNicolsonMethod)
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
			}
			else
			{
				matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
			}
		}

		//coeff_center += -1 / timeStep; // from p_n/p_t, p=partial differential

		indexCoefficient = equationMap[vert->GetID()];
		SCTM_ASSERT(indexCoefficient==indexEquation, 10012);
		//indexCoefficent = indexEquation = vertMap[currVert->GetID]
		if (useCrankNicolsonMethod)
		{
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center / 2;
		}
		else
		{
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
		}
		break;
	}
}

double DriftDiffusionSolver::getRhsBCVertex_DirectDiscretiztion(FDVertex *vert)
{
	double retVal = 0;
	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
	case FDBoundary::BC_Dirichlet:
		retVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
		break;
	case FDBoundary::BC_Cauchy:
		//for Cauchy boundary condition
		retVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
		break;	
	}
	
	return retVal;
}

double DriftDiffusionSolver::getRhsInnerVertex(FDVertex *vert)
{
	double rhsTime = 0; // the addend in rhs related to time step
	double rhsLastStepCurrent = 0; // the addend related to current of last time step, in Crank-Nicolson method

	double mobility = 0;
	double deltaX = 0;
	double deltaY = 0;
	
	double retVal = 0;

	//related to time step
	rhsTime = -1 * lastElecDensMap[vert->GetID()] / timeStep;
 
	getDeltaXYAtVertex(vert, deltaX, deltaY);

	if (useCrankNicolsonMethod)
	{
		//for Crank-Nicolson discretization method
		//related to east vertex
		mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
		double east = - mobility / vert->EastLength / deltaX *
			( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
			lastElecDensMap[vert->GetID()];

		east += - mobility / vert->EastLength / deltaX *
			( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) * 
			lastElecDensMap[vert->EastVertex->GetID()];

		//related to west vertex
		mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
		double west = - mobility / vert->WestLength / deltaX *
			( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
			lastElecDensMap[vert->GetID()];

		west += - mobility / vert->WestLength / deltaX *
			( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
			lastElecDensMap[vert->WestVertex->GetID()];

		//related to north vertex
		mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
		double north = - mobility / vert->NorthLength / deltaY *
			( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
			lastElecDensMap[vert->GetID()];

		north += - mobility / vert->NorthLength / deltaY *
			( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
			lastElecDensMap[vert->NorthVertex->GetID()];

		//related to south vertex
		mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
		double south = - mobility / vert->SouthLength / deltaY *
			( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
			lastElecDensMap[vert->GetID()];

		south += - mobility / vert->SouthLength / deltaY *
			( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
			lastElecDensMap[vert->SouthVertex->GetID()];
		
		//use Crank-Nicolson method
		rhsLastStepCurrent = ( east + west + south + north ) / 2;

	}
	
	retVal = rhsTime - rhsLastStepCurrent;
	return retVal;
}

double DriftDiffusionSolver::getRhsBCVertex_UsingCurrent(FDVertex *vert)
{
	double retVal = 0;

	double rhsTime = 0; // rhs addend related to time and electron density of this vertex of last simulation step
	double rhsBoundary = 0; // rhs addend related to the current boundary condition
	double rhsLastStepCurrent = 0; // rhs addend ralated to the current density from other vertices in last simulation step
	double rhsLastBoundary = 0; // rhs addend related to boundary condition in last simulation step

	double deltaX = 0;
	double deltaY = 0;

	bool notTrapping_NW = false;
	bool notTrapping_NE = false;
	bool notTrapping_SE = false;
	bool notTrapping_SW = false;

	double mobility = 0;
	double east = 0;
	double west = 0;
	double north = 0;
	double south = 0;

	static VertexMapDouble lastBCMap = VertexMapDouble(); // to store the overall current density of boundary conditions

	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
	case FDBoundary::BC_Dirichlet:
		//in BC_Dirichlet, the the rhs of boundary vertex is directly set.
		retVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
		break;
	case FDBoundary::BC_Cauchy:
		//for Cauchy boundary condition
		//calculation of the addend related current simulation time step
		rhsTime= -1 * lastElecDensMap[vert->GetID()] / timeStep;

		notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NorthwestElem);
		notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NortheastElem);
		notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SoutheastElem);
		notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SouthwestElem);

		getDeltaXYAtVertex(vert, deltaX, deltaY);

		//West
		if (!(notTrapping_NW && notTrapping_SW))
		{
			mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			west += - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 - 1 ) *
				lastElecDensMap[vert->WestVertex->GetID()];

			west += - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 + 1 ) *
				lastElecDensMap[vert->GetID()];
		}
		//East
		if (!(notTrapping_NE && notTrapping_SE))
		{
			mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			east += - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->EastVertex->GetID()];

			east += - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];
		}
		//South
		if (!(notTrapping_SE && notTrapping_SW))
		{
			mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			south += - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->SouthVertex->GetID()];

			south += - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];
		}
		//North
		if (!(notTrapping_NE && notTrapping_NW))
		{
			mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			north += - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->NorthVertex->GetID()];

			north += - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];
		}

		double bcVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
		double norm_alpha =  vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormX();
		double norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormY();
		//p_J / p_x = (Je - Jw) / dx + (Jn - Js) / dy
		if (norm_alpha > 0) //lack of east part
		{
			rhsBoundary += bcVal * norm_alpha / deltaX;
		}
		else if (norm_alpha < 0) //lack of west part
		{
			rhsBoundary += - bcVal * norm_alpha / deltaX;
		}
		if (norm_beta > 0) // lack of north part
		{
			rhsBoundary += bcVal * norm_beta / deltaY;
		}
		else if (norm_beta < 0) //  lack of south part
		{
			rhsBoundary += - bcVal * norm_beta / deltaY;
		}
		
		if (this->useCrankNicolsonMethod)
		{
			rhsBoundary = rhsBoundary / 2;
			//calculation of the addend related to the current density of last simulation time step
			rhsLastStepCurrent = (west + east + south + north) / 2;
			//this addend is not included in the first simulation step
			if (lastBCMap.find(vert->GetID()) != lastBCMap.end())
			{
				rhsLastBoundary = lastBCMap[vert->GetID()] / 2;
			}
			//to store the overall current density of boundary conditions
			lastBCMap[vert->GetID()] = rhsBoundary;
		}

		//the sign represents moving the symbol from right to left
		retVal = rhsTime - rhsBoundary - rhsLastStepCurrent - rhsLastBoundary;
		break;
	}
	return retVal;
}

double DriftDiffusionSolver::CalculateTotalLineDensity()
{
	double ret = 0;
	
	bool notTrapping_NW = false;
	bool notTrapping_NE = false;
	bool notTrapping_SE = false;
	bool notTrapping_SW = false;

	FDVertex *vert = NULL;
	double dens = 0;
	double area = 0;

	Normalization norm = Normalization();

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		vert = this->vertices.at(iVert);
		ret += vert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea)
			* vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
	}
	return norm.PullLineDensity(ret);
}

void DriftDiffusionSolver::getDeltaXYAtVertex(FDVertex *vert, double &dx, double &dy)
{
	//in the trapping boundary, the length of one or two direction may be zero
	//however, this vertex may not be a real boundary vertex, i.e. the corresponding length may not be zero
	bool notTrapping_NW = false;
	bool notTrapping_NE = false;
	bool notTrapping_SE = false;
	bool notTrapping_SW = false;

	notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NorthwestElem);
	notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->NortheastElem);
	notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SoutheastElem);
	notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, vert->SouthwestElem);

	//West
	if (!(notTrapping_NW && notTrapping_SW))
	{
		dx += vert->WestLength / 2;
	}
	//East
	if (!(notTrapping_NE && notTrapping_SE))
	{
		dx += vert->EastLength / 2;
	}
	//South
	if (!(notTrapping_SE && notTrapping_SW))
	{
		dy += vert->SouthLength / 2;
	}
	//North
	if (!(notTrapping_NE && notTrapping_NW))
	{
		dy += vert->NorthLength / 2;
	}
	SCTM_ASSERT(dx!=0 && dy!=0, 10015);
}

DDTest::DDTest(FDDomain *_domain) : DriftDiffusionSolver(_domain)
{
	
}

void DDTest::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPairVertex;
	std::pair<VertexMapDouble::iterator, bool> insertPairPrpty;
	int vertID = 0;
	int equationID = 0;
	double phyValue = 0;
	FDVertex *currVert = NULL;

	//Previously, this is only a temporary method
	//the following is used to set the mobility to uniform value, because the original setting method leads to incorrect calculation
	//of the mobility at the trapping layer interface.
	//This is now solved by using different method of setting vertex-related physical values.
	//Normalization norm = Normalization();
	//double mobility = MaterialDB::GetMatPrpty(&MaterialDB::Si3N4, MaterialDB::MatProperty::Mat_ElectronMobility);

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	Normalization norm = Normalization();
	double lastDensity = norm.PushDensity(1e12);
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();

		equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, lastDensity));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);
	}
}

void DDTest::setBndCurrent()
{
	FDVertex *currVert = NULL;
	double bcVal_in = 1e-6; //the magnitude of current density vector, the direction is the same with boundary normal direction
	double bcVal_out = 0;
	// current density at boundary for test, in [A/cm^2]
	// note that the current direction is the reversed direction of electron flow

	Normalization norm = Normalization();
	bcVal_in = norm.PushCurrDens(bcVal_in);

	//the sequence of the assignment is in accordance with the direction
	//refresh the boundary condition if needed, determined by the following direction
						bool inNorth = false;
	bool inWest = false;						bool inEast = false;
						bool inSouth = true;
	bool inNorthWest = false; bool inNorthEast = false;
	bool inSouthWest = true; bool inSouthEast = true;

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NorthwestElem);
		bool notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NortheastElem);
		bool notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SoutheastElem);
		bool notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SouthwestElem);

		bool valid_NW = SimpleONO::isValidElem(currVert->NorthwestElem);
		bool valid_NE = SimpleONO::isValidElem(currVert->NortheastElem);
		bool valid_SE = SimpleONO::isValidElem(currVert->SoutheastElem);
		bool valid_SW = SimpleONO::isValidElem(currVert->SouthwestElem);

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
	UtilsTimer.Set();
	prepareSolver(); //call method from base, DriftDiffusionSolver
	refreshCoeffMatrixDueToBC();
	this->matrixSolver.SolveMatrix(rhsVector, this->elecDensity);
	
	UtilsDebug.PrintSparseMatrix(matrixSolver.matrix);
	//UtilsDebug.PrintVector(this->rhsVector, "right hand side vector");
	//UtilsDebug.PrintVector(this->elecDensity, "electron density");
	
	fillBackElecDens();
	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());

	UtilsData.WriteDDResult(this->vertices);
	//SctmFileStream write = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\DDTest\\eDensity.txt", SctmFileStream::Write);
	//write.WriteDDResult(this->vertices, "electron density");
}

void DDTest::setBndDensity()
{
	FDVertex *currVert = NULL;

	//the sequence of the assignment is in accordance with the direction
						bool bindNorth = false;
	bool bindWest = false;						bool bindEast = false;
						bool bindSouth = true;
	//the value of the electron density in each directions, in [cm^-3]
						double valNorth = 0;
	double valWest = 0;							double valEast = 0;
						double valSouth = 1e12;

	Normalization norm = Normalization();
	valNorth = norm.PushDensity(valNorth);
	valSouth = norm.PushDensity(valSouth);
	valEast = norm.PushDensity(valEast);
	valWest = norm.PushDensity(valWest);

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NorthwestElem);
		bool notTrapping_NE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->NortheastElem);
		bool notTrapping_SE = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SoutheastElem);
		bool notTrapping_SW = SimpleONO::isNotElemOf(FDRegion::Trapping, currVert->SouthwestElem);

		//TODO: how to do with the situation when two adjacent boundaries are bound to certain values.
		//South side
		if (( notTrapping_SW && notTrapping_SE ) && bindSouth)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valSouth);
			continue;
		}

		//North side
		if (( notTrapping_NW && notTrapping_NE ) && bindNorth)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valNorth);
			continue;
		}

		//East side
		if (( notTrapping_NE && 
			notTrapping_SE) && bindEast)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valEast);
			continue;
		}
		
		//West side
		if (( notTrapping_NW &&
			notTrapping_SW) && bindWest)
		{
			currVert->BndCond.RefreshBndCond(true, FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valWest);
			continue;
		}

	}
}

void DDTest::refreshBoundary()
{
	//setBndDensity();
	setBndCurrent();
}

void DDTest::refreshCoeffMatrixDueToBC()
{
	FDVertex *vert = NULL;
	int equID = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		vert = this->vertices.at(iVert);
		
		if (vert->IsAtBoundary(FDBoundary::eDensity))
		{
			if (vert->BndCond.GetBCType(FDBoundary::eDensity) == FDBoundary::BC_Dirichlet)
			{
				equID = equationMap[vert->GetID()];
				this->matrixSolver.RefreshRowOfDirichletBC(equID);
			}
		}
	}
}
