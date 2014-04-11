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
#include "Normalization.h"
#include "SctmUtils.h"
#include "SctmMath.h"
using namespace SctmPhys;
using namespace SctmUtils;

DriftDiffusionSolver::DriftDiffusionSolver(FDDomain *_domain): domain(_domain), totalVertices(domain->GetVertices())
{
	this->bcMethod = UsingCurrentDensity; // this is the correct boundary condition method
	this->useCrankNicolsonMethod = false; // Crank-Nicolson method isn't complete currently.
	this->useScharfetterGummelMethod = true;
	this->lastTimeStep = 0;
	this->temperature = SctmGlobalControl::Get().Temperature;
	getDDVertices(domain);
	initializeSolver();
}

void DriftDiffusionSolver::SolveDD(VertexMapDouble &bc1, VertexMapDouble &bc2)
{
	//SctmTimer::GetInstance().Set();

	//set the simulation time step
	setTimeStep();
	
	//refresh the vertex map for building Coefficient matrix
	refreshVertexMap();
	
	//build and refresh the coefficient matrix
	buildCoefficientMatrix();
	setCoeffMatrixForTimestep();
	
	//handle the tunneling current. Update the coefficient matrix or refreshing BndCond value for building Rhs vector
	//IMPORTANT! this is done before building the rhs vector
	//for tunneling-in current the boundary condition is revised for building rhs vector
	//for tunneling-out current the coefficient matrix is updated using the tunneling coefficient
	handleBndTunnelCurrDens(bc1, bc2);
	
	//buildRhsVector and refreshRhsWithBC are called together, because for each simulation step, the initial building of Rhs is
	//different due to the difference in last time electron density
	buildRhsVector();

	//dealing the the trapping/detrapping mechanism
	//some of these considerations update the matrix coefficient and some update rhs vector
	//updateCoeffMatrixForTrapping();
	updateRhsForTrapping_ExplicitMethod();
	updateRhsForDetrapping();
	updateRhsForMFNTunneling();

	//solve the matrix
	this->matrixSolver.SolveMatrix(rhsVector, this->elecDensity);
	
	//fill back electron density to last time density, this is also done in refreshing vertex map
	fillBackElecDens();

	//SctmDebug::GetInstance().PrintSparseMatrix(matrixSolver.matrix);
	//SctmMessaging::GetInstance().PrintTimeElapsed(SctmTimer::GetInstance().SinceLastSet());
}

void DriftDiffusionSolver::initializeSolver()
{
	int vertSize = this->ddVertices.size();
	this->rhsVector.resize(vertSize);
	this->elecDensity.resize(vertSize);
	buildVertexMap(); //call the method in DriftDiffusionSolver (Base class), because at this time the derived class is not constructed.
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
	for (std::size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		vertID = currVert->GetID();

		//equationID = iVert; //equation index is also the vertex index in the vertices container
		insertPairVertex = this->equationMap.insert(VertexMapInt::value_type(vertID, equationID));
		SCTM_ASSERT(insertPairVertex.second==true, 10011);

		insertPairPrpty = this->mobilityMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eMobility)));
		SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		//insertPairPrpty = this->potentialMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		//SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		//insertPairPrpty = this->lastElecDensMap.insert(VertexMapDouble::value_type(vertID, currVert->Phys->GetPhysPrpty(PhysProperty::eDensity)));
		//SCTM_ASSERT(insertPairPrpty.second==true, 10011);

		equationID += 1;
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &in_current)
{
}

void DriftDiffusionSolver::buildCoefficientMatrix()
{
	int indexEquation = 0;
	int matrixSize = this->ddVertices.size();
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));

	FDVertex *currVert = NULL;

	for (std::size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);
		//boundary vertex
		if (currVert->IsAtBoundary(FDBoundary::eDensity))
		{
			switch (this->bcMethod)
			{
			case DirectDiscretization:
				setCoeffBCVertex_DirectDiscretization(currVert);
				break;
			case UsingCurrentDensity:
				setCoeffBCVertex_UsingCurrent(currVert);
				break;
			}
			continue;
		}
		//inner vertex
		setCoeffInnerVertex(currVert);
	}
}

void DriftDiffusionSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	int equationID = 0;

	double rhsVal = 0;
	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
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

void DriftDiffusionSolver::setCoeffMatrixForTimestep()
{
	//The boundary conditions are always BC_Cauchy, so there is no need to refresh the matrix for this reason.
	//However, due to that the time step is likely to be different in each simulation step, the time dependent item
	//in the coefficient matrix has to be add to the corresponding item.
	FDVertex *currVert = NULL;
	int indexEquation = 0;
	int indexCoefficient = 0;

	double coeffToAdd = 0; // coefficient for center vertex
	//if ( lastTimeStep == 0 )
	//	coeffToAdd = -1 / timeStep; // from p_n/p_t, p=partial differential
	//else
	//	coeffToAdd = -1 / timeStep - ( -1 / lastTimeStep );
	coeffToAdd = -1 / timeStep;

	for (std::size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);

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

		bool notTrapping_NW = FDDomain::isNotTrappingElem(currVert->NorthwestElem);
		bool notTrapping_NE = FDDomain::isNotTrappingElem(currVert->NortheastElem);
		bool notTrapping_SE = FDDomain::isNotTrappingElem(currVert->SoutheastElem);
		bool notTrapping_SW = FDDomain::isNotTrappingElem(currVert->SouthwestElem);

		if (!(notTrapping_NE && notTrapping_NW && notTrapping_SE && notTrapping_SW))
		{
			this->ddVertices.push_back(currVert);
		}
	}
}

void DriftDiffusionSolver::setTimeStep()
{
 	timeStep = SctmTimeStep::Get().TimeStep();
}

void DriftDiffusionSolver::UpdateElecDens()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double edens = 0;
	int VertID = 0;
	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		VertID = currVert->GetID();
		equationID = equationMap[VertID];
		edens = this->elecDensity.at(equationID);
		currVert->Phys->SetPhysPrpty(PhysProperty::eDensity, edens);
		
		//It is also essential to refresh property map. this is done at the preparation of the solver.
		//lastElecDensMap[VertID] = edens;
	}
}

void DriftDiffusionSolver::fillBackElecDens()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double edens = 0;
	int VertID = 0;
	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		VertID = currVert->GetID();
		equationID = equationMap[VertID];
		edens = this->elecDensity.at(equationID);
		//currVert->Phys->SetPhysPrpty(PhysProperty::eDensity, edens);
		//It is also essential to refresh property map.
		lastElecDensMap[VertID] = edens;
	}
}

void DriftDiffusionSolver::setCoeffBCVertex_DirectDiscretization(FDVertex *vert)
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
		{
			coeff_center = 1;
			indexEquation = equationMap[vert->GetID()];
			indexCoeff = equationMap[vert->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoeff) = coeff_center;
			break;
		}


		case FDBoundary::BC_Cauchy:
		{
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

}

void DriftDiffusionSolver::setCoeffInnerVertex(FDVertex *vert)
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
	if (useScharfetterGummelMethod)
	{
		coeff_adjacent = mobility / vert->EastLength / deltaX * 
			SctmMath::Bernoulli_Potential(potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]);
		coeff_center += - mobility / vert->EastLength / deltaX *
			SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->EastVertex->GetID()]);
	}
	else
	{
		coeff_adjacent = - mobility / vert->EastLength / deltaX *
			( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
		coeff_center += - mobility / vert->EastLength / deltaX *
			( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
	}
	
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
	if (useScharfetterGummelMethod)
	{
		coeff_adjacent = mobility / vert->WestLength / deltaX * 
			SctmMath::Bernoulli_Potential(potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]);
		coeff_center += - mobility / vert->WestLength /deltaX *
			SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->WestVertex->GetID()]);
	}
	else
	{
		coeff_adjacent = - mobility / vert->WestLength / deltaX *
			( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 - 1 );
		coeff_center += - mobility / vert->WestLength / deltaX *
			( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 + 1 );
	}
	
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
	if (useScharfetterGummelMethod)
	{
		coeff_adjacent = mobility / vert->NorthLength / deltaY * 
			SctmMath::Bernoulli_Potential(potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]);
		coeff_center += - mobility / vert->NorthLength / deltaY *
			SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->NorthVertex->GetID()]);
	}
	else
	{
		coeff_adjacent = - mobility / vert->NorthLength / deltaY *
			( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
		coeff_center += - mobility / vert->NorthLength / deltaY *
			( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
	}

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
	if (useScharfetterGummelMethod)
	{
		coeff_adjacent = mobility / vert->SouthLength / deltaY * 
			SctmMath::Bernoulli_Potential(potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]);
		coeff_center += - mobility / vert->SouthLength / deltaY *
			SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->SouthVertex->GetID()]);
	}
	else
	{
		coeff_adjacent = - mobility / vert->SouthLength / deltaY *
			( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
		coeff_center += - mobility / vert->SouthLength / deltaY *
			( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
	}

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

void DriftDiffusionSolver::processBndCond()
{
	//current nothing is done here.
}

void DriftDiffusionSolver::setCoeffBCVertex_UsingCurrent(FDVertex *vert)
{
	int indexEquation = 0;
	int indexCoefficient = 0;

	double mobility = 0;
	double coeff_adjacent = 0; // coefficient for adjacent vertex, could be west, east, south, north
	double coeff_center = 0; // coefficient for center vertex

	double deltaX = 0; // dJ/dx
	double deltaY = 0; // dJ/dy

	double bndNorm_alpha = 0;
	double bndNorm_beta = 0;

	//Notice. can not set the row of matrix according to the boundary normal direction, this is because the boundary condition
	//at the corner has been specially treated, which can not be applied to setting the matrix.
	//double norm_alpha =  vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormX();
	//double norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormY();

	//we should consider the normalized current equation
	//for electron, J = -u[ n * p_phi/p_x - p_n/p_x ]

	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
		case FDBoundary::BC_Dirichlet:
		{
			coeff_center = 1;
			indexEquation = equationMap[vert->GetID()];
			indexCoefficient = equationMap[vert->GetID()];
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_center;
			break;
		}
		
		case FDBoundary::BC_Cauchy:
		{
			indexEquation = equationMap[vert->GetID()];

			bndNorm_alpha = vert->BndCond.GetBndDirection(FDBoundary::eDensity).X();
			bndNorm_beta = vert->BndCond.GetBndDirection(FDBoundary::eDensity).Y();

			getDeltaXYAtVertex(vert, deltaX, deltaY);

			//has west vertex
			if ( bndNorm_alpha >= 0 )
			{
				mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				indexCoefficient = equationMap[vert->WestVertex->GetID()];
				if (useScharfetterGummelMethod)
				{
					coeff_adjacent = mobility / vert->WestLength / deltaX * 
						SctmMath::Bernoulli_Potential(potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]);
					coeff_center += - mobility / vert->WestLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->WestVertex->GetID()]);
				}
				else
				{
					coeff_adjacent = - mobility / vert->WestLength / deltaX *
						( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
					coeff_center += - mobility / vert->WestLength / deltaX *
						( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
				}

				if (useCrankNicolsonMethod)
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
				}
				else
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
				}

			}

			//has east vertex
			if ( bndNorm_alpha <= 0 )
			{
				mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				indexCoefficient = equationMap[vert->EastVertex->GetID()];
				if (useScharfetterGummelMethod)
				{
					coeff_adjacent = mobility / vert->EastLength / deltaX * 
						SctmMath::Bernoulli_Potential(potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]);
					coeff_center += - mobility / vert->EastLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->EastVertex->GetID()]);
				}
				else
				{
					coeff_adjacent = - mobility / vert->EastLength / deltaX *
						( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
					coeff_center += - mobility / vert->EastLength / deltaX *
						( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
				}

				if (useCrankNicolsonMethod)
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
				}
				else
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
				}
			}

			//has south vertex
			if ( bndNorm_beta >= 0 )
			{
				mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				indexCoefficient = equationMap[vert->SouthVertex->GetID()];
				if (useScharfetterGummelMethod)
				{
					coeff_adjacent = mobility / vert->SouthLength / deltaY * 
						SctmMath::Bernoulli_Potential(potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]);
					coeff_center += - mobility / vert->SouthLength  / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->SouthVertex->GetID()]);
				}
				else
				{
					coeff_adjacent = - mobility / vert->SouthLength / deltaY *
						( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
					coeff_center += - mobility / vert->SouthLength / deltaY *
						( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
				}

				if (useCrankNicolsonMethod)
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent / 2;
				}
				else
				{
					matrixSolver.matrix.insert(indexEquation, indexCoefficient) = coeff_adjacent;
				}
			}

			//has north vertex
			if ( bndNorm_beta <= 0 )
			{
				mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				indexCoefficient = equationMap[vert->NorthVertex->GetID()];
				if (useScharfetterGummelMethod)
				{
					coeff_adjacent = mobility / vert->NorthLength / deltaY * 
						SctmMath::Bernoulli_Potential(potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]);
					coeff_center += - mobility / vert->NorthLength / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->NorthVertex->GetID()]);
				}
				else
				{
					coeff_adjacent = - mobility / vert->NorthLength / deltaY *
						( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 );
					coeff_center += - mobility / vert->NorthLength / deltaY *
						( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 );
				}

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
	
	double east = 0;
	double west = 0;
	double south = 0;
	double north = 0;

	double retVal = 0;

	//related to time step
	rhsTime = -1 * lastElecDensMap[vert->GetID()] / timeStep;
 
	getDeltaXYAtVertex(vert, deltaX, deltaY);

	if (useCrankNicolsonMethod)
	{
		//for Crank-Nicolson discretization method
		//related to east vertex
		//Scharfetter-Gummel Method is considered only when Crank-Nicolson method is used.
		if (useScharfetterGummelMethod)
		{
			mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			east = - mobility / vert->EastLength / deltaX *
				SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->EastVertex->GetID()]) *
				lastElecDensMap[vert->GetID()];

			east += mobility / vert->EastLength / deltaX *
				SctmMath::Bernoulli_Potential(potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) * 
				lastElecDensMap[vert->EastVertex->GetID()];

			//related to west vertex
			mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			west = - mobility / vert->WestLength / deltaX *
				SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->WestVertex->GetID()]) *
				lastElecDensMap[vert->GetID()];

			west += mobility / vert->WestLength / deltaX *
				SctmMath::Bernoulli_Potential(potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) *
				lastElecDensMap[vert->WestVertex->GetID()];

			//related to north vertex
			mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			north = - mobility / vert->NorthLength / deltaY *
				SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->NorthVertex->GetID()]) *
				lastElecDensMap[vert->GetID()];

			north += mobility / vert->NorthLength / deltaY *
				SctmMath::Bernoulli_Potential(potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) *
				lastElecDensMap[vert->NorthVertex->GetID()];

			//related to south vertex
			mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			south = - mobility / vert->SouthLength / deltaY *
				SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->SouthVertex->GetID()]) *
				lastElecDensMap[vert->GetID()];

			south += mobility / vert->SouthLength / deltaY *
				SctmMath::Bernoulli_Potential(potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) *
				lastElecDensMap[vert->SouthVertex->GetID()];
		}
		else
		{
			mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			east = - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];

			east += - mobility / vert->EastLength / deltaX *
				( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) * 
				lastElecDensMap[vert->EastVertex->GetID()];

			//related to west vertex
			mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			west = - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];

			west += - mobility / vert->WestLength / deltaX *
				( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->WestVertex->GetID()];

			//related to north vertex
			mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			north = - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];

			north += - mobility / vert->NorthLength / deltaY *
				( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->NorthVertex->GetID()];

			//related to south vertex
			mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
			south = - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
				lastElecDensMap[vert->GetID()];

			south += - mobility / vert->SouthLength / deltaY *
				( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
				lastElecDensMap[vert->SouthVertex->GetID()];
		}
		
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

	double bndNorm_alpha = 0;
	double bndNorm_beta = 0;

	double mobility = 0;
	double east = 0;
	double west = 0;
	double north = 0;
	double south = 0;

	static VertexMapDouble lastBCMap = VertexMapDouble(); // to store the overall current density of boundary conditions

	switch (vert->BndCond.GetBCType(FDBoundary::eDensity))
	{
		case FDBoundary::BC_Dirichlet:
		{
			//in BC_Dirichlet, the the rhs of boundary vertex is directly set.
			retVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
			break;
		}
		
		case FDBoundary::BC_Cauchy:
		{
			//for Cauchy boundary condition
			//calculation of the addend related current simulation time step
			rhsTime= -1 * lastElecDensMap[vert->GetID()] / timeStep;

			bndNorm_alpha = vert->BndCond.GetBndDirection(FDBoundary::eDensity).X();
			bndNorm_beta = vert->BndCond.GetBndDirection(FDBoundary::eDensity).Y();

			getDeltaXYAtVertex(vert, deltaX, deltaY);

			//has west vertex
			if ( bndNorm_alpha >= 0 )
			{
				mobility = (mobilityMap[vert->WestVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				if (useScharfetterGummelMethod)
				{
					west = - mobility / vert->WestLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->WestVertex->GetID()]) *
						lastElecDensMap[vert->GetID()];

					west += mobility / vert->WestLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) *
						lastElecDensMap[vert->WestVertex->GetID()];
				}
				else
				{
					west += - mobility / vert->WestLength / deltaX *
						( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 - 1 ) *
						lastElecDensMap[vert->WestVertex->GetID()];

					west += - mobility / vert->WestLength / deltaX *
						( (potentialMap[vert->WestVertex->GetID()] - potentialMap[vert->GetID()]) /2 + 1 ) *
						lastElecDensMap[vert->GetID()];
				}
			}

			//has east vertex
			if ( bndNorm_alpha <= 0 )
			{
				mobility = (mobilityMap[vert->EastVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				if (useScharfetterGummelMethod)
				{
					east = - mobility / vert->EastLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->EastVertex->GetID()]) *
						lastElecDensMap[vert->GetID()];

					east += mobility / vert->EastLength / deltaX *
						SctmMath::Bernoulli_Potential(potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) * 
						lastElecDensMap[vert->EastVertex->GetID()];
				}
				else
				{
					east += - mobility / vert->EastLength / deltaX *
						( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
						lastElecDensMap[vert->EastVertex->GetID()];

					east += - mobility / vert->EastLength / deltaX *
						( (potentialMap[vert->EastVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
						lastElecDensMap[vert->GetID()];
				}
			}

			//has south vertex
			if ( bndNorm_beta >= 0 )
			{
				mobility = (mobilityMap[vert->SouthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				if (useScharfetterGummelMethod)
				{
					south = - mobility / vert->SouthLength / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->SouthVertex->GetID()]) *
						lastElecDensMap[vert->GetID()];

					south += mobility / vert->SouthLength / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) *
						lastElecDensMap[vert->SouthVertex->GetID()];
				}
				else
				{
					south += - mobility / vert->SouthLength / deltaY *
						( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
						lastElecDensMap[vert->SouthVertex->GetID()];

					south += - mobility / vert->SouthLength / deltaY *
						( (potentialMap[vert->SouthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
						lastElecDensMap[vert->GetID()];
				}
			}

			//has north vertex
			if ( bndNorm_beta <= 0 )
			{
				mobility = (mobilityMap[vert->NorthVertex->GetID()] + mobilityMap[vert->GetID()]) / 2;
				if (useScharfetterGummelMethod)
				{
					north = - mobility / vert->NorthLength / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->GetID()] - potentialMap[vert->NorthVertex->GetID()]) *
						lastElecDensMap[vert->GetID()];

					north += mobility / vert->NorthLength / deltaY *
						SctmMath::Bernoulli_Potential(potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) *
						lastElecDensMap[vert->NorthVertex->GetID()];
				}
				else
				{
					north += - mobility / vert->NorthLength / deltaY *
						( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 - 1 ) *
						lastElecDensMap[vert->NorthVertex->GetID()];

					north += - mobility / vert->NorthLength / deltaY *
						( (potentialMap[vert->NorthVertex->GetID()] - potentialMap[vert->GetID()]) / 2 + 1 ) *
						lastElecDensMap[vert->GetID()];
				}
			}

			double bcVal = vert->BndCond.GetBCValue(FDBoundary::eDensity);
			double norm_alpha =  vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormX();
			double norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormY();
			//p_J / p_x + p_J / p_y = (Je - Jw) / dx + (Jn - Js) / dy
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
		
	}
	return retVal;
}

double DriftDiffusionSolver::CalculateTotalLineDensity()
{
	double ret = 0;

	FDVertex *vert = NULL;

	Normalization norm = Normalization(this->temperature);

	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		vert = this->ddVertices.at(iVert);
		ret += vert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea)
			* vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
	}
	return norm.PullLineDensity(ret);
}

void DriftDiffusionSolver::getDeltaXYAtVertex(FDVertex *vert, double &dx, double &dy)
{
	//in the trapping boundary, the length of one or two direction may be zero
	//however, this vertex may not be a real boundary vertex, i.e. the corresponding length may not be zero
	dx = 0;
	dy = 0;

	double bndNorm_alpha = 0;
	double bndNorm_beta = 0;
	if ( vert->IsAtBoundary(FDBoundary::eDensity) )
	{
		bndNorm_alpha = vert->BndCond.GetBndDirection(FDBoundary::eDensity).X();
		bndNorm_beta = vert->BndCond.GetBndDirection(FDBoundary::eDensity).Y();
	}
	//so for inner vertex, bndNorm_alpha = 0, bndNorm_beta = 0

	//has west vertex
	if ( bndNorm_alpha >= 0 )
	{
		dx += vert->WestLength / 2;
	}
	//has east vertex
	if ( bndNorm_alpha <= 0)
	{
		dx += vert->EastLength / 2;
	}
	//has south vertex
	if ( bndNorm_beta >=0 )
	{
		dy += vert->SouthLength / 2;
	}
	//has north vertex
	if ( bndNorm_beta <= 0 )
	{
		dy += vert->NorthLength / 2;
	}
	
	SCTM_ASSERT(dx!=0 && dy!=0, 10015);
}

void DriftDiffusionSolver::handleCurrDensBC_in(FDVertex *vert, double currdens)
{
	int vertID = 0;

	SCTM_ASSERT(vert->IsAtBoundary(FDBoundary::eDensity), 10022);
	SCTM_ASSERT(vert->BndCond.GetBCType(FDBoundary::eDensity) == FDBoundary::BC_Cauchy, 10022);

	vert->BndCond.RefreshBndCond(FDBoundary::eDensity, currdens);
	//it should be noticed that for tunneling-in electron current, the current density direction is the 
	//same with the boundary condition, so currDens should be positive value.
	//and the method for building Rhs vector will handle the in-tunneling current density
}

void DriftDiffusionSolver::refreshVertexMap()
{
	int vertID = 0;
	FDVertex *currVert = NULL;
	double density = 0;
	double pot = 0;
	for (VertexMapInt::iterator it = equationMap.begin(); it != equationMap.end(); ++it)
	{
		vertID = it->first;
		currVert = this->domain->GetVertex(vertID);

		density = currVert->Phys->GetPhysPrpty(PhysProperty::eDensity);
		lastElecDensMap[vertID] = density;

		pot = currVert->Phys->GetPhysPrpty(PhysProperty::ElectrostaticPotential);
		potentialMap[vertID] = pot;
	}
}

void DriftDiffusionSolver::handleCurrDensBC_out(FDVertex *vert, double tunCoeff)
{
	//the method for getting rhs value of bc vertex did nothing with tunneling-out boundary condition, because the tunneling-out current value
	//is 0 in that method.
	//the value of tunneling coefficient is in [A*cm]

	int vertID = 0;
	int equID = 0;
	double deltaX = 0;
	double deltaY = 0;
	double norm_alpha = 0;
	double norm_beta = 0;
	double coeffToAdd = 0;


	SCTM_ASSERT(vert->IsAtBoundary(FDBoundary::eDensity), 10022);
	SCTM_ASSERT(vert->BndCond.GetBCType(FDBoundary::eDensity) == FDBoundary::BC_Cauchy, 10022);

	//currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, tunCoeff);
	vertID = vert->GetID();
	equID = equationMap[vertID];
	getDeltaXYAtVertex(vert, deltaX, deltaY);

	//use the boundary condition direction, not the boundary direction
	norm_alpha =  vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormX();
	norm_beta = vert->BndCond.GetBCNormVector(FDBoundary::eDensity).NormY();

	//p_J / p_x = (Je - Jw) / dx 
	//p_J / p_y = (Jn - Js) / dy
	//here, the vector value of the tunneling-out current density is tunCoeff*(norm_alpha, norm_beta) with negative tunCoeff
	//for positive value, coeffToAdd = tunCoeff * norm_alpha (beta)
	//for negative value, coeffToAdd = (-tunCoeff) * (-norm_alpha)
	if (norm_alpha > 0) //lack of east part of the above equation
	{
		coeffToAdd +=  tunCoeff * norm_alpha / deltaX;
	}
	else if (norm_alpha < 0) //lack of west part of the above equation
	{
		//this sign is derived from (-Js)
		coeffToAdd += - tunCoeff * norm_alpha / deltaX;
	}

	if (norm_beta > 0) // lack of north part of the above equation
	{
		coeffToAdd += tunCoeff * norm_beta / deltaY;
	}
	else if (norm_beta < 0) //  lack of south part of the above equation
	{
		coeffToAdd += - tunCoeff * norm_beta / deltaY;
	}

	this->matrixSolver.RefreshMatrixValue(equID, equID, coeffToAdd, SctmSparseMatrixSolver::Add);
}

void DriftDiffusionSolver::handleBndTunnelCurrDens(VertexMapDouble &bc1, VertexMapDouble &bc2)
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	VertexMapDouble bcMap;
	bcMap.insert(bc1.begin(), bc1.end());
	bcMap.insert(bc2.begin(), bc2.end());
	for (VertexMapDouble::iterator it = bcMap.begin(); it != bcMap.end(); ++it)
	{
		vertID = it->first;
		currVert = domain->GetVertex(vertID);
		FDBoundary::TunnelTag tunTag = currVert->BndCond.GetBCTunnelTag();
		SCTM_ASSERT(tunTag!=FDBoundary::noTunnel, 10027);

		switch (tunTag)
		{
			case FDBoundary::eTunnelIn:
			{
				handleCurrDensBC_in(currVert, it->second);
				break;
			}
			case FDBoundary::eTunnelOut:
			{
				handleCurrDensBC_out(currVert, it->second);
				break;
			}

		}
	}
}

void DriftDiffusionSolver::updateCoeffMatrixForTrapping()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	int indexEqu = 0;
	int indexCoeff= 0;
	double coeff_trapping = 0;

	static string captureModel = SctmGlobalControl::Get().TrapCaptureModel;

	for (size_t iVert = 0; iVert != ddVertices.size(); ++iVert)
	{
		currVert = ddVertices.at(iVert);
		vertID = currVert->GetID();

		indexEqu = equationMap[vertID];
		SCTM_ASSERT(indexEqu==iVert, 10012);
		indexCoeff = indexEqu;

		//notice the negative sign
		if (captureModel == "J-Model")
		{
			coeff_trapping = -currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_J_Model) *
				currVert->Trap->GetTrapPrpty(TrapProperty::eEmptyTrapDens);
		}
		else if (captureModel == "V-Model")
		{
			coeff_trapping = -currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_V_Model) *
				currVert->Trap->GetTrapPrpty(TrapProperty::eEmptyTrapDens);
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10036);
		}

		matrixSolver.RefreshMatrixValue(indexEqu, indexCoeff, coeff_trapping, SctmSparseMatrixSolver::Add);
	}
}

void DriftDiffusionSolver::updateRhsForDetrapping()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	int equID = 0;

	double eEmission_SRH = 0;
	double eEmission_PF = 0;
	double eTrappedDens = 0;
	
	double rhs_detrapping = 0;

	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		vertID = currVert->GetID();
		equID = equationMap[vertID];

		eEmission_SRH = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_BasicSRH);

		if (SctmGlobalControl::Get().PhysicsPFModel == "Frequency")
		{
			//consider Poole-Frenkel detrapping
			eEmission_PF = currVert->Trap->GetTrapPrpty(TrapProperty::eEmissionCoeff_PF);
		}

		eTrappedDens = currVert->Trap->GetTrapPrpty(TrapProperty::eTrapped);

		rhs_detrapping = (eEmission_SRH + eEmission_PF) * eTrappedDens;
		// the negative sigh symbolizes moving the addend from right to left of the equation.
		rhsVector.at(equID) += -rhs_detrapping;
	}
}

void DriftDiffusionSolver::updateRhsForMFNTunneling()
{
	double deltaX = 0;
	double deltaY = 0;

	double eCurrDens_MFN_X = 0;
	double eCurrDens_MFN_Y = 0;
	double rhs_mfn = 0;

	FDVertex *currVert = NULL;
	int vertID = 0;
	int equID = 0;
	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		vertID = currVert->GetID();
		equID = equationMap[vertID];

		getDeltaXYAtVertex(currVert, deltaX, deltaY);
		eCurrDens_MFN_X = currVert->Phys->GetPhysPrpty(PhysProperty::eCurrDensMFN_X);
		eCurrDens_MFN_Y = currVert->Phys->GetPhysPrpty(PhysProperty::eCurrDensMFN_Y);

		//if the electron current density is negative, the electrons flow into the vertex
		//TODO: A temporary method is used here, because the sign is related to the current flow direction.
		rhs_mfn = -(eCurrDens_MFN_X / deltaX + eCurrDens_MFN_Y / deltaY);

		// the negative sigh symbolizes moving the addend from right to left of the equation.
		rhsVector.at(equID) += -rhs_mfn;
	}
}

void DriftDiffusionSolver::updateRhsForTrapping_ExplicitMethod()
{
	FDVertex *currVert = NULL;
	int vertID = 0;
	int equIndex = 0;
	double eDensLastTime = 0;
	double rhs_trapping = 0;

	static string captureModel = SctmGlobalControl::Get().TrapCaptureModel;

	for (size_t iVert = 0; iVert != ddVertices.size(); ++iVert)
	{
		currVert = ddVertices.at(iVert);
		vertID = currVert->GetID();

		equIndex = equationMap[vertID];
		SCTM_ASSERT(equIndex == iVert, 10012);
		
		eDensLastTime = this->lastElecDensMap[vertID];

		//notice the negative sign
		if (captureModel == "J-Model")
		{
			rhs_trapping = - currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_J_Model) *
				currVert->Trap->GetTrapPrpty(TrapProperty::eEmptyTrapDens) *
				eDensLastTime;
		}
		else if (captureModel == "V-Model")
		{
			rhs_trapping = - currVert->Trap->GetTrapPrpty(TrapProperty::eCaptureCoeff_V_Model) *
				currVert->Trap->GetTrapPrpty(TrapProperty::eEmptyTrapDens) *
				eDensLastTime;
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10036);
		}
		
		//the negative sign means moving the addend from right to left of the equation
		rhsVector.at(equIndex) = -rhs_trapping;
	}
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
	//double mobility = MaterialDB::GetMatPrpty(&MaterialDB::Si3N4, MaterialDB::MatProperty::Mat_ElectronMobility);

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	Normalization norm = Normalization(this->temperature);
	double lastDensity = norm.PushDensity(1e12);
	for (std::size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
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
	double bcVal_in = 1e-3; //the magnitude of current density vector, the direction is the same with boundary normal direction
	double bcVal_out = 0;
	// current density at boundary for test, in [A/cm^2]
	// note that the current direction is the reversed direction of electron flow

	Normalization norm = Normalization(this->temperature);
	bcVal_in = norm.PushCurrDens(bcVal_in);

	//the sequence of the assignment is in accordance with the direction
	//refresh the boundary condition if needed, determined by the following direction
						bool inNorth = false;
	bool inWest = false;						bool inEast = false;
						bool inSouth = true;
	bool inNorthWest = false; bool inNorthEast = false;
	bool inSouthWest = true; bool inSouthEast = true;

	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = FDDomain::isNotTrappingElem(currVert->NorthwestElem);
		bool notTrapping_NE = FDDomain::isNotTrappingElem(currVert->NortheastElem);
		bool notTrapping_SE = FDDomain::isNotTrappingElem(currVert->SoutheastElem);
		bool notTrapping_SW = FDDomain::isNotTrappingElem(currVert->SouthwestElem);

		bool valid_NW = FDDomain::isValidElem(currVert->NorthwestElem);
		bool valid_NE = FDDomain::isValidElem(currVert->NortheastElem);
		bool valid_SE = FDDomain::isValidElem(currVert->SoutheastElem);
		bool valid_SW = FDDomain::isValidElem(currVert->SouthwestElem);

		//Southeast corner
		if (( !notTrapping_NW && notTrapping_NE && 
			notTrapping_SW && notTrapping_SE ) && inSouthEast)
		{
			if (			!valid_NE &&
				valid_SW)
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if (			valid_NE &&
				!valid_SW)
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Southwest corner
		if (( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE) && inSouthWest)
		{
			if ( valid_NW &&
				!valid_SE )
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if ( !valid_NW &&
				valid_SE )
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Northwest corner
		if (( notTrapping_NW && notTrapping_NE &&
			notTrapping_SW && !notTrapping_SE ) && inNorthWest)
		{
			if (              valid_NE &&
				!valid_SW )
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if (              !valid_NE &&
				valid_SW )
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//Northeast corner
		if (( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE ) && inNorthEast)
		{
			if ( valid_NW &&
				!valid_SE)
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			if ( !valid_NW &&
				valid_SE)
			{
				currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
				continue;
			}
			//when the two adjacent neighbors are both valid (other region) or invalid, the current density is considered to be along the diagonal
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//North side
		if (( notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && !notTrapping_SE ) && inNorth)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//South side
		if (( !notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && notTrapping_SE ) && inSouth)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//East side
		if (( !notTrapping_NW && notTrapping_NE && 
			!notTrapping_SW && notTrapping_SE ) && inEast)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

		//West side
		if (( notTrapping_NW && !notTrapping_NE && 
			notTrapping_SW && !notTrapping_SE) && inWest)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, bcVal_in);
			continue;
		}

	}
}

void DDTest::SolveDD()
{
	SctmTimer::Get().Set();
	//prepareSolver(); //call method from base, DriftDiffusionSolver
	refreshCoeffMatrixDueToBC();
	this->matrixSolver.SolveMatrix(rhsVector, this->elecDensity);
	
	SctmDebug::Get().PrintSparseMatrix(matrixSolver.matrix);
	//SctmDebug::GetInstance().PrintVector(this->rhsVector, "right hand side vector");
	//SctmDebug::GetInstance().PrintVector(this->elecDensity, "electron density");
	
	UpdateElecDens();
	SctmMessaging::Get().PrintTimeElapsed(SctmTimer::Get().PopLastSet());

	SctmData::Get().WriteElecDens(this->ddVertices);
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

	Normalization norm = Normalization(this->temperature);
	valNorth = norm.PushDensity(valNorth);
	valSouth = norm.PushDensity(valSouth);
	valEast = norm.PushDensity(valEast);
	valWest = norm.PushDensity(valWest);

	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		currVert = this->ddVertices.at(iVert);
		//for the current tunneling from tunneling oxide

		bool notTrapping_NW = FDDomain::isNotTrappingElem(currVert->NorthwestElem);
		bool notTrapping_NE = FDDomain::isNotTrappingElem(currVert->NortheastElem);
		bool notTrapping_SE = FDDomain::isNotTrappingElem(currVert->SoutheastElem);
		bool notTrapping_SW = FDDomain::isNotTrappingElem(currVert->SouthwestElem);

		//TODO: how to do with the situation when two adjacent boundaries are bound to certain values.
		//South side
		if (( notTrapping_SW && notTrapping_SE ) && bindSouth)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valSouth);
			continue;
		}

		//North side
		if (( notTrapping_NW && notTrapping_NE ) && bindNorth)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valNorth);
			continue;
		}

		//East side
		if (( notTrapping_NE && 
			notTrapping_SE) && bindEast)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valEast);
			continue;
		}
		
		//West side
		if (( notTrapping_NW &&
			notTrapping_SW) && bindWest)
		{
			currVert->BndCond.RefreshBndCond(FDBoundary::eDensity, FDBoundary::BC_Dirichlet, valWest);
			continue;
		}

	}
}

void DDTest::processBndCond()
{
	//setBndDensity();
	setBndCurrent();
}

void DDTest::refreshCoeffMatrixDueToBC()
{
	FDVertex *vert = NULL;
	int equID = 0;
	for (size_t iVert = 0; iVert != this->ddVertices.size(); ++iVert)
	{
		vert = this->ddVertices.at(iVert);
		
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
