/**
* @file PoissonSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-14   20:10
* @note
* @todo
*/
#include "PoissonSolver.h"
#include "DomainDetails.h"
#include "SctmUtils.h"
#include "Material.h"

using MaterialDB::GetMatPrpty;
using MaterialDB::MatProperty;
using SctmUtils::SctmFileOperator;
using SctmUtils::UtilsTimer;
using SctmUtils::UtilsMsg;

TwoDimPoissonSolver::TwoDimPoissonSolver(FDDomain *domain) :vertices(domain->GetVertices())
{
	prepareSolver();
}

void TwoDimPoissonSolver::prepareSolver()
{
	int vertSize = this->vertices.size();
	this->potential.resize(vertSize);
	this->rhsVector.resize(vertSize);
	buildVertexMap();
	buildCoefficientMatrix();
	buildRhsVector();
	//the structure does not change, so the final coefficient matrix after refreshing does not change either
	refreshCoefficientMatrix();
}

void TwoDimPoissonSolver::buildCoefficientMatrix()
{
	int matrixSize = this->vertices.size();
	matrix.resize(matrixSize, matrixSize);
	matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));//reserver room for 5 non-zeros per column

	//the sequence of the equation and the sequence of the coefficient in each equation are in accordance with the 
	//sequence of the vertices container.
	int indexCoefficient = 0;//the column index of the matrix,  the index of the coefficient of each equation
	int indexEquation = 0;//the row index of the matrix, the index of the equation

	FDVertex *currVert = NULL;
	double epsilon = 0;
	double val = 0;

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = iVert;
		SCTM_ASSERT(indexEquation==equationMap[currVert->GetID()], 10008);

		//fill the coefficient related to west vertex
		if ( currVert->WestVertex != NULL )
		{
			val = 0;
			indexCoefficient = equationMap[currVert->WestVertex->GetID()];
			if ( currVert->NorthwestElem != NULL)
			{
				epsilon = GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * currVert->NorthLength / currVert->WestLength;
			}
			if ( currVert->SouthwestElem != NULL)
			{
				epsilon = GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * currVert->SouthLength / currVert->WestLength;
			}
			matrix.insert(indexEquation, indexCoefficient) = val;
		}

		//fill the coefficient related to east vertex
		if ( currVert->EastVertex != NULL )
		{
			val = 0;
			indexCoefficient = equationMap[currVert->EastVertex->GetID()];
			if ( currVert->NortheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * currVert->NorthLength / currVert->EastLength;
			}
			if ( currVert->SoutheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * currVert->SouthLength / currVert->EastLength;
			}
			matrix.insert(indexEquation, indexCoefficient) = val;
		}

		//fill the coefficient related to the current vertex
		if ( currVert != NULL )
		{
			val = 0;
			indexCoefficient = indexEquation; //indexCoefficient = indexEquation = vertMap[currVert->GetInternalID]
			if ( currVert->NorthwestElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += -0.5 * epsilon * ( currVert->NorthLength / currVert->WestLength + currVert->WestLength / currVert->NorthLength );
			}

			if ( currVert->SouthwestElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += -0.5 * epsilon * ( currVert->SouthLength / currVert->WestLength + currVert->WestLength / currVert->SouthLength ); 
			}

			if ( currVert->SoutheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += -0.5 * epsilon * ( currVert->SouthLength / currVert->EastLength + currVert->EastLength / currVert->SouthLength );
			}

			if ( currVert->NortheastElem != NULL)
			{
				epsilon = GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += -0.5 * epsilon * ( currVert->NorthLength / currVert->EastLength + currVert->EastLength / currVert->NorthLength );
			}
			matrix.insert(indexEquation, indexCoefficient) = val;
		}

		//fill in the coefficient related to the south vertex
		if ( currVert->SouthVertex != NULL )
		{
			val = 0;
			indexCoefficient = equationMap[currVert->SouthVertex->GetID()];
			if ( currVert->SouthwestElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * ( currVert->WestLength / currVert->SouthLength );
			}
			if ( currVert->SoutheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * ( currVert->EastLength / currVert->SouthLength);
			}
			matrix.insert(indexEquation, indexCoefficient) = val;
		}

		//fill in the coefficient related to the north vertex
		if ( currVert->NorthVertex != NULL )
		{
			val = 0;
			indexCoefficient = equationMap[currVert->NorthVertex->GetID()];
			if ( currVert->NortheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * ( currVert->EastLength / currVert->NorthLength);
			}
			if ( currVert->NorthwestElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * ( currVert->WestLength / currVert->NorthLength);
			}
			matrix.insert(indexEquation, indexCoefficient) = val;
		}
	}
}

void TwoDimPoissonSolver::buildVertexMap()
{
	std::pair<VertexMapInt::iterator, bool> insertPair;
	//map the vertices in order to get the specific vertex
	int vertID = 0;
	int equationIndex = 0; //equation index is also the vertex index in the vertices container
	FDVertex *currVert = NULL;

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetID();
		equationIndex = iVert;
		insertPair = this->equationMap.insert(VertexMapInt::value_type(vertID, equationIndex));
		SCTM_ASSERT(insertPair.second==true, 10007);//to check if the insertion is successful.
	}
}

void TwoDimPoissonSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	double charge = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		charge = currVert->Phys.GetPhysPrpty(PhysProperty::NetCharge);
		this->rhsVector.at(iVert) = charge;
	}
}

void TwoDimPoissonSolver::refreshCoefficientMatrix()
{
	//-------------------------------------------------------------------------------
	//the following method is used to iterate the non-zero coefficient of the matrix
	//for (int k=0; k<sparseMatrix.outerSize(); ++k)
	//	for (SparseMatrix<double>::InnerIterator it(sparseMatrix,k); it; ++it)
	//	{
	//		it.valueRef() = 9; // for get the reference of the coefficient
	//		it.row(); // get the row index
	//		it.col(); // get the column index (here it is equal to k)
	//		it.index(); // inner index, here it is equal to it.row()
	//	}
	//--------------------------------------------------------------------------------

	FDVertex *currVert = NULL;
	int equationIndexToSet= 0;
	int coefficientIndexToSet = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//This vertex is a boundary vertex and it is linked with Dirichlet boundary condition
		//As to the vertex linked with Neumann boundary condition, the equation in the matrix stays unchanged.
		if ( currVert->IsAtBoundary(FDBoundary::Potential) && ( currVert->BndCond.GetBCType(FDBoundary::Potential) == FDBoundary::BC_Dirichlet ) ) 
		{ 
			//the equation index is the same with the position of the vertex in the vertices list
			//so the equation index to set equals to iVert here. Use iVert alternatively.
			equationIndexToSet = equationMap[currVert->GetID()];
			coefficientIndexToSet = equationIndexToSet;
			//Although the method of refreshing coefficient value in SparseMatrixSolver class is ready,
			//it cannot be applied here, because refreshing the whole row of the matrix is need here.
			for (int k = 0; k < this->matrix.outerSize(); ++k)
				for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
				{
					if ( it.row() == equationIndexToSet )
					{
						if ( it.col() == coefficientIndexToSet )
							it.valueRef() = 1;
						else
							it.valueRef() = 0;
					}
				}
		}
	}
}

void TwoDimPoissonSolver::refreshRhs()
{
	FDVertex *currVert = NULL;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		if (currVert->IsAtBoundary(FDBoundary::Potential))
		{
			switch (currVert->BndCond.GetBCType(FDBoundary::Potential))
			{
			case FDBoundary::BC_Dirichlet:
				rhsVector.at(iVert) = currVert->BndCond.GetBCValue(FDBoundary::Potential);//for potential           
				break;
			case FDBoundary::BC_Neumann:
			case FDBoundary::BC_Artificial: //BC_Artificial is a special kind of BC_Neumann
				if (currVert->NorthwestElem == NULL)
				{
					if (currVert->NortheastElem != NULL)
						rhsVector.at(iVert) += -0.5 * currVert->BndCond.GetBCValueWestEast(FDBoundary::Potential) * currVert->NorthLength *
											GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					if (currVert->SouthwestElem != NULL)
						rhsVector.at(iVert) += +0.5 * currVert->BndCond.GetBCValueSouthNorth(FDBoundary::Potential) * currVert->WestLength *
											GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);	
				}

				if (currVert->NortheastElem == NULL)
				{
					if (currVert->NorthwestElem != NULL)
						rhsVector.at(iVert) += +0.5 * currVert->BndCond.GetBCValueWestEast(FDBoundary::Potential) * currVert->NorthLength *
											GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					if (currVert->SoutheastElem != NULL)
						rhsVector.at(iVert) += +0.5 * currVert->BndCond.GetBCValueSouthNorth(FDBoundary::Potential) * currVert->EastLength *
											GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}

				if (currVert->SouthwestElem == NULL)
				{
					if (currVert->SoutheastElem != NULL)
						rhsVector.at(iVert) += -0.5 * currVert->BndCond.GetBCValueWestEast(FDBoundary::Potential) * currVert->SouthLength *
											GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					if (currVert->NorthwestElem != NULL)
						rhsVector.at(iVert) += -0.5 * currVert->BndCond.GetBCValueSouthNorth(FDBoundary::Potential) * currVert->WestLength *
											GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}
				if (currVert->SoutheastElem == NULL)
				{
					if (currVert->SouthwestElem != NULL)
						rhsVector.at(iVert) += +0.5 * currVert->BndCond.GetBCValueWestEast(FDBoundary::Potential) * currVert->SouthLength *
											GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					if (currVert->NortheastElem != NULL)
						rhsVector.at(iVert) += -0.5 * currVert->BndCond.GetBCValueSouthNorth(FDBoundary::Potential) * currVert->EastLength *
											GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}
				break;
			}
		}
	}
}

void TwoDimPoissonSolver::SolvePotential()
{
	UtilsTimer.Set();
	refreshRhs();
	//SctmUtils::UtilsDebug.PrintSparseMatrix(this->matrix);
	//SctmUtils::UtilsDebug.PrintVector(this->rhsVector, "right-hand side");
	SolveMatrix(rhsVector, potential);
	//SctmUtils::UtilsDebug.PrintVector(this->potential, "potential");
	fillBackPotential();
	UtilsMsg.PrintTimeElapsed(UtilsTimer.SinceLastSet());
}

void TwoDimPoissonSolver::fillBackPotential()
{
	FDVertex *currVert = NULL;
	double pot = 0;
	int equationID = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		equationID = equationMap[currVert->GetID()];
		pot = potential.at(equationID); //iVert = EquationID
		currVert->Phys.SetPhysPrpty(PhysProperty::ElectrostaticPotential, pot);
	}

	SctmFileOperator write = SctmFileOperator("E:\\PhD Study\\SimCTM\\SctmTest\\PoissonTest\\potential.txt", SctmFileOperator::Write);
	write.WritePoissonResult(this->vertices);
}
