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
#include "SctmPhys.h"
#include "FDDomain.h"

using SctmPhys::PhysProperty;
using MaterialDB::GetMatPrpty;
using MaterialDB::MatProperty;
using SctmUtils::SctmFileStream;

TwoDimPoissonSolver::TwoDimPoissonSolver(FDDomain *domain) :vertices(domain->GetVertices())
{
	initializeSolver();
}

void TwoDimPoissonSolver::initializeSolver()
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
	matrixSolver.matrix.resize(matrixSize, matrixSize);
	matrixSolver.matrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));//reserver room for 5 non-zeros per column

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
		indexEquation = equationMap[currVert->GetID()];
		SCTM_ASSERT(indexEquation==iVert, 10008);

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
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = val;
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
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = val;
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
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = val;
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
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = val;
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
			matrixSolver.matrix.insert(indexEquation, indexCoefficient) = val;
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
		insertPair = this->equationMap.insert(VertexMapInt::value_type(vertID, equationIndex));
		SCTM_ASSERT(insertPair.second==true, 10007);//to check if the insertion is successful.
		equationIndex += 1;
	}
}

void TwoDimPoissonSolver::buildRhsVector()
{
	FDVertex *currVert = NULL;
	double charge = 0;
	int equationID = 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		equationID = equationMap[currVert->GetID()];
		SCTM_ASSERT(equationID==iVert, 10008);
		charge = currVert->Phys->GetPhysPrpty(PhysProperty::NetCharge);
		this->rhsVector.at(equationID) = charge;
	}
}

void TwoDimPoissonSolver::refreshCoefficientMatrix()
{
	FDVertex *currVert = NULL;
	int equationID= 0;
	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		//This vertex is a boundary vertex and it is linked with Dirichlet boundary condition
		//As to the vertex linked with Neumann boundary condition, the equation in the matrix stays unchanged.
		if ( currVert->IsAtBoundary(FDBoundary::Potential) && ( currVert->BndCond.GetBCType(FDBoundary::Potential) == FDBoundary::BC_Dirichlet ) ) 
		{ 
			//the equation index is the same with the position of the vertex in the vertices list
			//so the equation index to set equals to iVert here. Use iVert alternatively.
			equationID = equationMap[currVert->GetID()];
			matrixSolver.RefreshRowOfDirichletBC(equationID);
		}
	}
}

void TwoDimPoissonSolver::refreshRhs()
{
	FDVertex *currVert = NULL;
	int equationID = 0;
	double bcVal = 0;
	double norm_alpha = 0;
	double norm_beta = 0;

	double deltaX = 0;
	double daltaY = 0;
	double eps_dot_dL = 0; // epsilon * deltaX(Y)
	double deltaRhs = 0;

	for (size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = vertices.at(iVert);
		equationID = equationMap[currVert->GetID()];
		if (currVert->IsAtBoundary(FDBoundary::Potential))
		{
			switch (currVert->BndCond.GetBCType(FDBoundary::Potential))
			{
			case FDBoundary::BC_Dirichlet:
				rhsVector.at(equationID) = currVert->BndCond.GetBCValue(FDBoundary::Potential);//for potential           
				break;
			case FDBoundary::BC_Neumann:
				//CAUTION!! The case of BC_Neumann with certain boundary condition value has not been tested and verified.
				bcVal = currVert->BndCond.GetBCValue(FDBoundary::Potential);
				norm_alpha = currVert->BndCond.GetBCNormVector(FDBoundary::Potential).NormX();
				norm_beta = currVert->BndCond.GetBCNormVector(FDBoundary::Potential).NormY();

				if (norm_alpha > 0)
				{
					eps_dot_dL = (currVert->NorthwestElem == NULL) ? 0 :
								0.5 * currVert->NorthLength * GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					eps_dot_dL += (currVert->SouthwestElem == NULL) ? 0 :
								0.5 * currVert->SouthLength * GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				} 
				else if (norm_alpha < 0)
				{
					eps_dot_dL = (currVert->NortheastElem == NULL) ? 0 :
								0.5 * currVert->NorthLength * GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					eps_dot_dL += (currVert->SoutheastElem == NULL) ? 0 :
								0.5 * currVert->SouthLength * GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}
				deltaRhs = bcVal * norm_alpha * eps_dot_dL;

				if (norm_beta > 0)
				{
					eps_dot_dL = (currVert->SouthwestElem == NULL) ? 0 :
								0.5 * currVert->WestLength * GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
					eps_dot_dL = (currVert->SoutheastElem == NULL) ? 0 :
								0.5 * currVert->EastLength * GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}
				else if (norm_beta < 0)
				{
					eps_dot_dL = (currVert->NorthwestElem == NULL) ? 0 :
								0.5 * currVert->WestLength * GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				eps_dot_dL = (currVert->NortheastElem == NULL) ? 0 :
								0.5 * currVert->EastLength * GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_DielectricConstant);
				}
				deltaRhs += bcVal * norm_beta * eps_dot_dL;

				rhsVector.at(equationID) += deltaRhs;
				
				break;
			}
		}
	}
}

void TwoDimPoissonSolver::SolvePotential()
{
	SctmUtils::UtilsTimer.Set();
	//the density of different kinds of charge is updated
	refreshRhs();
	//SctmUtils::UtilsDebug.PrintSparseMatrix(this->matrix);
	//SctmUtils::UtilsDebug.PrintVector(this->rhsVector, "right-hand side");
	matrixSolver.SolveMatrix(rhsVector, potential);
	//SctmUtils::UtilsDebug.PrintVector(this->potential, "potential");
	fillBackPotential();
	SctmUtils::UtilsMsg.PrintTimeElapsed(SctmUtils::UtilsTimer.SinceLastSet());
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
		currVert->Phys->SetPhysPrpty(PhysProperty::ElectrostaticPotential, pot);
	}

	SctmFileStream write = SctmFileStream("E:\\PhD Study\\SimCTM\\SctmTest\\PoissonTest\\potential.txt", SctmFileStream::Write);
	write.WritePoissonResult(this->vertices);
}
