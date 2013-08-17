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

void TwoDimPoisson::buildCoefficientMatrix()
{
	int matrixSize = this->vertices.size();
	sparseMatrix.resize(matrixSize, matrixSize);
	sparseMatrix.reserve(Eigen::VectorXd::Constant(matrixSize, 5));//reserver room for 5 non-zeros per column

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
		
		//fill the coefficient related to west vertex
		if ( currVert->WestVertex != NULL )
		{
			val = 0;
			indexCoefficient = vertMap[currVert->WestVertex->GetInternalID()];
			if ( currVert->NorthwestElem != NULL)
			{
				epsilon = GetMatPrpty(currVert->NorthwestElem->Region->RegionMaterial, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * currVert->NorthLength / currVert->WestLength;
			}
			if ( currVert->SouthwestElem != NULL)
			{
				epsilon = GetMatPrpty(currVert->SouthwestElem->Region->RegionMaterial, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * currVert->SouthLength / currVert->WestLength;
			}
			sparseMatrix.insert(indexEquation, indexCoefficient) = val;
		}

		//fill the coefficient related to east vertex
		if ( currVert->EastVertex != NULL )
		{
			val = 0;
			indexCoefficient = vertMap[currVert->EastVertex->GetInternalID()];
			if ( currVert->NortheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->NortheastElem->Region->RegionMaterial, MatProperty::Mat_DielectricConstant);
				val = 0.5 * epsilon * currVert->NorthLength / currVert->EastLength;
			}
			if ( currVert->SoutheastElem != NULL )
			{
				epsilon = GetMatPrpty(currVert->SoutheastElem->Region->RegionMaterial, MatProperty::Mat_DielectricConstant);
				val += 0.5 * epsilon * currVert->SouthLength / currVert->EastLength;
			}
		}
		////////////////////////here////////////////////
	}
}

void TwoDimPoisson::prepareVertexMap()
{
	std::pair<MapForVertex::iterator, bool> insertPair;
	//map the vertices in order to get the specific vertex
	int vertID = 0;
	int equationIndex = 0; //equation index is also the vertex index in the vertices container
	FDVertex *currVert = NULL;

	//this map is filled in order to obtain the vertex index from the vertex internal id. This is useful
	//in setting up the equation, i.e. filling the matrix.
	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		vertID = currVert->GetInternalID();
		equationIndex = iVert;
		insertPair = this->vertMap.insert(MapForVertex::value_type(vertID, equationIndex));
		SCTM_ASSERT(insertPair.second==true, 7);//to check if the insertion is successful.
	}
}
