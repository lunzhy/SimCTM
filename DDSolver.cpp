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
#include "DomainDetails.h"
#include "Material.h"
#include "SctmPhys.h"
using namespace SctmPhys;

DriftDiffusionSolver::DriftDiffusionSolver(vector<FDVertex *> _vertices)
	:vertices(_vertices)
{
	prepareSolver();
}

void DriftDiffusionSolver::prepareSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");
	
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->eDensity.resize(vertSize);

	this->T = 300; //temporarily used here
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
		vertID = currVert->GetInternalID();

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

	double KT_div_q = SctmPhys::k0 * this->T / SctmPhys::q;

	for (std::size_t iVert = 0; iVert != this->vertices.size(); ++iVert)
	{
		currVert = this->vertices.at(iVert);
		indexEquation = vertMap[currVert->GetInternalID()];
		SCTM_ASSERT(indexEquation==iVert, 100012);

		coeff_center = 0;

		//the vertex is at current density boundary condition
		if ( currVert->BndCond.Valid(FDBoundary::eCurrentDensity) )
		{

		}
		else // the vertex is NOT at current boundary condition
		{
			//fill the coefficient related to current vertex


			if ( currVert->EastVertex != NULL )
			{
				mobility = (mobilityMap[currVert->EastVertex->GetInternalID()] + mobilityMap[currVert->GetInternalID()]) / 2;
				indexCoefficient = vertMap[currVert->EastVertex->GetInternalID()];
				coeff_adjacent = mobility / currVert->EastLength *
								( (potentialMap[currVert->EastVertex->GetInternalID()] - potentialMap[currVert->GetInternalID()]) / 2
								+ KT_div_q );
				coeff_center += mobility / currVert->EastLength *
								( (potentialMap[currVert->EastVertex->GetInternalID()] - potentialMap[currVert->GetInternalID()]) / 2
								- KT_div_q );

			}

			if ( currVert->WestVertex != NULL )
			{
				
			}

			if ( currVert->NorthVertex != NULL )
			{

			}

			if ( currVert->SouthVertex != NULL )
			{

			}
		}
	}
}
