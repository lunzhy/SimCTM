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

DriftDiffusionSolver::DriftDiffusionSolver(vector<FDVertex *> _vertices)
	:vertices(_vertices)
{
	prepareSolver();
}

void DriftDiffusionSolver::prepareSolver()
{
	SctmUtils::UtilsMsg.PrintHeader("Solving Drift-Diffusion");
	
	q = SctmPhys::q;
	int vertSize = this->vertices.size();
	this->rhsVector.resize(vertSize);
	this->eDensity.resize(vertSize);
}

void DriftDiffusionSolver::buildVertexMap()
{
	std::pair<MapForVertex::iterator, bool> insertPairVertex;
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
	}
}

void DriftDiffusionSolver::setBndCondCurrent(vector<double> &current)
{

}
