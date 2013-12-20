/**
* @file DomainDetails.cpp
* @brief This file contains the detailed implementation of the methods declared in the header file.
*
*
*
* @author
* @version 
* @date 2013-7-4   19:41
* @note
* @todo
*/

#include "DomainDetails.h"
#include "SctmMath.h"
#include "SctmUtils.h"
#include "SctmPhys.h"
using SctmPhys::PhysProperty;

double FDVertex::Distance( FDVertex *vertex1, FDVertex *vertex2 )
{
	return SctmMath::sqrt(SctmMath::square(vertex1->X - vertex2->X) + SctmMath::square(vertex1->Y - vertex2->Y));
}

bool FDVertex::IsAtBoundary(FDBoundary::BCName bcName)
{
	return this->BndCond.Valid(bcName);
}

FDVertex::FDVertex(unsigned _id, double _x, double _y) : X(_x), Y(_y), id(_id)
{
	//These pointers will be set to NULL outside when the setting of adjacency of a vertex 
	//and initializing the contact.
	EastVertex = NULL;
	WestVertex = NULL;
	NorthVertex = NULL;
	SouthVertex = NULL;
	NortheastElem = NULL;
	NorthwestElem = NULL;
	SoutheastElem = NULL;
	SouthwestElem = NULL;
	Contact = NULL;
	Phys = new PhysProperty(this);
	Trap = NULL; // initially no trap property is attached to the vertex, trap property is set in the specific method in FDDomain
}

void FDBoundary::SetBnd(BCName bcName, BCType bcType, VectorValue bndVec, double bcValue /*= 0*/)
{
	SCTM_ASSERT( this->bnd_valid.find(bcName) == this->bnd_valid.end(), 10017 ); //make sure bcName not exists

	if (bndVec.DirectionValid())
	{
		bnd_valid[bcName] = true;
		
		bndVec.Normalize();
		//by default, the bc direction is the same with boundary direction
		bnd_normVec[bcName] = bndVec;
		bc_normVec[bcName] = bndVec;
		
		bc_types[bcName] = bcType;
		bc_values[bcName] = bcValue;
	}
	else
	{
		SCTM_ASSERT(SCTM_ERROR, 10023);
	}
}

FDBoundary::BCType FDBoundary::GetBCType(BCName bcName)
{
	//map<BCName, BCType>::iterator iter;
	//iter = this->bc_types.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_types.end(), 10010);
	//return iter->second;

	//the same with above method
	SCTM_ASSERT(bc_types.find(bcName)!=bc_types.end(), 10010);
	return bc_types[bcName];

}

double FDBoundary::GetBCValue(BCName bcName)
{
	//map<BCName, double>::iterator iter;
	//iter = this->bc_values.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_values.end(), 10010);
	//return iter->second;

	SCTM_ASSERT(bc_values.find(bcName)!=bc_values.end(), 10010);
	return bc_values[bcName];
}

bool FDBoundary::Valid(BCName bcName)
{
	return bnd_valid[bcName];
	//map<BCName, bool>::iterator iter;
	//iter = this->bc_valid.find(bcName);
	//SCTM_ASSERT(iter!=this->bc_valid.end(), 10010);
	//return iter->second;
}

void FDBoundary::RefreshBndCond(BCName bcName, double newValue, VectorValue bcNormVec /*= VectorValue(0, 0)*/)
{
	SCTM_ASSERT(bnd_valid.find(bcName)!=bnd_valid.end(), 10014); //make sure bcName exists

	bc_values[bcName] = newValue;
}

void FDBoundary::RefreshBndCond(BCName bcName, BCType bcType, double bcVal /*= 0*/, VectorValue bcNormVec /*= VectorValue(0, 0)*/)
{
	SCTM_ASSERT(bnd_valid.find(bcName)!=bnd_valid.end(), 10014); //make sure bcName exists
	SCTM_ASSERT( bcType == BC_Dirichlet || bcNormVec.DirectionValid(), 10016 );

	bc_types[bcName] = bcType;
	bc_values[bcName] = bcVal;
	if (bcNormVec.DirectionValid())
	{
		bcNormVec.Normalize();
	}
	bc_normVec[bcName] = bcNormVec;
}

void FDBoundary::RefreshBndCond(BCName bcName, VectorValue bcNormVec)
{
	SCTM_ASSERT(bnd_valid.find(bcName)!=bnd_valid.end(), 10014); //make sure bcName exists
	SCTM_ASSERT( bc_types[bcName] == BC_Dirichlet || bcNormVec.DirectionValid(), 10016 );

	if (bcNormVec.DirectionValid())
	{
		bcNormVec.Normalize();
	}
	bc_normVec[bcName] = bcNormVec;
}

VectorValue & FDBoundary::GetBCNormVector(BCName bcName)
{
	SCTM_ASSERT(bc_normVec.find(bcName)!=bc_normVec.end(), 10010);
	return bc_normVec[bcName];
}

VectorValue & FDBoundary::GetBndDirection(BCName bcName)
{
	SCTM_ASSERT(bnd_normVec.find(bcName)!=bnd_normVec.end(), 10010);
	return bnd_normVec[bcName];
}

void FDBoundary::SetTunnelTag(TunnelTag tag)
{
	SCTM_ASSERT(bnd_valid[eDensity], 10026);
	SCTM_ASSERT(bc_types[eDensity]==BC_Cauchy || bc_types[eDensity]==BC_Neumann, 10026);

	tunTag = tag;
}

FDBoundary::TunnelTag FDBoundary::GetBCTunnelTag()
{
	SCTM_ASSERT(bnd_valid[eDensity], 10026);
	SCTM_ASSERT(bc_types[eDensity]==BC_Cauchy || bc_types[eDensity]==BC_Neumann, 10026);
	return tunTag;
}

FDBoundary::FDBoundary()
{
	tunTag = noTunnel;
}

FDElement::FDElement(unsigned int _id, FDVertex *_swVertex, FDVertex *_seVertex, FDVertex *_neVertex, FDVertex *_nwVertex)
	:id(_id), SouthwestVertex(_swVertex), SoutheastVertex(_seVertex), NortheastVertex(_neVertex), NorthwestVertex(_nwVertex)
{
	WestLength =FDVertex::Distance(NorthwestVertex, SouthwestVertex);
	SouthLength = FDVertex::Distance(SouthwestVertex, SoutheastVertex);
	EastLength = FDVertex::Distance(NortheastVertex, SoutheastVertex);
	NorthLength =FDVertex::Distance(NorthwestVertex, NortheastVertex);
	Area = WestLength * SouthLength;
	///TODO: judge if the corresponding lengths equal to each other
	double relError = (WestLength - EastLength) / WestLength;
	SCTM_ASSERT(relError < 0.01, 10003);
	relError = (SouthLength - NorthLength) / SouthLength;
	SCTM_ASSERT(relError < 0.01, 10003);
}

std::vector<FDVertex *> & FDContact::GetContactVerts()
{
	return vertices;
}
