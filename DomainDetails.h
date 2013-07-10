/**
* @file DomainDetails.h
* @brief This file describe the details of the finite differential domain.
* 
* The details of the finite differential domain include the description of 
* vertex, edge, element, region, etc.
*
* @author lunzhy
* @version 
* @date 2013-7-1   15:11
* @note
* @todo
*/
#ifndef _DOMAINDETAIL_H_
#define _DOMAINDETAIL_H_

#include "GeneralMath.h"
#include <vector>
#include "Normalization.h"
using namespace Utility;

class FDElement;
class FDVertex
{
public:
	FDVertex(unsigned _id, double _x, double _y): X(_x), Y(_y), id(_id) {}
	double X;
	double Y;
	double EastLength;
	double WestLength;
	double NorthLength;
	double SouthLength;
	FDVertex *EastVertex;
	FDVertex *WestVertex;
	FDVertex *NorthVertex;
	FDVertex *SouthVertex;
	FDElement *NortheastElem;
	FDElement *NorthwestElem;
	FDElement *SoutheastElem;
	FDElement *SouthwestElem;

	int GetInternalID() { return id; }
	static double Distance(FDVertex *vertex1, FDVertex *vertex2);
protected:
	unsigned int id;
};

/*
class FDEdge
{
public:
	FDEdge(unsigned _id, FDVertex *_vertex1, FDVertex *_vertex2)
		: id(_id), vertex1(_vertex1), vertex2(_vertex2)
	{
		length = GeneralMath::sqrt(GeneralMath::square(_vertex1->x - _vertex2->x) + GeneralMath::square(_vertex1->y - _vertex2->y));
	}
	FDVertex *vertex1;
	FDVertex *vertex2;
	double length;
protected:
	unsigned int id;
};
*/

class FDRegion;
class FDElement
{
public:
	//temporary construction method
	FDElement( unsigned int _id, FDVertex *_swVertex, FDVertex *_seVertex, FDVertex *_neVertex, FDVertex *_nwVertex )
		:id(_id), SouthwestVertex(_swVertex), SoutheastVertex(_seVertex), NortheastVertex(_neVertex), NorthwestVertex(_nwVertex)
	{
		WestLength =FDVertex::Distance(NorthwestVertex, SouthwestVertex);
		SouthLength = FDVertex::Distance(SouthwestVertex, SoutheastVertex);
		EastLength = FDVertex::Distance(NortheastVertex, SoutheastVertex);
		NorthLength =FDVertex::Distance(NorthwestVertex, NortheastVertex);
		//TODO: judge if the corresponding length equal to each other
	}
	
	double WestLength;
	double EastLength;
	double NorthLength;
	double SouthLength;
	FDVertex *NorthwestVertex;
	FDVertex *NortheastVertex;
	FDVertex *SouthwestVertex;
	FDVertex *SoutheastVertex;
	FDRegion *region;

	int GetInternalID() { return id; }
	void SetRegion( FDRegion *region ) {this->region = region;}
protected:
	unsigned int id;
};

class FDRegion
{
public:
	enum RegionType
	{
		TunnelingOxide,
		TrappingLayer,
		BlockingOxide
	};
	FDRegion(unsigned int _id, RegionType _type)
		:id(_id), Type(_type) {}
	RegionType Type;

	void AddElement( FDElement *elem ) {elements.push_back(elem);}
protected:
	std::vector<FDElement *> elements;
	unsigned int id;
};

#endif // !_DOMAINDETAIL_H_