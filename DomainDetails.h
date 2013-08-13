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

#include <vector>
#include "Normalization.h"
#include "SctmPhys.h"
#include "Material.h"
#include "SctmUtils.h"

using SctmPhys::PhysProperty;
using MaterialDB::Material;

class FDElement;
/// @brief FDVertex is the class describing the vertex in finite differential method
///
/// This class contains the information that will be used around specified vertex
class FDVertex
{
public:
	/// @brief FDVertex is the construction method for this class
	/// 
	/// The vertex is initialized with internal id, x and y coordinates of the vertex.
	/// The id of the vertex is given explicitly when it is called.
	///
	/// @param unsigned _id
	/// @param double _x
	/// @param double _y
	/// @pre
	/// @return 
	/// @note
	FDVertex(unsigned _id, double _x, double _y): X(_x), Y(_y), id(_id) {}
	double X; ///< coordinate of x direction
	double Y; ///< coordinate of y direction
	double EastLength; ///< the edge length to the east of the current vertex
	double WestLength; ///< the edge length to the west of the current vertex
	double NorthLength; ///< the edge length to the north of the current vertex
	double SouthLength; ///< the edge length to the south of the current vertex
	FDVertex *EastVertex; ///< the pointer to east vertex
	FDVertex *WestVertex; ///< the pointer to west vertex
	FDVertex *NorthVertex; ///< the pointer to north vertex
	FDVertex *SouthVertex; ///< the pointer to south vertex
	FDElement *NortheastElem; ///< the pointer to northeast element
	FDElement *NorthwestElem; ///< the pointer to northwest element
	FDElement *SoutheastElem; ///< the pointer to southeast element
	FDElement *SouthwestElem; ///< the pointer to southwest element
	PhysProperty Phys; ///< the physical values attached to current vertex
	
	/// @brief GetInternalID returns the internal id of specified vertex
	/// 
	///
	/// 
	/// @pre
	/// @return int
	/// @note
	int GetInternalID() { return id; }
	/// @brief Distance can get the distance between two given vertices
	/// 
	/// This is a static method and can be called directly with the class name.
	/// 
	/// @param FDVertex * vertex1
	/// @param FDVertex * vertex2
	/// @pre
	/// @return double
	/// @note This is a static method.
	static double Distance(FDVertex *vertex1, FDVertex *vertex2);
protected:
	unsigned int id; ///< internal id of the vertex
};


class FDRegion;
/// @brief FDElement is the class describing the element used in finite differential method
///
/// FDElement involves the information near current element.
class FDElement
{
public:
	//temporary construction method
	/// @brief FDElement is the construction method of this class.
	/// 
	/// The element in finite differential method is initialized with internal id 
	/// and four vertices of the rectangle element.
	/// 
	/// @param unsigned int _id
	/// @param FDVertex * _swVertex
	/// @param FDVertex * _seVertex
	/// @param FDVertex * _neVertex
	/// @param FDVertex * _nwVertex
	/// @pre
	/// @return 
	/// @note the sequence of vertex in the construction of element is sw -> se -> ne -> nw, as is shown below
	///       4----------------3
	///       |                |
	///       |                |
	///       |                |
	///       1----------------2
	FDElement( unsigned int _id, FDVertex *_swVertex, FDVertex *_seVertex, FDVertex *_neVertex, FDVertex *_nwVertex )
		:id(_id), SouthwestVertex(_swVertex), SoutheastVertex(_seVertex), NortheastVertex(_neVertex), NorthwestVertex(_nwVertex)
	{
		WestLength =FDVertex::Distance(NorthwestVertex, SouthwestVertex);
		SouthLength = FDVertex::Distance(SouthwestVertex, SoutheastVertex);
		EastLength = FDVertex::Distance(NortheastVertex, SoutheastVertex);
		NorthLength =FDVertex::Distance(NorthwestVertex, NortheastVertex);
		Area = WestLength * SouthLength;
		///TODO: judge if the corresponding lengths equal to each other
		double relError = (WestLength - EastLength) / WestLength;
		SCTM_ASSERT(relError > 0.01, 3);
		relError = (SouthLength - NorthLength) / SouthLength;
		SCTM_ASSERT(relError > 0.01, 3);
	}
	
	double WestLength; ///< the length of west edge of this element
	double EastLength; ///< the length of east edge of this element
	double NorthLength; ///< the length of north edge of this element
	double SouthLength; ///< the length of south edge of this element
	double Area; ///< the area of this element, the unit is in accordance with length
	FDVertex *NorthwestVertex; ///< the northwest vertex of this element
	FDVertex *NortheastVertex; ///< the northeast vertex of this element
	FDVertex *SouthwestVertex; ///< the southwest vertex of this element
	FDVertex *SoutheastVertex; ///< the southeast vertex of this element
	FDRegion *Region; ///< the related region of this element

	/// @brief GetInternalID returns the internal id of the specified element.
	/// 
	///
	/// 
	/// @pre
	/// @return int
	/// @note
	int GetInternalID() { return id; }

	/// @brief SetRegion sets the region of this element.
	/// 
	///
	/// 
	/// @param FDRegion * region
	/// @pre
	/// @return void
	/// @note
	void SetRegion( FDRegion *region ) {this->Region = region;}
protected:
	unsigned int id; ///< the internal id of the element
};

/// @brief FDRegion describes the region used in finite differential method
///
/// This class contains the information of the region.
class FDRegion
{
public:
	/// @brief The type of the region 
	enum RegionType
	{
		Tunneling, ///< tunneling oxide
		Trapping, ///< trapping layer
		Blocking ///< blocking oxide
	};
	/// @brief FDRegion is the construction method of the class
	/// 
	/// FDRegion is constructed with the internal id and type of the region.
	/// 
	/// @param unsigned int _id
	/// @param RegionType _type
	/// @pre
	/// @return 
	/// @note
	FDRegion(unsigned int _id, RegionType _type)
		:id(_id), Type(_type) {}

	RegionType Type; ///< type of the region, in enum RegionType
	Material * RegionMaterial; ///< the material of current region, a pointer to const material

	/// @brief AddElement adds element in current region
	/// 
	///
	/// 
	/// @param FDElement * elem
	/// @pre
	/// @return void
	/// @note
	void AddElement( FDElement *elem ) {elements.push_back(elem);}
protected:
	std::vector<FDElement *> elements; ///< the elements that are involved in current region
	unsigned int id; ///< internal id
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
#endif // !_DOMAINDETAIL_H_