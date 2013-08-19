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
using std::string;


/// @brief FDBoundary is a the class describing the boundary conditions in finite differential method
///
/// The object of FDBoundary is a member in FDVertex to store the boundary information(valid/invalid).
/// FDBoundary is determined by the location of the vertex, either at the contact (governed by FDContanct) or 
/// at the boundary of the domain (artificial boundary)
class FDBoundary
{
public:
	
	/// @brief BndCond is the enum of different kind of boundary conditions
	enum BndCond
	{
		BC_Dirichlet, ///< the first type of boundary condition
		BC_Neumann, ///< the second type of boundary condition
		BC_Artificial, ///< a special Neumann boundary condition, used when the vertex is at the boundary of the domain
	};
	/// @brief FDBoundary is the construction method of this class
	/// 
	/// The object of BndCond is constructed with the construction of the specified vertex, because it is a member of the 
	/// vertex object. So the validity of the boundary condition is set to false, meaning that by default, the vertex is not
	/// at the boundary of the domain. 
	/// 
	/// @pre
	/// @return 
	/// @note
	FDBoundary():valid(false){}
	/// @brief SetBndCond is called to set the boundary condition of the specific vertex
	/// 
	/// The boundary condition is set with boundary type and value. 
	/// With regard to BC_Dirichlet, only the boundary value 1 is used.
	/// With regard to BC_Neumann, both boundary value 1 and 2 are used.
	/// With respect to the BC_Artificial, no boundary value need to be set.
	/// 
	/// @param BndCond bndtype
	/// @param double bndvalue1
	/// @param double bndvalue2
	/// @pre
	/// @return void
	/// @note
	void SetBndCond(BndCond bndtype, double bndvalue1, double bndvalue2);
	/// @brief Valid is used to return the validity of the boundary condition
	/// 
	/// If the validity of the boundary condition is true, it means that the vertex is indeed a boundary vertex.
	/// 
	/// @pre
	/// @return bool
	/// @note It is important to consider the direction when setting the boundary condition
	bool Valid() const { return valid; }
	/// @brief BndType returns the boundary type of the vertex
	/// 
	/// Because each of vertices has a property of the boundary condition, so this method is use to check if the
	/// vertex is at the boundary.
	/// 
	/// @pre
	/// @return FDBoundary::BndCond
	/// @note if the vertex is not at boundary, this method also returns value 0.
	BndCond BndType() const { return bndType; }
	/// @brief BndValue2 is used to obtain the value of the boundary condition.
	/// 
	/// When the boundary condition of the vertex is BC_Dirichlet, the value is the potential of the vertex.
	/// When the boundary condition of the vertex is BC_Neumnn, the value is the electric field in x direction.
	/// 
	/// @pre
	/// @return double
	/// @note The value is normalized.
	double BndValue1() const;
	/// @brief BndValue2 is used to obtain the value of electric field in y direction
	/// 
	/// BndValue2 is valid only in the case of BC_Neumann. It represents the electric field.
	/// 
	/// @pre
	/// @return double
	/// @note
	double BndValue2() const;
	/// @brief BndValuePotential is used to obtain the potential value of the boundary condition.
	/// 
	/// This method can be called only in the condition of BC_Dirichlet.
	/// 
	/// @pre
	/// @return double
	/// @note
	double BndValuePotential() const;
	/// @brief BndValueElecFieldWestEast is used to obtain the electric field value of the boundary condition.
	/// 
	/// This method can be only called in the condition of BC_Neumann.
	/// This method returns the value of electric field in X direction, i.e. the direction from west to east.
	/// 
	/// @pre
	/// @return double
	/// @note
	double BndValueElecFieldWestEast() const;
	/// @brief BndValueElecFieldSouthNorth is used to obtain the electric field value of the boundary condition.
	/// 
	/// This method can be only called in the condition of BC_Neumann.
	/// This method returns the value of electric field in Y direction, i.e. the direction from south to north
	/// 
	/// @pre
	/// @return double
	/// @note
	double BndValueElecFieldSouthNorth() const;
protected:
	BndCond bndType; ///< the type of the boundary condition
	double bndValue1; ///< the value of the  boundary condition, represents potential or electric field in direction from west to east (normalized)
	double bndValue2; ///< the electric field in direction from south to north, needed with respect to the BC_Neumann
	bool valid; ///< the validity of the boundary condition. It is a token to indicate a boundary vertex.
	
};



class FDElement;
class FDContact;
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
	FDVertex(unsigned _id, double _x, double _y): X(_x), Y(_y), id(_id)
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
	}
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
	FDContact *Contact; ///< the pointer to the contact the vertex belongs

	PhysProperty Phys; ///< the physical values attached to current vertex
	FDBoundary BoundaryCond; ///< the boundary condition of current vertex
	
	/// @brief IsAtBoundary is used to check if the vertex is a boundary vertex
	/// 
	///
	/// 
	/// @pre
	/// @return bool
	/// @note
	bool IsAtBoundary();
	/// @brief IsAtContact is used to check if the vertex belongs to a contact.
	/// 
	///
	/// 
	/// @pre
	/// @return bool
	/// @note
	bool IsAtContact() const { return (Contact!=NULL); }
	/// @brief GetInternalID returns the internal id of specified vertex
	/// 
	/// This is used in solving two dimensional Poisson equation.
	/// 
	/// @pre
	/// @return int
	/// @note
	int GetInternalID() const { return id; }
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
	/// @brief SetContact is used to set the belonging contact of this vertex
	/// 
	///
	/// 
	/// @param FDContact * contact
	/// @pre
	/// @return void
	/// @note
	void SetContact(FDContact *contact) { this->Contact = contact; }
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
	void AddElement( FDElement *elem ) { elements.push_back(elem); }
protected:
	std::vector<FDElement *> elements; ///< the elements that are involved in current region
	unsigned int id; ///< internal id
};


/// @brief FDContact is the class describing the contact in finite differential domain.
class FDContact
{
public:
	string ContactName; ///< contact name
	double Voltage; ///< the contact voltage
	/// @brief FDContact is the construction method of the class.
	/// 
	/// A FDContact is initialized with given contact internal id, contact name and the applied voltage.
	/// 
	/// @param unsigned int _id
	/// @param string _name
	/// @param double _voltage
	/// @pre
	/// @return 
	/// @note
	FDContact(unsigned int _id, string _name, double _voltage):id(_id), ContactName(_name), Voltage(_voltage){}
	/// @brief AddVertex is used to add the vertex that belongs to the specific contact.
	/// 
	///
	/// 
	/// @param FDVertex * vert
	/// @pre
	/// @return void
	/// @note
	void AddVertex(FDVertex *vert) { vertices.push_back(vert); }
protected:
	std::vector<FDVertex *> vertices; ///< the vertices belong to this contact
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