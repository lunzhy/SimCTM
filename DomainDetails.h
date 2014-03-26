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
#include "SctmMath.h"
#include <map>

namespace MaterialDB
{
	class Material;
}
namespace SctmPhys
{
	class PhysProperty;
	class TrapProperty;
}
namespace SctmUtils
{
	class SctmData;
}

using MaterialDB::Material;
using std::string;
using std::map;
using SctmMath::VectorValue;
using SctmPhys::PhysProperty;
using SctmPhys::TrapProperty;
/// @brief FDBoundary is a the class describing the boundary conditions in finite differential method
///
/// The object of FDBoundary is a member in FDVertex to store the boundary information(valid/invalid).
/// FDBoundary is determined by the location of the vertex, either at the contact (governed by FDContanct) or 
/// at the boundary of the domain (artificial boundary)
class FDBoundary
{
public:
	/// @brief BCName is the enum of the name of the boundary conditions
	enum BCName
	{
		Potential, ///< potential boundary condition, note that BC_Neumann and BC_Artificial for potential is electric field.
		eDensity ///< electron density boundary condition. The BC_Cauchy of eDensity describes the electron current at the boundary.
	};

	/// @brief BndCond is the enum of different kind of boundary conditions.
	enum BCType
	{
		BC_Dirichlet, ///< the first type of boundary condition, giving the value to the boundary.
		BC_Neumann, ///< the second type of boundary condition, giving the value to the normal derivative at the boundary.
		BC_Cauchy, ///< the third type of boundary condition, giving the value to the normal derivative and the variable itself.
		///[Notice!]
		///It should be noticed that the specific form of BC_Neumnn or BC_Cauchy, in terms of the coefficient of the derivative,
		///is not defined here. It is determined by the physical expression that defines the boundary condition, which is an additional
		///information of the boundary condition.
		BC_Artificial ///< a special Neumann boundary condition, used when the vertex is at the artificial boundary of the domain.
	};

	/// @brief CurrentTag is the token to distinguish the tunneling-in and tunneling-out boundary
	enum TunnelTag
	{
		noTunnel,
		eTunnelOut,
		eTunnelIn,
		hTunnelIn,
		hTunnelOut
	};
	/// @brief FDBoundary is the construction method of this class
	/// 
	/// The object of BndCond is constructed with the construction of the specified vertex, because it is a member of the 
	/// vertex object. So the validity of the boundary condition is set to false, meaning that by default, the vertex is not
	/// at the boundary of the domain. 
	/// By default, the calling of non-existed map index will get the return value of 0 (false).
	/// So the construction method without parameters is needed.
	/// 
	/// @pre
	/// @return 
	/// @note
	FDBoundary();
	/// @brief RefreshBndCondValue is called to refresh the boundary condition value with given boundary condition name.
	/// 
	/// Both values of the boundary condition must be given in this method, otherwise the value of 0 is applied 
	/// in the refreshing.
	/// 
	/// @param BCName bcName
	/// @param double bcValue1
	/// @param double bcValue2
	/// @pre
	/// @return void
	/// @note
	void RefreshBndCond(BCName bcName, VectorValue bcNormVec);
	void RefreshBndCond(BCName bcName, double newValue, VectorValue bcNormVec = VectorValue(0, 0));
	void RefreshBndCond(BCName bcName, BCType bcType, double bcValue = 0, VectorValue bcNormVec = VectorValue(0, 0));
	/// @brief SetBnd is used to set the boundary type and its direction.
	/// 
	/// In this method, the boundary condition direction is set with the same direction as the boundary direction
	/// In some cases, for example, the current density boundary condition, the bc direction may be different with
	/// the boundary direction.
	/// 
	/// @param BCName bcName
	/// @param BCType bcType
	/// @param VectorValue bndVec
	/// @param double bcValue
	/// @pre
	/// @return void
	/// @note
	void SetBnd(BCName bcName, BCType bcType, VectorValue bndVec, double bcValue = 0);
	void SetTunnelTag(TunnelTag tag);
	/// @brief Valid is used to return the validity of the boundary condition with given specified BC name.
	/// 
	/// Both non-existent boundary condition and boundary condition with false validity will return false.
	/// 
	/// @param BCName bcName
	/// @pre
	/// @return bool
	/// @note
	bool Valid(BCName bcName);
	/// @brief GetBCType is used to obtain the boundary condition type of given name of boundary condition.
	/// 
	/// 
	/// @param BCName bcName
	/// @pre
	/// @return FDBoundary::BCType
	/// @note
	BCType GetBCType(BCName bcName);
	/// @brief GetBCValue is used to obtain the boundary condition value of given name of boundary condition.
	/// 
	/// This method is used to get the value of BC_Dirichlet boundary condition. This value is stored in
	/// the map bc_values.
	/// 
	/// @param BCName bcName
	/// @pre
	/// @return double
	/// @note
	double GetBCValue(BCName bcName);
	TunnelTag GetBCTunnelTag();
	VectorValue &GetBndDirection(BCName bcName);
	VectorValue &GetBCNormVector(BCName bcName);
protected:
	//bool valid; ///< the validity of the boundary condition. It is a token to indicate a boundary vertex
	map<BCName, bool> bnd_valid; ///< the validity of the boundary condition with given boundary condition name
	//Finding the value of non-existed key will return false.
	map<BCName, VectorValue> bnd_normVec; ///< the map to store normal vector of the boundary (not boundary condition)
	//boundary direction is used to determine whether the corresponding vertex exists.
	map<BCName, BCType> bc_types; ///< the map to store the types of different boundary conditions
	map<BCName, double> bc_values; ///< the map to store the values of different boundary conditions.
	//When BC is a vector value, if the BC has the same direction with the normal vector, this value is positive. reversed direction, negative 
	map<BCName, VectorValue> bc_normVec; ///< the map to store normal vector of the boundary condition
	TunnelTag tunTag; ///< the map to store the tunneling tag for the boundary
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
	FDVertex(unsigned _id, double _x, double _y);
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

	PhysProperty *Phys; ///< the physical values attached to current vertex
	TrapProperty *Trap; ///< the trap property attached to current vertex
	FDBoundary BndCond; ///< the boundary condition of current vertex
	
	/// @brief IsAtBoundary is used to check if the vertex is a boundary vertex with specified boundary name
	/// 
	///
	/// 
	/// @param FDBoundary::BCName bcName
	/// @pre
	/// @return bool
	/// @note
	bool IsAtBoundary(FDBoundary::BCName bcName);
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
	int GetID() const { return id; }
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
	FDElement( unsigned int _id, FDVertex *_swVertex, FDVertex *_seVertex, FDVertex *_neVertex, FDVertex *_nwVertex );
		
	
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
	int GetID() { return id; }

	/// @brief SetRegion sets the region of this element.
	/// 
	///
	/// 
	/// @param FDRegion * region
	/// @pre
	/// @return void
	/// @note
	void SetRegion( FDRegion *region ) {this->Region = region;}
	void SetVertexAdjacent();
protected:
	unsigned int id; ///< the internal id of the element
};

/// @brief FDRegion describes the region used in finite differential method
///
/// This class contains the information of the region.
class FDRegion
{
	friend class TripleCells;
	friend class SctmUtils::SctmData;
public:
	string RegName;
	/// @brief FDRegion is the construction method of the class
	/// 
	/// FDRegion is constructed with the internal id and type of the region.
	/// 
	/// @param unsigned int _id
	/// @param RegionType _type
	/// @pre
	/// @return 
	/// @note
	FDRegion(unsigned int _id, string _name, Material *_mat)
		:id(_id), RegName(_name), Mat(_mat) {}

	Material *Mat; ///< the material of current region, a pointer to const material

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
	double Workfunction;
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
	FDContact(unsigned int _id, string _name, double _voltage)
		: id(_id), ContactName(_name), Voltage(_voltage){}
	FDContact(unsigned int _id, string _name, double _voltage, double _workfunction)
		: id(_id), ContactName(_name), Voltage(_voltage), Workfunction(_workfunction){}
	/// @brief AddVertex is used to add the vertex that belongs to the specific contact.
	/// 
	///
	/// 
	/// @param FDVertex * vert
	/// @pre
	/// @return void
	/// @note
	void AddVertex(FDVertex *vert) { vertices.push_back(vert); }
	/// @brief GetContactVerts is used to get all vertices in the contact
	/// 
	/// This method is used in one dimensional problem to get or set the properties
	/// of the vertex at specific contact.
	/// 
	/// @pre
	/// @return std::vector<FDVertex *> &
	/// @note TODO: this could be enhanced.
	std::vector<FDVertex *> & GetContactVerts();
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