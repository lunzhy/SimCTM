/**
* @file FDDomain.h
* @brief This file deals with the domain definition with finite differential method
*
* FDDomain declares the finite differential domain used in the simulation. In addition,
* a helper class for FDDomian is declares, which is used in initialization of  rectangle 
* domain to position each point in the inner iteration with given structures 
*
* @author
* @version 
* @date 2013-7-2   11:16
* @note
* @todo
*/
#pragma once
#include <vector>
#include "DomainDetails.h"
#include <iostream>


/// @brief FDDomain describes the simulation domain in finite differential method.
///
/// This class is used as a base class for FD domain. It only contains the general properties in FD domain.
class FDDomain
{
public:

	/// @brief BuildDomain builds the specified domain structures, setting vertices, elements and regions.
	/// 
	/// This class in a virtual class. The detailed implementation of this class, i.e. detailed information
	/// of the target structures and method to build the domain, is defined in the derived class. 
	/// 
	/// @pre
	/// @return void
	/// @note
	virtual void BuildDomain() = 0;
protected:
	std::vector<FDVertex *> vertices; ///< the vertices contained in the domain
	std::vector<FDElement *> elements; ///< the elements contained in the domain
	std::vector<FDRegion *> regions; ///< the regions contained in the domain

	void setVertexPhysics();

public:
	/// @brief getVertex can get the vertex object with given id
	/// 
	/// This method returns the pointer of specified vertex object. In practice, a pointer with same type is
	/// always ready to get the returned pointer and do the following process.
	/// 
	/// @param unsigned int id
	/// @pre
	/// @return FDVertex *
	/// @note
	FDVertex * getVertex(unsigned int id);
	/// @brief getElement can get the element object with given id
	/// 
	/// This method returns the pointer of specified element object. In practice, a pointer with same type is
	/// always ready to get the returned pointer and do the following process.
	/// @param unsigned int id
	/// @pre
	/// @return FDElement *
	/// @note
	FDElement * getElement(unsigned int id);
	/// @brief getRegion can get the region object with given id
	/// 
	/// This method returns the pointer of specified region object. In practice, a pointer with same type is
	/// always ready to get the returned pointer and do the following process.
	/// @param unsigned int id
	/// @pre
	/// @return FDRegion *
	/// @note
	FDRegion * getRegion(unsigned int id);
};


/// @brief FDDomianHelper is a helper class in building rectangle structure.
///
/// When initializing a rectangle structure, this helper class is used to deal with the positioning of a vertex/grid.
/// It provide transformation between inner iteration id and vertex/grid position in coordinates (X, Y). 
/// In addition, the inner id sequence is shown as follow (an exmaple of a 6*6 vertex mesh)
///
///                    30..31..32..33..34..35
///                    24..25..26..27..28..29
///                    18..19..20..21..22..23
///                    12..13..14..15..16..17
///                    6...7...8...9...10..11
///                    0...1...2...3...4...5
///
class FDDomainHelper
{
public:
	/// @brief FDDomainHelper is the constructor method for this class
	/// 
	/// FDDomainHelper must be constructed with give total number in X and direction
	/// 
	/// @param int _cntX
	/// @param int _cntY
	/// @pre
	/// @return 
	/// @note
	FDDomainHelper( int _cntX, int _cntY )
		:cntX(_cntX), cntY(_cntY) {}
	/// @brief IdAt returns the inner id number with given X and Y.
	/// 
	///
	/// 
	/// @param int x
	/// @param int y
	/// @pre
	/// @return int
	/// @note
	int IdAt(int x, int y) { return cntX * y + x; }
	/// @brief GetX returns the X coordinate with given inner id.
	/// 
	///
	/// 
	/// @param int id
	/// @pre
	/// @return int
	/// @note
	int GetX(int id) { return id % cntX; }
	/// @brief GetY returns the Y coordinate with given inner id.
	/// 
	///
	/// 
	/// @param int id
	/// @pre
	/// @return int
	/// @note
	int GetY(int id) { return id / cntX; }
	/// @brief GetMaxX returns the maximum X coordinate
	/// 
	///
	/// 
	/// @pre
	/// @return int
	/// @note
	int GetMaxX() { return cntX - 1; }
	/// @brief GetMaxY returns the maximum Y coordinate
	/// 
	///
	/// 
	/// @pre
	/// @return int
	/// @note
	int GetMaxY() { return cntY - 1; }
protected:
	int cntX; ///< the number/count of points in X direction
	int cntY; ///< the number/count of points in Y direction
};
