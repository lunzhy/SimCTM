/**
* @file FDDomain.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-2   11:16
* @note
* @todo
*/

#include <vector>
#include "DomainDetails.h"
#include <iostream>

class FDDomain
{
public:
	virtual void BuildDomain() = 0;
protected:
	std::vector<FDVertex *> vertices;
	std::vector<FDElement *> elements;
	std::vector<FDRegion *> regions;

	FDVertex * getVertex(unsigned int id);
	FDElement * getElement(unsigned int id);
	FDRegion * getRegion(unsigned int id);
};

class FDDomainHelper
{
public:
	FDDomainHelper( int _cntX, int _cntY )
		:cntX(_cntX), cntY(_cntY) {}
	int IdAt(int x, int y) { return cntX * y + x; }
	int GetX(int id) { return id % cntX; }
	int GetY(int id) { return id / cntX; }
	int GetMaxX() { return cntX - 1; }
	int GetMaxY() { return cntY - 1; }
protected:
	int cntX;
	int cntY;
};
