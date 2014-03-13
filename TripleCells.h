/**
* @file TripleCells.h
* @brief This is the header file for the class used to build a domain of triple cells.
*
*
*
* @author
* @version 
* @date 2014-3-3   16:05
* @note
* @todo
*/

#ifndef _TRIPLECELLS_H_
#define _TRIPLECELLS_H_
#include "FDDomain.h"
#include <string>

using std::string;

class TripleCells : public FDDomain
{
public:
	TripleCells();
protected:
	double voltageGate1;
	double workfunctionGate1;
	double voltageGate2;
	double workfunctionGate2;
	double voltageGate3;
	double workfunctionGate3;

	string matNameStrIso;
	string matNameStrTunnel;
	string matNameStrTrap;
	string matNameStrBlock;

	double widthGate1;
	double widthIso2; //no isolation 1 in SimCTM structure
	double widthGate2;
	double widthIso3;
	double widthGate3;

	int gridWidthGate1;
	int gridWidthIso2;
	int gridWidthGate2;
	int gridWidthIso3;
	int gridWidthGate3;

	double thickIso;
	double thickBlock;
	double thickTrap;
	double thickTunnel;

	int gridThickIso;
	int gridThickBlock;
	int gridThickTrap;
	int gridThickTunnel;

	double temperature;

protected:
	void buildStructure();
	void setParametersFromParamParser();
	void setDomainDetails();
	void setAdjacency();
	void setTrapDistribution();
	void postProcessOfDomain();

	//the methods below are used to calculate the coordinates of specific vertex.
	double getVertCoordX(int idX, int idY);
	double getVertCoordY(int idX, int idY);
	int getVertIdAt(int idX, int idY);
	bool isValidVertex(int idX, int idY);
	void setSingleElement(int &idElem, FDRegion *region, int xbegin, int xend, int ybegin, int yend);

};

#endif