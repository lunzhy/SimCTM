/**
* @file FDDomainTest.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-7-5   11:15
* @note
* @todo
*/

#include "FDDomain.h"
class FDRegion;
class FDDomainTest : FDDomain
{
public:
	void BuildDomain();
protected:
	double xLength;
	double yLengthTunnel;
	double yLengthTrap;
	double yLengthBlock;

	int xCntVertex;
	int yCntVertexTunnel;
	int yCntVertexTrap;
	int yCntVertexBlock;
	int yCntTotalVertex;

	double xGrid;
	double yGridTunnel;
	double yGridTrap;
	double yGridBlock;

private:
	void setParameters();
	void setDomainDetails();
	void setAdjacency();
	void prepareStructures();
	void printStructure();
	
	double yNextGridLength( int iy );
	double xNextGridLength( int ix );
	FDRegion * thisRegion( int elemY);
};