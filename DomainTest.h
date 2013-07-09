/**
* @file DomainTest.h
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
class DomainTest : FDDomain
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

	double xGrid;
	double yGridTunnel;
	double yGridTrap;
	double yGridBlock;

private:
	void setStructures();
	void prepareStructures();
	void PrintStructure();
};