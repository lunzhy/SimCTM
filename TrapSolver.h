/**
* @file TrapSolver.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-12-6   15:40
* @note
* @todo
*/
#ifndef _TRAPSOLVER_H_
#define _TRAPSOLVER_H_

#include <map>
#include <vector>

class FDDomain;
class FDVertex;
using std::vector;
typedef std::map<int, double> VertexMapDouble; // <vertID, value>, map with vertex index as the key

class TrapSolver
{
public:
	TrapSolver(FDDomain *_domain);
	void SolveTrap();
	void UpdateTrapped();
protected:
	void initializeSolver();
	void refreshSolver();
	void setSolverTrapping();
	void setSolverDetrapping_BasicSRH();
	void setSolverBandToTrap();
	void setSolverTrapToBand();
	void solveEachVertex();

protected:
	double temperature;
	FDDomain *domain;
	vector<FDVertex *> &vertices;

	VertexMapDouble eXsectionMap;
	VertexMapDouble eMobilityMap;
	VertexMapDouble eTrapDensMap;
	VertexMapDouble coeffMap;
	VertexMapDouble rhsMap;
	VertexMapDouble eTrappedMap;
};

#endif
