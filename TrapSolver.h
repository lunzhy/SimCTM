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
	enum TrapType
	{
		eTrap,
		hTrap,
	};
	TrapSolver(FDDomain *_domain, TrapType _traptype);
	void SolveTrap();
	void UpdateTrapped();
protected:
	void eSolveTrap();
	void hSolveTrap();

	void initializeSolver();
	void refreshSolver();
	void setSolverTrapping();
	void setSolverDetrapping_BasicSRH();
	void setSolverBandToTrap();
	void setSolverTrapToBand();
	void setSolverPooleFrenkel_Frequency();
	void solveEachVertex();

protected:
	TrapType trapType;
	double temperature;
	FDDomain *domain;
	vector<FDVertex *> &vertices;

	VertexMapDouble mapTrapDensity;
	VertexMapDouble coeffMap;
	VertexMapDouble rhsMap;
	VertexMapDouble mapTrappedSolved;
};

class HoleTrapSolver
{
public:
	HoleTrapSolver(FDDomain *_domain);
};

class HoleConserveTrapSolver
{
public:
	HoleConserveTrapSolver(FDDomain *_domain);
	void SolveTrap();
	void UpdateTrapped();

protected:
	void refreshSolver();
	void setSolverTrap();
	void solveEachVertex();

protected:
	double temperature;
	double timestep;
	FDDomain *domain;
	vector<FDVertex *> &vertices;

	VertexMapDouble maphTrapDens;

	VertexMapDouble maphDensLasttime;
	VertexMapDouble maphTrappedLasttime;

	VertexMapDouble maphTrappedSolved;

	VertexMapDouble mapQuadEqu_A;
	VertexMapDouble mapQuadEqu_B;
	VertexMapDouble mapQuadEqu_C;
};

#endif
