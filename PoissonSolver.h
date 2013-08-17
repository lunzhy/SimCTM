/**
* @file PoissonSolver.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-14   20:09
* @note
* @todo
*/
#pragma once
#include "MatrixSolver.h"
#include <vector>
#include <map>


using std::vector;
class FDVertex;
typedef std::map<int, int, std::less<int>> MapForVertex;


class TwoDimPoisson : public SctmMath::SparseMatrixSolver
{
public:
	TwoDimPoisson(vector<FDVertex *> &_vertices):vertices(_vertices)
	{}
private:
	vector<FDVertex *> &vertices;
	vector<double> potential;
	MapForVertex vertMap;
protected:
	void prepareVertexMap();
	void buildCoefficientMatrix();
	void buildRHS();
	void parseBoundaryCondition();
	void refreshCoefficientMatrix();
	void refreshRHS();
	void solvePotential();
};