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


using std::vector;
class FDVertex;
class TwoDimPoisson : public SctmMath::SparseMatrixSolver
{
public:
	TwoDimPoisson(vector<FDVertex *> &_vertices):vertices(_vertices)
	{}
private:
	vector<FDVertex *> &vertices;
	vector<double> potential;
protected:
	void buildCoefficientMatrix();
	void buildRHS();
	void parseBoundaryCondition();
	void refreshCoefficientMatrix();
	void refreshRHS();
	void solvePotential();
};