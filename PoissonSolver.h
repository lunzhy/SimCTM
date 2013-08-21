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
#ifndef _POISSONSOLVER_H_
#define _POISSONSOLVER_H_
#include "MatrixSolver.h"
#include <vector>
#include <map>

using namespace SctmMath;
using std::vector;
class FDVertex;
typedef std::map<int, int, std::less<int>> MapForVertex;

class TwoDimPoisson : public SctmSparseMatrixSolver
{
public:
	TwoDimPoisson(vector<FDVertex *> &_vertices):vertices(_vertices)
	{}
private:
	vector<FDVertex *> &vertices;
	vector<double> potential;
	vector<double> rhsVector;
	MapForVertex vertMap;
protected:
	void prepareSolver();
	void buildVertexMap();
	void buildCoefficientMatrix();
	void buildRhsVector();
	void refreshCoefficientMatrix();
	void refreshRHS();
	void solvePotential();
	void fillBackPotential();
};

#endif