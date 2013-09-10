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
#include "FDDomain.h"
#include <vector>
#include <map>

using namespace SctmMath;
using std::vector;
typedef std::map<int, int, std::less<int>> MapForVertex;

class FDDomain;
class TwoDimPoissonSolver : public SctmSparseMatrixSolver
{
public:
	friend class SctmUtils::SctmDebug;
	TwoDimPoissonSolver(FDDomain *domain);
	void SolvePotential();
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
	void fillBackPotential();
};

#endif