/**
* @file DDSolver.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-29   10:57
* @note
* @todo
*/

#ifndef _DDSOLVER_H_
#define _DDSOLVER_H_

#include <vector>
#include <map>
#include "SctmUtils.h"
#include "MatrixSolver.h"

using std::vector;
using SctmMath::SctmSparseMatrixSolver;

class FDDomain;
class DriftDiffusionSolver
{
	//friend class SctmUtils::SctmDebug;
	typedef std::map<int, int> MapForVertex; // <equationID, vertID>
	typedef std::map<int, double> MapForPrpty; // <vertID, property value>
public:
	DriftDiffusionSolver(FDDomain *domain);
private:
	vector<FDVertex *> vertices;
	vector<FDVertex *> &totalVertices;

	double temperature;
	double timeStep;
	//the material and physical properties
	MapForPrpty mobilityMap; // mobility is used, so diffusion coefficient is derived
	MapForPrpty potentialMap;
	MapForPrpty lastDensityMap; // the density of last time step
	MapForVertex vertMap;

	vector<double> rhsVector;
	vector<double> eDensity;

	SctmSparseMatrixSolver matrixSolver;
protected:
	void prepareSolver();
	void getDDVertices(FDDomain *domain);
	void buildVertexMap();
	void setBndCondCurrent(vector<double> &current);
	void buildCoefficientMatrix();
	void buildRhsVector();
	void refreshCoefficientMatrix();
	void refreshRhs();
	void setTimeStep();
};

class DDTest : public DriftDiffusionSolver
{
public:
	DDTest(FDDomain *_domain);
};
#endif