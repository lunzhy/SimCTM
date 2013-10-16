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
	void SolveDD();
protected:
	vector<FDVertex *> vertices;
	vector<FDVertex *> &totalVertices;

	double temperature;
	double timeStep;
	//the material and physical properties
	MapForPrpty mobilityMap; // mobility is used, so diffusion coefficient is derived
	MapForPrpty potentialMap;
	MapForPrpty lastElecDensMap; // the electron density of last time step
	MapForVertex vertMap;

	vector<double> rhsVector;
	vector<double> elecDensity;

	SctmSparseMatrixSolver matrixSolver;
protected:
	void prepareSolver();
	void getDDVertices(FDDomain *domain);
	void buildVertexMap();
	/// @brief setBndCondCurrent is used to set the current boundary condition of the drift-diffusion area
	/// 
	/// Currently, this method is restricted in the simple structures.
	/// The boundary condition current consists the tunneling-in current from tunneling oxide and the tunneling-out current
	/// to the flow through the blocking layer.
	/// 
	/// @param vector<double> & in_current
	/// @pre
	/// @return void
	/// @note
	void setBndCondCurrent(vector<double> &in_current);
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
	void SolveDD();
protected:
	void prepareSolver(); // use method with the same name with base class
	void buildVertexMap();
	void setBndCondCurrent();
};
#endif