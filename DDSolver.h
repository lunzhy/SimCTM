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
#include "MatrixSolver.h"

class FDDomain;
class FDVertex;

using std::vector;
using SctmMath::SctmSparseMatrixSolver;

class DriftDiffusionSolver
{
	//friend class SctmUtils::SctmDebug;
	typedef std::map<int, int> VertexMapInt; // <VertID, equationID> store map for equation ID
	typedef std::map<int, double> VertexMapDouble; // <vertID, property value>, store the map for physical property
public:
	enum BCMethod
	{
		DirectDiscretization,
		UsingCurrentDensity,
	};
	DriftDiffusionSolver(FDDomain *domain);
	virtual void SolveDD();
	double CalculateTotalLineDensity();
protected:
	BCMethod bcMethod;
	bool useCrankNicolsonMethod;

	vector<FDVertex *> vertices;
	vector<FDVertex *> &totalVertices;

	double temperature;
	double timeStep;
	double lastTimeStep; // to revert the coefficients because they have been modified in last time step 
	//the material and physical properties
	VertexMapDouble mobilityMap; ///< mobility is used, so diffusion coefficient is derived
	VertexMapDouble potentialMap;
	VertexMapDouble lastElecDensMap; ///< the electron density of last time step
	VertexMapInt equationMap;

	//the index of this kind of vector is equation id, which is obtained from equation map.
	vector<double> rhsVector;
	vector<double> elecDensity;

	SctmSparseMatrixSolver matrixSolver;
protected:
	void initializeSolver();
	void prepareSolver();
	void getDDVertices(FDDomain *domain);
	virtual void buildVertexMap();
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
	/// @brief buildCoefficientMatrix
	/// 
	/// The filling in of the coefficient matrix utilizes the potential of vertex at last time step, which make the equation
	/// explicit in potential and velocity and implicit in carrier density.
	/// IMPORTANT Another important notice is that this method only fills the main part of the coefficient matrix. Time dependent addition
	/// of the coefficient is not done here.
	/// 
	/// @pre
	/// @return void
	/// @note
	void buildCoefficientMatrix();
	void setCoefficientBCVertex_UsingCurrent(FDVertex *vert);
	void setCoefficientBCVertex_DirectDiscretization(FDVertex *vert);
	void setCoefficientInnerVertex(FDVertex *vert);
	/// @brief refreshCoefficientMatrix
	/// 
	/// The boundary conditions are always BC_Cauchy, so there is no need to refresh the matrix for this reason.
	/// However, due to that the time step is likely to be different in each simulation step, the time dependent item
	/// in the coefficient matrix has to be add to the corresponding item.
	/// 
	/// @pre
	/// @return void
	/// @note
	void refreshCoefficientMatrix();
	/// @brief buildRhsVector
	/// 
	/// This method is called at each simulation step. Because in each time step, the density of each vertex has to be
	/// refreshed. It is called before refreshRhsWithBC .
	/// 
	/// @pre
	/// @return void
	/// @note
	void buildRhsVector();
	void buildRhsVector_UsingCurrent();
	double getRhsBCVertex_DirectDiscretiztion(FDVertex *vert);
	double getRhsBCVertex_UsingCurrent(FDVertex *vert);
	double getRhsInnerVertex(FDVertex *vert);
	void setTimeStep();
	void fillBackElecDens();
	virtual void refreshBoundary();
	void getDeltaXYAtVertex(FDVertex *vert, double &dx, double &dy);
};

class DDTest : public DriftDiffusionSolver
{
public:
	DDTest(FDDomain *_domain);
	void SolveDD();
protected:
	void buildVertexMap();
	void refreshBoundary();
	void setBndCurrent();
	void setBndDensity();
	void refreshCoeffMatrixDueToBC();
};
#endif