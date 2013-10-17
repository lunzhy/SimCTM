/**
* @file PoissonSolver.h
* @brief This file defines all the necessary classes and methods for solving two-dimensional Poisson problem.
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
typedef std::map<int, int, std::less<int>> VertexMapInt; // the map used for get the equation index with given vertex internal id

class FDDomain;

/// @brief TwoDimPoissonSolver is a Poisson solver used to solve two-dimensional Poisson problem
///
/// The solvers utilizes finite-difference method to solve two-dimensional Poisson problem.
/// Currently this class is inherited from the sparse matrix solver. But in other conditions, the solver serves
/// as a member of the specific class.
class TwoDimPoissonSolver : public SctmSparseMatrixSolver
{
public:
	friend class SctmUtils::SctmDebug; // so SctmDebug can get the private and protected member of this solver.
	/// @brief TwoDimPoissonSolver is the construction method of this solver.
	/// 
	/// This Poisson solver reads in the domain to be solved. Internally, all the vertices concluded in the domain
	/// is read in the Poisson solver. The Poisson solver is often used within the whole domain of the structure.
	/// The solver is prepared in the construction method.
	/// 
	/// @param FDDomain * domain
	/// @pre
	/// @return 
	/// @note
	TwoDimPoissonSolver(FDDomain *domain);
	/// @brief SolvePotential is called to solve the specific problem.
	/// 
	/// This method is the only method that is called outside this class.
	/// 
	/// @pre
	/// @return void
	/// @note
	void SolvePotential();
private:
	vector<FDVertex *> &vertices; ///< the vertices of the domain to be solved
	vector<double> potential; ///< the potential of the vertices with same sequence in the vertices vector (same sequence of the equations)
	vector<double> rhsVector; ///< the right-hand side of the equations
	VertexMapInt equationMap; ///< map the equation index to the internal id of the vertices. (equation index is the same with the index of vertices vector)
protected:
	/// @brief prepareSolver is used to prepare the Poisson solver.
	/// 
	/// Several preparation is finished in this method. Vertex map is build. The coefficient matrix is filled in and refreshed in consideration of the
	/// boundary conditions. The right-hand side vector is filled in with initial consideration of the boundary conditions.
	/// 
	/// @pre
	/// @return void
	/// @note
	void prepareSolver();
	/// @brief buildVertexMap is used to build the vertex map.
	/// 
	/// The vertex map is used to get the equation index, also the index of the specific vertex in vertices list of this class, with given
	/// vertex internal id. This is frequently used in filling the coefficient matrix.
	/// 
	/// @pre
	/// @return void
	/// @note
	void buildVertexMap();
	/// @brief buildCoefficientMatrix is used to build coefficient matrix without consideration of the boundary conditions.
	/// 
	/// This methods deal with the sparse matrix member inherited from the sparse matrix solver. This is the initial fill-in of
	/// the coefficient matrix.
	/// 
	/// @pre
	/// @return void
	/// @note
	void buildCoefficientMatrix();
	/// @brief buildRhsVector is used to build the right-hand side vector without consideration of the boundary conditions.
	/// 
	/// The right-hand side vector is build with the total net charge in the vertex.
	/// 
	/// @pre
	/// @return void
	/// @note
	void buildRhsVector();
	/// @brief refreshCoefficientMatrix is used to refresh the coefficient matrix with the boundary conditions.
	/// 
	/// Only equation of the vertex with Dirichlet BC (first type BC) is revised. In terms of the other types of BC, only
	/// right-hand side is revised.
	/// 
	/// @pre
	/// @return void
	/// @note
	void refreshCoefficientMatrix();
	/// @brief refreshRhs is used to refresh the right-hand side vector of the matrix equation.
	/// 
	/// The method to refresh right-hand side of the equation with BC_Neumann and BC_Artificial is a little complicated.
	/// When considering one element of the specific vertex, the revision of the right-hand side has to be done with the consideration
	/// of two adjacent elements of this element. This is because the effect of electric field in BC_Neumann in both directions depends
	/// on the existence of the above-mentioned adjacent elements.
	/// 
	/// @pre
	/// @return void
	/// @note
	void refreshRhs();
	/// @brief fillBackPotential is used fill the the potential property of the vertex with obtained calculation result.
	/// 
	///
	/// 
	/// @pre
	/// @return void
	/// @note
	void fillBackPotential();
};

#endif