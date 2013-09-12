/**
* @file MatrixSolver.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-14   19:14
* @note
* @todo
*/

#include "MatrixSolver.h"
#include "SctmUtils.h"

namespace SctmMath
{
	void SctmSparseMatrixSolver::SolveMatrix(std::vector<double> &rhs, std::vector<double> &solution)
	{
		//This is very IMPORTANT
		this->matrix.makeCompressed();
		
		//TODO: check the sparse matrix
		SCTM_ASSERT(rhs.size() == solution.size(), 10005);
		int vectorSize = rhs.size();
		//Map is used because internal type of Eigen is required when solving the matrix problem.
		Eigen::Map<Eigen::VectorXd> rhsOfEigen(rhs.data(), vectorSize, 1);
		Eigen::Map<Eigen::VectorXd> solutionOfEigen(solution.data(), vectorSize, 1);

		//use SparseLU to solve sparse matrix problem. SparseLU supports general square sparse systems
		Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> sparseSolver; //or int
		sparseSolver.analyzePattern(this->matrix);
		sparseSolver.factorize(this->matrix);
		solutionOfEigen = sparseSolver.solve(rhsOfEigen);

		if (sparseSolver.info() != Eigen::Success)
		{
			SCTM_ASSERT(SCTM_ERROR, 10006);
		}
	}
}

