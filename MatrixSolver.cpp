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
	void SparseMatrixSolver::SolveMatrix(std::vector<double> &rhs, std::vector<double> &solution)
	{
		//TODO: check the sparse matrix

		SCTM_ASSERT(rhs.size() == solution.size(), 5);
		int vectorSize = rhs.size();
		//Map is used because internal type of Eigen is required when solving the matrix problem.
		Eigen::Map<Eigen::VectorXd> rhsOfEigen(rhs.data(), vectorSize, 1);
		Eigen::Map<Eigen::VectorXd> solutionOfEigen(solution.data(), vectorSize, 1);

		//use SparseLU to solve sparse matrix problem. SparseLU supports general square sparse systems
		Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> sparseSolver; //or int
		sparseSolver.analyzePattern(this->sparseMatrix);
		sparseSolver.factorize(this->sparseMatrix);
		solutionOfEigen = sparseSolver.solve(rhsOfEigen);

		if (sparseSolver.info() != Eigen::Success)
		{
			SCTM_ASSERT(true, 6);
		}
	}
}

