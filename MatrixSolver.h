/**
* @file MatrixSolver.h
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-14   15:22
* @note
* @todo
*/

#pragma once
#include <Eigen/Sparse>
#include <vector>

namespace SctmMath
{
	class SparseMatrixSolver
	{
	protected:
		Eigen::SparseMatrix<double> sparseMatrix;
		void SolveMatrix(std::vector<double> &rhs, std::vector<double> &solution);
	};
}