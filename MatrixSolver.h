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

#ifndef _MATRIXSOLVER_H_
#define _MATRIXSOLVER_H_

#include <Eigen/Sparse>
#include <vector>

namespace SctmMath
{
	class SctmSparseMatrixSolver
	{
	protected:
		Eigen::SparseMatrix<double> sparseMatrix;
		void SolveMatrix(std::vector<double> &rhs, std::vector<double> &solution);
	};
}

#endif