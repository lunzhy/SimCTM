/**
* @file MatrixSolver.h
* @brief This file defines the class for solving matrix equation with sparse matrix
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
#include <Eigen/OrderingMethods>
#include <vector>

namespace SctmMath
{
	
	/// @brief SctmSparseMatrixSolver is the matrix solver in SimCTM, designed for sparse matrix.
	///
	/// Previously, this solver is inherited in the class that needs matrix solver.
	/// However, this solver serves as a member in the class that uses matrix solver. It is preferable to use this.
	class SctmSparseMatrixSolver
	{
	public:
		
		/// @brief
		enum RefreshMode
		{
			Add,///< add to the previous value
			Cover///< cover the previous value
		};
		Eigen::SparseMatrix<double> matrix; ///< the sparse matrix of the coefficients
		/// @brief SolveMatrix is used to solve the Ax=b problem with determined coefficient sparse matrix
		/// 
		/// The engine of SparseLU is used to solve the sparse matrix. Before solving the matrix equation, the coefficient sparse
		/// matrix has to be set. The matrix equation is solved with input right-hand side b and the solution vector x.
		/// 
		/// @param std::vector<double> & rhs the right-hand side the matrix equation
		/// @param std::vector<double> & solution the solution of the matrix equation
		/// @pre
		/// @return void
		/// @note
		void SolveMatrix(std::vector<double> &rhs, std::vector<double> &solution);
		void RefreshMatrixValue(int _row, int _col, double _value, RefreshMode _mode);
		void RefreshRowOfDirichletBC(int _row);

		void verifySolution(std::vector<double> &rhs, std::vector<double> &solution);
	};
}

#endif