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
//#include <Eigen/UmfPackSupport>

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
		Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int>> sparseSolver; //or int
		sparseSolver.analyzePattern(this->matrix);
		sparseSolver.factorize(this->matrix);
		solutionOfEigen = sparseSolver.solve(rhsOfEigen);

		/*
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> sparseSolver(this->matrix);
		sparseSolver.analyzePattern(this->matrix);
		sparseSolver.setTolerance(1e-30);
		sparseSolver.setMaxIterations(1000);
		solutionOfEigen = sparseSolver.solve(rhsOfEigen);
		std::cout << "#iterations:      " << sparseSolver.iterations() << std::endl;
		std::cout << "#estimated error: " << sparseSolver.error() << std::endl;
		*/

		if (sparseSolver.info() != Eigen::Success)
		{
			SCTM_ASSERT(SCTM_ERROR, 10006);
		}
		verifySolution(rhs, solution);
	}

	void SctmSparseMatrixSolver::RefreshMatrixValue(int _row, int _col, double _value, RefreshMode _mode)
	{
		//-------------------------------------------------------------------------------
		//the following method is used to iterate the non-zero coefficient of the matrix
		//for (int k=0; k<sparseMatrix.outerSize(); ++k)
		//	for (SparseMatrix<double>::InnerIterator it(sparseMatrix,k); it; ++it)
		//	{
		//		it.valueRef() = 9; // for get the reference of the coefficient
		//		it.row(); // get the row index
		//		it.col(); // get the column index (here it is equal to k)
		//		it.index(); // inner index, here it is equal to it.row()
		//	}
		//--------------------------------------------------------------------------------
		for (int k = 0; k < this->matrix.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
			{
				if ( (it.row() == _row) && (it.col() == _col))
				{
					if ( _mode == Add )
						it.valueRef() += _value;
					else if ( _mode == Cover )
						it.valueRef() = _value;
				}
			}
	}

	void SctmSparseMatrixSolver::RefreshRowOfDirichletBC(int _row)
	{
		//-------------------------------------------------------------------------------
		//the following method is used to iterate the non-zero coefficient of the matrix
		//for (int k=0; k<sparseMatrix.outerSize(); ++k)
		//	for (SparseMatrix<double>::InnerIterator it(sparseMatrix,k); it; ++it)
		//	{
		//		it.valueRef() = 9; // for get the reference of the coefficient
		//		it.row(); // get the row index
		//		it.col(); // get the column index (here it is equal to k)
		//		it.index(); // inner index, here it is equal to it.row()
		//	}
		//--------------------------------------------------------------------------------
		for (int k = 0; k < this->matrix.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
			{
				if ( it.row() == _row )
				{
					if ( it.col() == _row )
					{
						it.valueRef() = 1;
					}
					else
					{
						it.valueRef() = 0;
					}
				}
			}
	}

	void SctmSparseMatrixSolver::verifySolution(std::vector<double> &rhs, std::vector<double> &solution)
	{
		int vectorSize = rhs.size();
		Eigen::Map<Eigen::VectorXd> rhsOfEigen(rhs.data(), vectorSize, 1);
		Eigen::Map<Eigen::VectorXd> solutionOfEigen(solution.data(), vectorSize, 1);
		Eigen::VectorXd residue;

		using namespace SctmMath;
		residue = this->matrix * solutionOfEigen - rhsOfEigen;
		double modulus = SctmMath::sqrt(residue.dot(residue));

		using namespace SctmUtils;
		string msg = "modulus : " + SctmConverter::DoubleToString(modulus);
		SctmMessaging::Get().PrintMessageLine(msg);
	}
}

