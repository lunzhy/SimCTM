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

extern "C" void pardisoinit(void   *, int    *, int *, int *, double *, int *);
extern "C" void pardiso(void   *, int    *, int *, int *, int *, int *,
	double *, int    *, int *, int *, int *, int *,
	int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec(int *, int *, double *, int *);
extern "C" void pardiso_printstats(int *, int *, double *, int *, int *, int *, double *, int *);

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
		//verifySolution(rhs, solution);
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

	int SctmSparseMatrixSolver::solveWithPardiso(std::vector<double> &rhs, std::vector<double> &solution)
	{
		Eigen::SparseMatrix<double, Eigen::RowMajor> rowMajorMatrix = this->matrix;

		int n = rowMajorMatrix.rows(); //rows==cols
		int *ia = rowMajorMatrix.outerIndexPtr();
		int *ja = rowMajorMatrix.innerIndexPtr();
		double *a = rowMajorMatrix.valuePtr();

		/* RHS and solution vectors. */
		rhs.resize(n);
		solution.resize(n);

		double *b = rhs.data();
		double *x = solution.data();

		/* Set right hand side to one. */
		for (int i = 0; i < n; ++i) {
			b[i] = i;
		}

		int nnz = ia[n];
		int mtype = 11;        /* Real unsymmetric matrix */
		int nrhs = 1;          /* Number of right hand sides. */

		/* Internal solver memory pointer pt,                  */
		/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
		/* or void *pt[64] should be OK on both architectures  */
		void *pt[64];

		/* Pardiso control parameters. */
		int iparm[64];
		double dparm[64];
		int solver;
		int maxfct, mnum, phase, error, msglvl;

		/* Number of processors. */
		int num_procs;

		/* Auxiliary variables. */
		char *var;
		int i;

		double ddum;              /* Double dummy */
		int idum;              /* Integer dummy. */

		/* -------------------------------------------------------------------- */
		/* ..  Setup Pardiso control parameters and initialize the solvers      */
		/*     internal address pointers. This is only necessary for the FIRST   */
		/*     call of the PARDISO solver.                                      */
		/* ---------------------------------------------------------------------*/

		error = 0;
		solver = 0; /* use sparse direct solver */
		pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

		if (error != 0)
		{
			if (error == -10)
				printf("No license file found \n");
			if (error == -11)
				printf("License is expired \n");
			if (error == -12)
				printf("Wrong username or hostname \n");
			return 1;
		}
		else
			printf("[PARDISO]: License check was successful ... \n");


		/* Numbers of processors, value of OMP_NUM_THREADS */
		var = std::getenv("OMP_NUM_THREADS");
		if (var != NULL)
			sscanf(var, "%d", &num_procs);
		else {
			printf("Set environment OMP_NUM_THREADS to 1");
			exit(1);
		}
		iparm[2] = num_procs;


		maxfct = 1;         /* Maximum number of numerical factorizations.  */
		mnum = 1;         /* Which factorization to use. */

		msglvl = 1;         /* Print statistical information  */
		error = 0;         /* Initialize error flag */


		/* -------------------------------------------------------------------- */
		/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
		/*     notation.                                                        */
		/* -------------------------------------------------------------------- */
		for (i = 0; i < n + 1; ++i) {
			ia[i] += 1;
		}
		for (i = 0; i < nnz; ++i) {
			ja[i] += 1;
		}

		//solve the system
		pardiso(pt, &maxfct, &mnum, &mtype, &phase,
			&n, a, ia, ja, &idum, &nrhs,
			iparm, &msglvl, b, x, &error, dparm);

		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			exit(3);
		}

		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i]);
		}
		printf("\n");
	}

	void SctmSparseMatrixSolver::PardisoTest()
	{
		Eigen::SparseMatrix<double> aColMajor;
		aColMajor.resize(5, 5);
		aColMajor.insert(0, 1) = 1;
		aColMajor.insert(0, 4) = 2;
		aColMajor.insert(1, 3) = 3;
		aColMajor.insert(2, 2) = 4;
		aColMajor.insert(3, 1) = 5;
		aColMajor.insert(4, 3) = 6;
		aColMajor.insert(4, 4) = 7;

		this->matrix = aColMajor;

		std::vector<double> rhs, solution;
		solveWithPardiso(rhs, solution);
	}

}

