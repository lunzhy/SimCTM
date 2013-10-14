/**
* @file SctmUtils.h 
* @brief This file contains the utilities used in the simulation
*
*
*
* @author
* @version 
* @date 2013-7-21   16:16
* @note
* @todo
*/
#ifndef _SCTMUTILS_H_
#define _SCTMUTILS_H_

#define DEBUG
#define SCTM_DEBUG_ENABLE true
#include <string>
#include <ctime>
#include <iostream>
#include "DomainDetails.h"
#include <Eigen/Sparse>

//use macro DEBUG to determine if SCTM_ASSERT is defined
#ifdef DEBUG
	#define SCTM_ERROR false
	#define SCTM_ASSERT(cond, err_code) if (!(cond)) { SctmUtils::SctmDebug::ErrorCodeParser(err_code); }
#else
	#define SCTM_ERROR
	#define SCTM_ASSERT(cond, err_code)
#endif // DEBUG


using std::string;
using std::vector;
using std::cout;
using std::endl;
class FDDomain;
class FDBoundary;

namespace SctmUtils
{
	/// @brief This class is used to deal with the timing problems in the simulation
	class SctmTimer
	{
		static const int clockPerSecond = CLOCKS_PER_SEC;
	public:
		SctmTimer(): start_time(0), end_time(0), set_time(0) {}
		void Start();
		void Set();
		void End();
		double SinceLastSet();
		double TotalTime();
	protected:
		std::clock_t start_time;
		std::clock_t end_time;
		std::clock_t set_time;
	};
	
	
	/// @brief SctmTimeStep is used to control the time step in the simulation.
	class SctmTimeStep
	{
	public:
		double NextTimeStep();
	};

	/// @brief The methods used in debugging are defined in this class.
	///
	/// This class is also used to observe the intermediate results during the simulation
	class SctmDebug
	{
	public:
		/// @brief SctmDebug is the construction method for this class
		/// 
		/// The member of enable is initialized in this method.
		/// 
		/// @pre
		/// @return 
		/// @note
		SctmDebug(): enable(SCTM_DEBUG_ENABLE) {}
		/// @brief PrintErrorInfo is used to output message to console with given message.
		/// 
		///
		/// 
		/// @param string msg
		/// @pre
		/// @return void
		/// @note
		static void PrintErrorInfo(string msg);
		/// @brief ErrorCodeParser is used to parse the error code generated by assertion.
		/// 
		///
		/// 
		/// @param int err_code
		/// @pre
		/// @return void
		/// @note
		static void ErrorCodeParser(int err_code);
		void PrintDomainDetails(FDDomain &domain);
		void PrintSparseMatrix(Eigen::SparseMatrix<double> &matrix);
		void PrintSparseMatrixRow(Eigen::SparseMatrix<double> &matrix, int rowIndex);
		void PrintVector(const std::vector<double> &vec, const char *title = "");
		void PrintValue(int i) { std::cout << i << " ";}
		void PrintValue(double d) { std::cout << d << " "; }
		void PrintValue(bool b) { std::cout << (b ? "true" : "false") << " ";}
		void PrintValue(std::string &s) { std::cout << s << " ";}
		void PrintBCType(FDBoundary &bc);
	private:
		bool enable;
	};

	//TODO: a common message class is needed to output the process of the computation
	class SctmMessaging
	{
	public:
		void PrintWelcomingInformation();
		void PrintHeader(const char *header);
		void PrintTimeElapsed(double time);
		void PrintFileError(const char *filename);
	protected:
		void printLine(string &line);
		void printLine(const char *line);
	};

	class SctmFileOperator
	{
	public:
		enum Mode
		{
			Write,
			Read
		};
		SctmFileOperator(string _filename, Mode _mode);
		void WriteVector(vector<double> &vec, const char *title);
		void Write2DVectorForOrigin(vector<double> &vecX, vector<double> &vecY, vector<vector<double>> &vector2D, const char *title);
		void ReadTunnelParameter(vector<double> &cbedges, vector<double> &elecfields);
	private:
		string fileName;
	};

	extern SctmMessaging UtilsMsg;
	extern SctmTimer UtilsTimer;
	extern SctmDebug UtilsDebug;
}

#endif