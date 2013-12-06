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
#include <fstream>
#include "SctmMath.h"
#include <map>
#include <vector>
using SctmMath::VectorValue;

//use macro DEBUG to determine if SCTM_ASSERT is defined
#ifdef DEBUG
	#define SCTM_ERROR false // for making SCTM_ASSERT fail
	#define SCTM_ASSERT(cond, err_code) if (!(cond)) { SctmUtils::SctmDebug::ErrorCodeParser(err_code); }
#else
	#define SCTM_ERROR
	#define SCTM_ASSERT(cond, err_code)
#endif // DEBUG


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::fstream;

class FDDomain;

typedef std::map<int, double> VertexMapDouble; // <vertID, phyValue>

namespace SctmUtils
{
	/// @brief This class is used to deal with the timing problems in the simulation
	class SctmTimer
	{
		static const int clockPerSecond = CLOCKS_PER_SEC; ///< the number of inner clock times per second. CLOCKS_PER_SEC is in ctime
	public:
		/// @brief SctmTimer is the construction method of the SctmTimer.
		/// 
		///
		/// 
		/// @pre
		/// @return 
		/// @note
		SctmTimer(): start_time(0), end_time(0), set_time(0) {}
		/// @brief Start is used to start or restart the timer.
		/// 
		///
		/// 
		/// @pre
		/// @return void
		/// @note
		void Start();
		/// @brief Set is used to set the current time in the timer.
		/// 
		/// Set is called to save current time in the beginning of executing certain procedures.
		/// After the procedures are finished. SinceLastSet is used to obtain total time of the procedures.
		/// 
		/// @pre
		/// @return void
		/// @note
		void Set();
		/// @brief End is used to end the timer.
		/// 
		///
		/// 
		/// @pre
		/// @return void
		/// @note
		void End();
		/// @brief SinceLastSet returns the time elapsed from last set of the timer.
		/// 
		/// The timer will still go on when this method is called. It only returns the time elapsed from last set of the timer.
		/// 
		/// @pre
		/// @return double
		/// @note
		double SinceLastSet();
		/// @brief TotalTime is used to return the total time elapsed.
		/// 
		/// This method should be called after End() is called. 
		/// 
		/// @pre
		/// @return double
		/// @note
		double TotalTime();
	protected:
		std::clock_t start_time; ///< the start time of the timer, in inner clock unit
		std::clock_t end_time; ///< the end time of the timer, in inner clock unit
		std::clock_t set_time; ///< the time of setting the timer, in inner clock unit
	};
	
	
	/// @brief SctmTimeStep is used to control the time step in the simulation.
	///
	/// the time used in the simulation is stored in the normalized value.
	class SctmTimeStep
	{
	public:
		SctmTimeStep();
		void GenerateNext();
		double ElapsedTime() const;
		int StepNumber() const;
		double TimeStep() const;
		bool End() const;
	protected:
		double currElapsedTime; /// current time of the simulation
		int currStepNumber; /// current step of the simulation, starting with 1
		double currTimeStep; /// current simulation time step
		
		vector<double> timeSequence;
	protected:
		double getTimeStep();
		double getTimeStep_old();

		void generateTimeSequence();
		bool isEndTime(double time, double endTime);
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
		void PrintDomainDetails(FDDomain *domain);
		void PrintSparseMatrix(const Eigen::SparseMatrix<double> &matrix);
		void PrintSparseMatrixRow(Eigen::SparseMatrix<double> &matrix, int rowIndex);
		void PrintVector(const std::vector<double> &vec, const char *title = "");
		void PrintValue(int i) { std::cout << i << " ";}
		void PrintValue(double d) { std::cout << d << " "; }
		void PrintValue(bool b) { std::cout << (b ? "true" : "false") << " ";}
		void PrintValue(std::string &s) { std::cout << s << " ";}
		void PrintNewLine() { std::cout << std::endl; }
		void PrintBCType(FDBoundary::BCType bcType);
		void PrintDirectionVector(VectorValue &dv);

		void WritePoisson(FDDomain *domain);
		void WriteBandInfo(FDDomain *domain);
		void WriteDensity(FDDomain *domain);
	private:
		bool enable;
	};

	
	/// @brief SctmMessaging deals with the output message and information of the computation process
	class SctmMessaging
	{
	public:
		void PrintWelcomingInformation();
		void PrintHeader(const char *header);
		void PrintTimeElapsed(double time);
		void PrintFileError(const char *filename);
		void PrintDirectoryError();
		void PrintValue(double);
	protected:
		void printLine(string &line);
		void printLine(const char *line);
	};

	
	/// @brief SctmFileOperator provides methods to read and write file.
	///
	/// Currently, the files of input parameters and output results for testing is manipulated using this class.
	/// Writing to one file needs to using the same object. New object will generate new file.
	class SctmFileStream
	{
	public:
		enum FileMode
		{
			Write,
			Read,
			Append,
		};
		SctmFileStream(string _filename, FileMode _mode);
		
		void WriteVector(vector<double> &vec, const char *title = "title not assigned");
		void WriteVector(vector<double> &vec1, vector<double> &vec2, const char *title = "title not assigned");
		void WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, const char *title = "title not assigned");
		void WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> vec4, const char *title = "title not assigned");

		void WriteLine(string &line);
	private:
		string fileName;
	};

	
	/// @brief
	///
	///
	class SctmData
	{
	public:
		SctmData();
		void ReadTunnelParamter();
		void WriteElecCurrDens(vector<FDVertex *> &vertices);
		void WriteElecDens(vector<FDVertex *> &vertices);
		void WritePotential(vector<FDVertex *> &vertices);
		void WriteBandInfo(vector<FDVertex *> &vertices);
		void WriteElecField(vector<FDVertex *> &vertices);
		void WriteTunnelCurrentFromSubs(FDDomain *domain, VertexMapDouble &currDensity);
		void WriteTotalElecDens(vector<FDVertex *> &vertices);
		void WriteTunnelCoeff(FDDomain *domain, VertexMapDouble &inCurrDens, VertexMapDouble &outCurrCoeff);
		void WriteTrapOccupation(vector<FDVertex *> &vertices);
	protected:
		string fileName;
		string directoryName;
		string generateFileSuffix();
	};

	class ConvertToString
	{
	public:
		static string Int(int num);
		static string Double(double num, bool useScientific = true, int numAfterPoionts = 3);
	};

	extern SctmMessaging UtilsMsg;
	extern SctmTimer UtilsTimer;
	extern SctmDebug UtilsDebug;
	extern SctmTimeStep UtilsTimeStep;
	extern SctmData UtilsData;
}

#endif