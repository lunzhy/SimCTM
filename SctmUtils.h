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
#include "Material.h"

using SctmMath::VectorValue;

//use macro DEBUG to determine if SCTM_ASSERT is defined
#ifdef DEBUG
	#define SCTM_ERROR false // for making SCTM_ASSERT fail
	#define SCTM_ASSERT(cond, err_code) if (!(cond)) { SctmUtils::SctmDebug::ErrorCodeParser(err_code); }
#else
	#define SCTM_ERROR
	#define SCTM_ASSERT(cond, err_code)
#endif // DEBUG

#ifdef WIN32
	#define SCTM_ENV "Windows"
#else
	#define SCTM_ENV "Linux" 
#endif // WIN32

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::fstream;
using MaterialDB::Mat;

class FDDomain;
class OneDimSubsSolver;
class SolverPack;

typedef std::map<int, double> VertexMapDouble; // <vertID, phyValue>
typedef std::map<string, double> KeywordsTimerMap; // keywords timer

namespace SctmUtils
{
	class SctmData;
	/// @brief This class is used to deal with the timing problems in the simulation
	class SctmTimer
	{
		friend class SctmData;
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
		static SctmTimer& Get();
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
		double PopLastSet();
		/// @brief TotalTime is used to return the total time elapsed.
		/// 
		/// This method should be called after End() is called. 
		/// 
		/// @pre
		/// @return double
		/// @note
		double TotalTime();
		void Timeit(string keywords, double time);
	protected:
		std::clock_t start_time; ///< the start time of the timer, in inner clock unit
		std::clock_t end_time; ///< the end time of the timer, in inner clock unit
		std::clock_t set_time; ///< the time of setting the timer, in inner clock unit

		KeywordsTimerMap keywordTimer;
		vector<double> setList;
	};
	
	
	/// @brief SctmTimeStep is used to control the time step in the simulation.
	///
	/// the time used in the simulation is stored in the real value.
	class SctmTimeStep
	{
	public:
		SctmTimeStep();
		static SctmTimeStep& Get();
		void GenerateNext();
		void Reset();
		double ElapsedTime() const;
		int StepNumber() const;
		bool End() const;
		bool IsMajorTime();
		bool IsGateVoltageChanged();
		bool IsStepWriteData();

		double TimeStep() const;
		double VoltageCellA() const;
		double VoltageCellB() const;
		double VoltageCellC() const;
	protected:
		double temperature;
		double currElapsedTime; /// current time of the simulation
		int currStepNumber; /// current step of the simulation, 0 for initial condition, solver pack starting with 1
		double currTimeStep; /// current simulation time step
		
		vector<double> timeSequence; ///< simulation time sequence, in real value
		vector<double> VgSequenceCellA; ///< gate voltage sequence of cell 1, gate voltage in Single cell, in real value
		vector<double> VgSequenceCellB; ///< gate voltage sequence of cell 2, in real value
		vector<double> VgSequenceCellC; ///< gate voltage sequence of cell 3, in real value
	protected:
		double getTimeStep();
		bool isEndTime(double time, double endTime);
		void readTimestep();
		void generateTimeSequence();
		void fillTimeStepInsideSpan(double start, double end, double vg1, double vg2 = 0, double vg3 = 0);
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
		SctmDebug();
		static SctmDebug& Get();
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
		void PrintValue(std::string s) { std::cout << s << " "; }
		void PrintNewLine() { std::cout << std::endl; }
		void PrintBCType(FDBoundary::BCType bcType);
		void PrintDirectionVector(VectorValue &dv);

		void WritePoisson(FDDomain *domain);
		void WriteBandInfo(FDDomain *domain);
		void WriteDensity(FDDomain *domain);
		void WritePooleFrenkel(FDDomain *domain);
		void WriteMatrixEquation(Eigen::SparseMatrix<double> &matrix, std::vector<double> &rhs, std::vector<double> &solution);
	private:
		bool enable;
		double temperature;
	};

	
	/// @brief SctmMessaging deals with the output message and information of the computation process
	class SctmMessaging
	{
	public:
		static SctmMessaging& Get();
		void PrintWelcomingInformation();
		void PrintHeader(const char *header);
		void PrintTimeElapsed(double time);
		void PrintFileError(const char *filename, const char *msg = "");
		void PrintDirectoryError();
		void PrintInvalidLineInParFile(const char *filename, int lineNum);
		void PrintInvalidParameterName(string& name);
		void PrintValue(double);
		void PrintMessageLine(string line);
	protected:
		void printLine(string line);
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
		void WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> vec4, vector<double> vec5, const char *title = "title not assigned");
		void WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> vec4, vector<double> vec5, vector<double> vec6, const char *title = "title not assigned");
		void WriteVector(vector<int> &vec1, vector<double> &vec2, vector<double> &vec3, const char *title = "title not assigned");
		void WriteVector(vector<int> &vec1, vector<int> &vec2, vector<double> &vec3, const char *title = "title not assigned");
		void WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<string> vec4, const char *title = "title not assigned");
		void WriteLine(string &line);

		void ReadVector(vector<int> &vec1, vector<double> &vec2, vector<double> &vec3);
		void ReadVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> &vec4);
		void ReadVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> &vec4, vector<double> &vec5, vector<double> &vec6);

		static bool FileExisted(string _filename);

	private:
		string fileName;
	};

	
	/// @brief
	///
	///
	class SctmData
	{
	public:
		enum ehInfo
		{
			eInfo,
			hInfo
		};
		SctmData();
		static SctmData& Get();
		void WriteCurrDens(vector<FDVertex *> &vertices, ehInfo ehinfo);
		void WriteCarrierDens(vector<FDVertex *> &vertices, ehInfo ehinfo);
		void WritePotential(vector<FDVertex *> &vertices);
		void WriteBandInfo(vector<FDVertex *> &vertices);
		void WriteElecField(vector<FDVertex *> &vertices);
		void WriteTrappedInfo(vector<FDVertex *> &vertices, ehInfo ehinfo);
		void WriteTrapped(vector<FDVertex *> &vertices);

		void WriteTotalCarrierDens(vector<FDVertex *> &vertices);
		void WriteFlatBandVoltageShift(FDDomain *domain);
		void WriteSubstrateResult(OneDimSubsSolver *subsSolver, SolverPack *solverPack, bool singlefile = false);
		void WriteTrapDensity(vector<FDVertex *> &vertices);
		void WriteTunnelInfo(FDDomain *domain, VertexMapDouble &tnnlOxide, VertexMapDouble &blckOxide, ehInfo ehinfo);
		void WriteTimeConstantTAT(vector<FDVertex *> vertices);
		void WriteCurrDensTAT(vector<FDVertex *> vertices);
		void WriteCurrDensTAT2B(vector<FDVertex *> vertices);

		void WriteSubsVertices(vector<FDVertex *> &vertices);
		void WriteVfbShiftEachInterface(FDDomain *domain);
		void WriteTimerInfo(SctmTimer &timer);
		void WritePooleFrenkelDecrease(vector<FDVertex *> &vertices);
		void WritePooleFrenkelInfo();
		void WriteSubstrateFromInput();
		void WriteTrappedDensRegionwise(FDDomain *domain);
		void WriteTunnelOutDensity(FDDomain *domain, VertexMapDouble &tunToSubs, VertexMapDouble &tbToSubs, VertexMapDouble &tunToGate, VertexMapDouble &tbToGate);
		
		//unused
		void WriteTunnelFromSubs(FDDomain *domain, VertexMapDouble &currDensCoeff, ehInfo ehinfo);
		void WriteTunnelCoeff(FDDomain *domain, VertexMapDouble &inCurrDens, VertexMapDouble &outCurrCoeff);

		void ReadSubsInfoFromFile(VertexMapDouble &fermiAboveMap, VertexMapDouble &channelPotMap);
		void ReadTimestep(vector<double> &timestep, vector<double> &vg1, vector<double> &vg2, vector<double> &vg3);
		void ReadTrappedOcc(vector<double>& eOcc, vector<double>& hOcc);
	protected:
		double temperature;
		string fileName;
		string pathSep;
		string directoryName;
		string generateFileSuffix();
	};


	class SctmConverter
	{
	public:
		static string IntToString(int num);
		static string DoubleToString(double num, bool useScientific = true, int numAfterPoionts = 6);
		static double StringToDouble(const string &strVal);
		static int StringToInt(const string &strVal);
		static bool StringToBool(const string &strVal);
	};


	
	/// @brief SctmGlobalControl stores and manages the global parameters used in the simulation
	///
	/// All the parameters stores in this class are in input value. (not real value)
	class SctmGlobalControl
	{
	public:
		static SctmGlobalControl& Get();
		SctmGlobalControl();
		static void SetGlobalControl(string defaultParPath, string prjpath);
	public:
		//simulation settings
		string ProjectDirectory;
		string UserParFile;
		string DefaulParFile;

		//simulation parameters
		string Structure; ///< the simulation structure
		string Coordinate; ///< the coordinate system in the simulation
		string Solver; ///< the matrix solver

		double Temperature; ///< temperature of the simulation, in [K]
		double SimStartTime; ///< simulation start time
		double SimEndTime; ///< simulation end time
		int SimStepsPerDecade; ///< simulation steps per time decade
		string SimTimeStepMode; ///< simulation time step mode
		string SimTimeStepScale; ///< simulation time step scale
		double SimTimeStepMax; ///< minimum of simulation time step
		int SimStepWriteData; ///< interval steps for writing data

		//density parameters
		double ChannelRadius; ///< the radius of the channel when using cylindrical coordinate
		double SubstrateDoping; ///< substrate doping, positive for N-type, negative for P-type
		double ElecUniTrapDens; ///< the uniform electron trap density, in [cm^-3]
		double HoleUniTrapDens; ///< the uniform hole trap density, in [cm^-3]
		double UniTrapDens; ///< the uniform trap density for both carriers
		string TrapDistribution; ///< the distribution 

		//physics models
		string Carriers; ///< carrier type solved in the simulation
		string TrapCaptureModel;
		bool PhysicsMFN; ///< Modified Fowler-Nordheim tunneling
		bool PhysicsB2T; ///< Band-to-Trap tunneling in
		bool PhysicsT2B; ///< Trap-to-Band tunneling out
		bool PhysicsTAT; ///< TAT and TATB tunneling
		string PhysicsPFModel; ///< Poole-Frenkel model

		//for debug
		double TrapOccupation; ///< trap occupation status
		string TrappedCell;
		bool LateralTunneling; ///< lateral tunneling (non-orthogonal) around the gate
		string CallPytaurus; ///< the mode for calling Pytaurus in Linux
		bool UpdateSubstrate; ///< solve and update substrate potential when calling Pytaurus
		bool ClearCarrier;
		bool RetentionAfterPrgrm; ///< Retention after program
		double RetentionEndTime; ///< Retention end time after program
		bool ReadTrappedDist; ///< Read trapped electron/hole distribution from input
		string SubstrateMethod; ///< Read or solve the substrate potential and fermi energy


		//device structure
		//this part should be in the class derived from FDDomain, i.e. SimpleONO class
		//the length are in [nm]
		double GateVoltage; ///< gate voltage, in [V]
		double GateWorkFunction; ///< the work function of gate material, in [eV]
		double XLength;
		double YLengthTunnel;
		double YLengthTrap;
		double YLengthBlock;
		int XGridNum;
		int YGridNumTunnel;
		int YGridNumTrap;
		int YGridNumBlock;
		Mat::Name TunnelMaterial;
		Mat::Name TrapMaterial;
		Mat::Name BlockMaterial;
	protected:
		static void setGlobalCntrl_Directly();
		static void setGlobalCntrl_FromParFile();
	};


	class ParamBase;
	class SctmParameterParser
	{
	public:
		enum ParName
		{
			structure,
			coordinate,
			solver,
			temperature,
			time_start,
			time_end,
			time_stepPerDecade,
			time_stepMode,
			time_stepScale,
			time_stepMax,
			step_write_data,

			subs_radius,
			subs_type,
			subs_doping,
			trap_eDensity,
			trap_hDensity,
			trap_density,
			trap_distribution,
			
			carriers,
			trap_capture,
			physics_mfn,
			physics_b2t,
			physics_t2b,
			physics_pf,
			physics_tat,

			//for debugging
			debug_trap_occupy,
			debug_trap_cell,
			debug_lateral_tunnel,
			debug_call_pytaurus,
			debug_update_substrate,
			debug_clear_carrier,
			debug_rAfterP,
			debug_rEndTime,
			debug_read_trapped,
			debug_substrate_method,

			//for trap-assisted tunneling
			tat_trap_maxdens,
			tat_trap_pos,
			tat_trap_sig,
			tat_trap_xsection,
			tat_trap_energy,
			tat_t2b_frequency,

			//single cell structure
			sc_gate_voltage,
			sc_gate_workfunction,

			sc_width_value,
			sc_width_grid,

			sc_tunnel_thick,
			sc_tunnel_grid,
			sc_tunnel_material,

			sc_trap_thick,
			sc_trap_grid,
			sc_trap_material,
			
			sc_block_thick,
			sc_block_grid,
			sc_block_material,

			//triple cell structure
			tc_gate1_voltage,
			tc_gate1_workfunction,
			tc_gate1_width,
			tc_gate1_width_grid,

			tc_gate2_voltage,
			tc_gate2_workfunction,
			tc_gate2_width,
			tc_gate2_width_grid,

			tc_gate3_voltage,
			tc_gate3_workfunction,
			tc_gate3_width,
			tc_gate3_width_grid,

			tc_iso_material,
			tc_iso_thick,
			tc_iso_thick_grid,
			tc_iso1_width,
			tc_iso1_width_grid,
			tc_iso2_width,
			tc_iso2_width_grid,
			tc_iso3_width,
			tc_iso3_width_grid,
			tc_iso4_width,
			tc_iso4_width_grid,
			
			tc_tunnel_thick,
			tc_tunnel_thick_grid,
			tc_tunnel_material,

			tc_trap_thick,
			tc_trap_thick_grid,
			tc_trap_material,

			tc_block_thick,
			tc_block_thick_grid,
			tc_block_material,

			//material parameters
			Si_bandgap,
			Si_dielectricConstant,
			Si_electronAffinity,
			Si_eMass,
			Si_eDOSMass,
			Si_eMobility,
			Si_hMass,
			Si_hDOSMass,
			Si_hMobility,

			SiO2_bandgap,
			SiO2_dielectricConstant,
			SiO2_electronAffinity,
			SiO2_eMass,
			SiO2_hMass,

			Al2O3_bandgap,
			Al2O3_dielectricConstant,
			Al2O3_electronAffinity,
			Al2O3_eMass,
			Al2O3_hMass,

			Si3N4_bandgap,
			Si3N4_dielectricConstant,
			Si3N4_highFrequencyDielConst,
			Si3N4_electronAffinity,
			Si3N4_eMass,
			Si3N4_eDOSMass,
			Si3N4_eMobility,
			Si3N4_eTrapEnergy,
			Si3N4_eXsection,
			Si3N4_eTrapXsection,
			Si3N4_eFrequencyT2B,
			Si3N4_eFrequencyPF,
			Si3N4_hMass,
			Si3N4_hDOSMass,
			Si3N4_hMobility,
			Si3N4_hXsection,
			Si3N4_hTrapXsection,
			Si3N4_hTrapEnergy,
			Si3N4_hFrequencyT2B,
			Si3N4_hFrequencyPF,

			HfO2_bandgap,
			HfO2_dielectricConstant,
			HfO2_highFrequencyDielConst,
			HfO2_electronAffinity,
			HfO2_eMass,
			HfO2_eDOSMass,
			HfO2_eMobility,
			HfO2_eTrapEnergy,
			HfO2_eXsection,
			HfO2_eTrapXsection,
			HfO2_eFrequencyT2B,
			HfO2_eFrequencyPF,
			HfO2_hDOSMass,
			HfO2_hMass,
			HfO2_hMobility,
			HfO2_hXsection,
			HfO2_hTrapXsection,
			HfO2_hTrapEnergy,
			HfO2_hFrequencyT2B,
			HfO2_hFrequencyPF
		};

		SctmParameterParser();
		static SctmParameterParser& Get();
		ParamBase *GetPar(ParName name);

	protected:
		void ReadParFile(string &file, std::map<ParName, ParamBase*> &mapToSet);
		std::map<ParName, ParamBase*> defaultParMap;
		std::map<ParName, ParamBase*> userParMap;

		string defaultParFile;
		string userParFile;

		void parseParValue(std::map<ParName, ParamBase*> &mapToSet, string &name, string &valStr);
		static void checkParValue();

		bool isCommentOrSpaceLine(string &line);
		bool isValid(string &line);
		string getParToken(string &line);
		string getParValStr(string &line);
		string trimSpace(string &line);
		string trimComment(string &line);
		static bool isStringValidChoice(string& val, vector<string>& choices);
	};


	class ParamBase
	{
	public:
		ParamBase(SctmParameterParser::ParName name) : parName(name) { }
	protected:
		SctmParameterParser::ParName parName;
		// for the requirement of a polymorphic class in dynamic_cast
		virtual void dummy(){}
	};


	template <typename T>
	class Param : public ParamBase
	{
	public:
		Param(SctmParameterParser::ParName name, const T &val)
			:ParamBase(name), parVal(val) { }
		T Value() { return parVal; }
	private:
		T parVal;
	};


	class BadParConversion : public std::runtime_error
	{
	public:
		BadParConversion(const string &s) : std::runtime_error(s) { }
	};


	
	/// @brief Environment-dependent variables
	///
	///
	class SctmEnv
	{
	public:
		SctmEnv();
		static SctmEnv& Get();
		static bool IsLinux();
		static bool IsWindows();

		string DebugPrjPath;
		string DefaultParamPath;
		string ClearPrjPyPath;

		string PathSep; ///< path separator
		string PytaurusPath;
	};


	class SctmPyCaller
	{
	public:
		static void PyPrepare(string folderPath);
		static void PyClean(string folderPath);
		static void PySolve();
		static void PyBuildStructure();
		static void PyParseAvgVfb();

	protected:
		string pyPath;
	};

}

#endif