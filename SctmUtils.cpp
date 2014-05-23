/**
* @file SctmUtils.cpp
* @brief
*
*
*
* @author
* @version 
* @date 2013-8-9   1:20
* @note
* @todo
*/
#include "SctmUtils.h"
#include "FDDomain.h"
#include "DomainDetails.h"
#include "Normalization.h"
#include <iomanip>
#include <cstring>
#include <fstream>
#include "SctmPhys.h"
#include "SubstrateSolver.h"
#include "Material.h"
#include <sstream>

using std::cout;
using std::endl;
using std::fstream;
using std::stringstream;
using SctmPhys::PhysProperty;
using SctmPhys::TrapProperty;
using MaterialDB::Mat;

namespace SctmUtils
{
	void SctmTimer::Start()
	{
		start_time = clock();
		set_time = 0;
		end_time = 0;
		return;
	}

	void SctmTimer::Set()
	{
		//set_time = clock();
		clock_t current_time = clock();
		setList.push_back(current_time);
		return;
	}

	void SctmTimer::End()
	{
		end_time = clock();
		return;
	}

	double SctmTimer::PopLastSet()
	{
		double time = 0;
		clock_t current_time = clock();
		if (setList.size() == 0)
		{
			SCTM_ASSERT(SCTM_ERROR, 10046);
		}
		double lastSet = setList.back();
		setList.pop_back();
		//TODO: is this conversion correct?
		//time = (double)((current_time - set_time) / (clock_t)clockPerSecond);
		time = (double)(current_time - lastSet) / (clock_t)clockPerSecond;
		return time;
	}

	double SctmTimer::TotalTime()
	{
		double time = 0;
		time = (double)(end_time - start_time ) / (clock_t)clockPerSecond;
		return time;
	}

	SctmTimer& SctmTimer::Get()
	{
		static SctmTimer instance;
		return instance;
	}

	void SctmTimer::Timeit(string keywords, double time)
	{
		keywordTimer[keywords] += time;
	}




	void SctmDebug::PrintErrorInfo(string msg)
	{
		cout << msg << endl;
		exit(0);
	}

	void SctmDebug::ErrorCodeParser(int err_code)
	{
		string msg;
		switch (err_code)
		{
		case 10001:
			msg = "[SctmPhys.cpp] Can not get the specified physical property.";
			break;
		case 10002:
			msg = "[Material.cpp] Non-existed material property.";
			break;
		case 10003:
			msg = "[DomainDetails.cpp] Not rectangular element in constructing elements.";
			break;
		case 10004:
			msg = "[SctmPhys.cpp] Error in calculating vertex-related physical value using material-base property.";
			break;
		case 10005:
			msg = "[MatrixSolver.cpp] The size of right-hand side vector and solution vector are not equal in matrix solver.";
			break;
		case 10006:
			msg = "[MatrixSolver.cpp] The solver of SparseLU fails to solver matrix equation.";
			break;
		case 10007:
			msg = "[PoissonSolver.cpp] Unsuccessful insertion of pair into vertex map occurred.";
			break;
		case 10008:
			msg = "[PoissonSolver.cpp] Error found in vertex map";
			break;
		//case 10009:
		//	msg = "[DomainDetails.cpp] Error occurred in setting boundary condition.";
		//	break;
		case 10010:
			msg	= "[DomainDetails.cpp] Could not find the boundary condition name or the required boundary condition is not set.";
			break;
		case 10011:
			msg = "[DDSolver.cpp] Unsuccessful insertion of pair into vertex or physical parameter map occurred.";
			break;
		case 10012:
			msg = "[DDSolver.cpp] Error found in vertex map when filling coefficient matrix.";
			break;
		case 10013:
			msg = "[DDSolver.cpp] Bad boundary condition (not Dirichlet BC) configuration at the boundary.";
			break;
		case 10014:
			msg = "[DomainDetails.cpp] Refreshing boundary conditions meet non-existed boundary condition.";
			break;
		case 10015:
			msg = "[DDSolver.cpp] Zero value of delta X or Y occurred when filling the DD matrix.";
			break;
		case 10016:
			msg = "[DomainDetails.cpp] Neumann and Cauchy boundary condition encounters normal vector of (0,0).";
			break;
		case 10017:
			msg = "[DomainDetails.cpp] SetBndCond encounters existed boundary conditions.";
			break;
		case 10018:
			msg = "[SctmMath.h] Try to normalize a zero vector or to get the norm X or Y of a zero vector.";
			break;
		case 10019:
			msg = "[SctmPhy.h] Fail to set the physical property value that depends on other value.";
			break;
		case 10020:
			msg = "[TunnelSolver.cpp] Error occurs with the match between input interface vertices and the tunneling beginning vertices.";
			break;
		case 10021:
			msg = "[TunnelSolver.cpp] Emin/Emax error in tunneling solver.";
			break;
		case 10022:
			msg = "[DDSolver.cpp] Errors occurred when passing the boundary current density from tunneling layer to drift-diffusion solver.";
			break;
		case 10023:
			msg = "[DomainDetails.cpp] Invalid boundary direction.";
			break;
		case 10024:
			msg = "[SctmPhys.cpp] Error occurs in calculating electric field or current density. Zero value of length is obtained.";
			break;
		case 10025:
			msg = "[SctmPhys.cpp] Calculating electric field or current density meets non-existent vertex. The corresponding grid number should be larger than 3.";
			break;
		case 10026:
			msg = "[DomainDetails.cpp] Can not set or get tunneling tag to this vertex.";
			break;
		case 10027:
			msg = "[DDSolver] Error occurred when reading tunneling current for boundary condition.";
			break;
		case 10028:
			msg = "[SctmPhys.cpp] Can not set the required trap property.";
			break;
		case 10029:
			msg = "[SctmPhys.cpp] Can not get the required trap property.";
			break;
		case 10030:
			msg = "[FDDomain.cpp] Can not find the required contact.";
			break;
		case 10031:
			msg = "[SubstrateSolver.cpp] The denominator is too small when using Newton Method.";
			break;
		case 10032:
			msg = "[SubstrateSolver.cpp] Solving surface potential meets maximum iteration.";
			break;
		case 10033:
			msg = "[SolverPack.cpp] Error occurs in processing substrate solver result.";
			break;
		case 10034:
			msg = "[Parameter file] Invalid line in parameter file or invalid parameter value.";
			break;
		case 10035:
			msg = "[SctmUtils.cpp] Non-existed parameter is required.";
			break;
		case 10036:
			msg = "[Parameter file] Invalid trap capture model.";
			break;
		case 10037:
			msg = "[Parameter file] Default Parameter file does not exist.";
			break;
		case 10038:
			msg = "[SctmPhys.cpp] This kind of vertex physics does not have multi properties";
			break;
		case 10039:
			msg = "[SctmPhy.cpp] This material name is not included in the multi properties map.";
			break;
		case 10040:
			msg = "[SctmPhy.cpp] Error occurs when setting multi properties for vertices.";
			break;
		case 10041:
			msg = "[TunnelSolver.cpp] This direction has not been considered.";
			break;
		case 10042:
			msg = "[SimpleONO.cpp] Error occurred in setting trap distribution.";
			break;
		case 10043:
			msg = "[Parameter file] Invalid type of trap distribution.";
			break;
		case 10044:
			msg = "[Parameter file] Invalid trap occupation status (occupation not in [0, 1]).";
			break;
		case 10045:
			msg = "[Parameter file] Invalid Poole-Frenkel model name.";
			break;
		case 10046:
			msg = "[SctmUtils.cpp] Time solver meets unpaired set/pop methods.";
			break;
		case 10047:
			msg = "[SctmUtils.cpp] Wrong parameter name in param file.";
			break;
		case 10048:
			msg = "[FDDomain.cpp] Can not find the required region.";
			break;
		case 10049:
			msg = "[Material.cpp] Invalid material name to parse";
			break;
		case 10050:
			msg = "[TunnelSolver.cpp] The last vertex is not at contact in TrapToGateTunnelSolver";
			break;
		case 10051:
			msg = "[TunnelSolver.cpp] The specific vertex does not correspond to substrate vertex";
			break;
		case 10052:
			msg = "[TunnelSolver.cpp] Wrong trapping region name in SlopingTunnelTrapToGate class";
			break;
		case 10053:
			msg = "[Parameter] Invalid parameter value for calling pytaurus.";
			break;
		case 10054:
			msg = "[Parameter] Invalid parameter value for time step mode.";
			break;
		case 10055:
			msg = "[SctmUtils.cpp] The required time step file does not exist.";
			break;
		case 10056:
			msg = "[TunnelSolver.cpp] Invalid tunneling tag in electron tunneling solver.";
			break;
		default:
			msg = "Untracked error";
		}
		PrintErrorInfo(msg);
	}

	void SctmDebug::PrintDomainDetails(FDDomain *domain)
	{
		if (!this->enable)
			return;
		using namespace MaterialDB;
		FDVertex *currVert = NULL;
		Normalization norm = Normalization(this->temperature);
		for (size_t iVert = 0; iVert != domain->vertices.size(); ++iVert)
		{
			currVert = domain->GetVertex(iVert);
			PrintValue(currVert->GetID()); cout << " -- ";
			//if (currVert->Trap == NULL)
			//{
				
			//	cout << "not related to trapping layer";
			//}
			//else
			//{
			//	PrintValue(norm.PullArea(currVert->Trap->GetTrapPrpty(TrapProperty::eCrossSection)));
			//	cout << " -- ";
			//	PrintValue(norm.PullEnergy(currVert->Trap->GetTrapPrpty(TrapProperty::EnergyFromCondBand)));
			//	cout << " -- ";
			//	PrintValue(norm.PullDensity(currVert->Trap->GetTrapPrpty(TrapProperty::eTrapDensity)));
			//}
			//PrintValue(currVert->Phys->GetPhysPrpty(PhysProperty::DielectricConstant));
			//PrintValue(currVert->IsAtContact());
			//PrintValue(currVert->IsAtBoundary(FDBoundary::eCurrentDensity));
			//PrintValue(currVert->BndCond.Valid(FDBoundary::eCurrentDensity));
			PrintValue(currVert->WestVertex == NULL ? "N/A" : SctmConverter::IntToString(currVert->WestVertex->GetID()));
			PrintValue(currVert->NorthVertex == NULL ? "N/A" : SctmConverter::IntToString(currVert->NorthVertex->GetID()));
			PrintValue(currVert->EastVertex == NULL ? "N/A" : SctmConverter::IntToString(currVert->EastVertex->GetID()));
			PrintValue(currVert->SouthVertex == NULL ? "N/A" : SctmConverter::IntToString(currVert->SouthVertex->GetID()));
			cout << " -- ";
			PrintValue(currVert->IsAtBoundary(FDBoundary::Potential));
			if (currVert->IsAtBoundary(FDBoundary::Potential))
			{
				PrintDirectionVector(currVert->BndCond.GetBndDirection(FDBoundary::Potential));
				PrintBCType(currVert->BndCond.GetBCType(FDBoundary::Potential));
				if (currVert->BndCond.GetBCType(FDBoundary::Potential) == FDBoundary::BC_Dirichlet)
				{
					PrintValue(norm.PullPotential(currVert->BndCond.GetBCValue(FDBoundary::Potential)));
				}
				else
				{
					PrintValue(norm.PullPotential(currVert->BndCond.GetBCValue(FDBoundary::Potential)));
					PrintDirectionVector(currVert->BndCond.GetBCNormVector(FDBoundary::Potential));
				}
			}
			cout << " -- -- -- ";
			PrintValue(currVert->IsAtBoundary(FDBoundary::eDensity));
			if (currVert->IsAtBoundary(FDBoundary::eDensity))
			{
				PrintDirectionVector(currVert->BndCond.GetBndDirection(FDBoundary::eDensity));
				PrintBCType(currVert->BndCond.GetBCType(FDBoundary::eDensity));
				if (currVert->BndCond.GetBCType(FDBoundary::eDensity) == FDBoundary::BC_Dirichlet)
				{
					PrintValue(norm.PullCurrDens(currVert->BndCond.GetBCValue(FDBoundary::eDensity)));
				}
				else
				{
					PrintValue(norm.PullCurrDens(currVert->BndCond.GetBCValue(FDBoundary::eDensity)));
					PrintDirectionVector(currVert->BndCond.GetBCNormVector(FDBoundary::eDensity));
				}
			}
			//PrintValue(norm.PullLength(currVert->WestLength));
			//PrintValue(norm.PullLength(currVert->NorthLength));
			//PrintValue(norm.PullLength(currVert->EastLength));
			//PrintValue(norm.PullLength(currVert->SouthLength));
			cout << " -- ";
			PrintValue(currVert->NorthwestElem == NULL ? "N/A" : SctmConverter::IntToString(currVert->NorthwestElem->GetID()));
			PrintValue(currVert->NortheastElem == NULL ? "N/A" : SctmConverter::IntToString(currVert->NortheastElem->GetID()));
			PrintValue(currVert->SoutheastElem == NULL ? "N/A" : SctmConverter::IntToString(currVert->SoutheastElem->GetID()));
			PrintValue(currVert->SouthwestElem == NULL ? "N/A" : SctmConverter::IntToString(currVert->SouthwestElem->GetID()));
			cout << " -- ";
			//PrintValue(currVert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea));
			cout << " -- ";
			//PrintValue(currVert->Phys->GetPhysPrpty(PhysProperty::eMobility));
			//PrintValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->Phys->GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->Phys->GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->SouthVertex==NULL ? -1 : currVert->Phys->GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->NorthVertex==NULL ? -1 : currVert->Phys->GetPhysPrpty(PhysProperty::ElectronAffinity));
			//cout << " -- ";
			//PrintValue(currVert->NorthwestElem==NULL ? -1 : GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			//PrintValue(currVert->NortheastElem==NULL ? -1 : GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_Bandgap));;
			//PrintValue(currVert->SouthwestElem==NULL ? -1 : GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			//PrintValue(currVert->SoutheastElem==NULL ? -1 : GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_Bandgap));
			cout << endl;
		}
		cout << endl;
	}

	void SctmDebug::PrintBCType(FDBoundary::BCType bcType)
	{
		if (!this->enable)
			return;
		string typestring;
		switch (bcType)
		{
		case FDBoundary::BC_Dirichlet:
			typestring = "Dirichlet";
			break;
		case FDBoundary::BC_Neumann:
			typestring = "Neumann";
			break;
		case FDBoundary::BC_Cauchy:
			typestring = "Cauchy";
			break;
		case FDBoundary::BC_Artificial:
			typestring = "Artificial";
			break;
		default:
			break;
		}
		PrintValue(typestring);
	}

	void SctmDebug::PrintSparseMatrix(const Eigen::SparseMatrix<double> &matrix)
	{
		if (!this->enable)
			return;
		cout << matrix << endl;
	}

	void SctmDebug::PrintVector(const std::vector<double> &vec, const char *title)
	{
		if (!this->enable)
			return;
		if (title != "")
		{
			cout << "=================================== The vector of " << title 
				<< " ==================================" << endl;
		}
		cout.setf(std::ios::fixed);
		cout.precision(2);
		for (size_t iVec = 0; iVec != vec.size(); ++iVec)
		{
			if (iVec != 0) { cout << "  "; }
			cout <<  "(" << iVec << ")" << vec.at(iVec);
		}
		cout << endl;
		cout.setf(std::ios::fixed);
	}

	void SctmDebug::PrintSparseMatrixRow(Eigen::SparseMatrix<double> &matrix, int rowIndex)
	{
		if (!this->enable)
			return;
		cout << "======================== Row number = " << rowIndex << " of the sparse matrix"
			<< " ========================" << endl;
		std::vector<double> vec(matrix.cols());
		for (int k = 0; k < matrix.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
			{
				if	(it.row() == rowIndex)
				{
					vec.at(it.col()) = it.valueRef();
				}
			}
		PrintVector(vec);
		cout << "========================================================================================" << endl << endl;
		return;
	}

	void SctmDebug::PrintDirectionVector(VectorValue &dv)
	{
		cout << dv << " ";
	}

	void SctmDebug::WritePoisson(FDDomain *domain)
	{
		SctmData::Get().WritePotential(domain->GetVertices());
	}

	void SctmDebug::WriteBandInfo(FDDomain *domain)
	{
		SctmData::Get().WriteBandInfo(domain->GetVertices());
	}

	void SctmDebug::WriteDensity(FDDomain *domain)
	{
		SctmData::Get().WriteElecDens(domain->GetDDVerts());
	}

	SctmDebug::SctmDebug() : enable(SCTM_DEBUG_ENABLE)
	{
		this->temperature = SctmGlobalControl::Get().Temperature;
	}

	SctmDebug& SctmDebug::Get()
	{
		static SctmDebug instance;
		return instance;
	}

	void SctmDebug::WritePooleFrenkel(FDDomain *domain)
	{
		SctmData::Get().WritePooleFrenkelDecrease(domain->GetDDVerts());
	}

	void SctmDebug::WriteMatrixEquation(Eigen::SparseMatrix<double> &matrix, std::vector<double> &rhs, std::vector<double> &solution)
	{
		static bool isWritten = false;
		if (!isWritten)
		{
			isWritten = true;
		}
		else
		{
			return;
		}

		string fileName = SctmEnv::Get().DebugPrjPath + SctmEnv::Get().PathSep + "Miscellaneous" +
			SctmEnv::Get().PathSep + "matrix.txt";
		std::ofstream file(fileName);
		
		Eigen::MatrixXd densMat(matrix);
		file << densMat << endl << endl << endl;

		for (size_t iVec = 0; iVec != rhs.size(); ++iVec)
		{
			file << rhs.at(iVec) << " ";
		}
		file << endl;

		for (size_t iVec = 0; iVec != solution.size(); ++iVec)
		{
			file << solution.at(iVec) << " ";
		}
		file << endl;
	}








	void SctmMessaging::printLine(string line)
	{
		std::cout << line << std::endl;
	}

	void SctmMessaging::printLine(const char *line)
	{
		cout << line << endl;
	}

	void SctmMessaging::PrintWelcomingInformation()
	{
		printLine("==================================================================================");
		printLine("                                                                                  ");
		printLine("                _______. __  .___  ___.   ______ .___________..___  ___.          ");
		printLine("               /       ||  | |   \\/   |  /      ||           ||   \\/   |          ");
		printLine("              |   (----`|  | |  \\  /  | |  ,----'`---|  |----`|  \\  /  |          ");
		printLine("               \\   \\    |  | |  |\\/|  | |  |         |  |     |  |\\/|  |          ");
		printLine("           .----)   |   |  | |  |  |  | |  `----.    |  |     |  |  |  |          ");
		printLine("           |_______/    |__| |__|  |__|  \\______|    |__|     |__|  |__|          ");
		printLine("                                                                                  ");
		printLine("                                                                                  ");
		printLine("==================================================================================");
		printLine("\n");
		return;
	}

	void SctmMessaging::PrintTimeElapsed(double time)
	{
		cout << "Simulation time step: ";
		cout.setf(std::ios::fixed);
		cout.precision(8);
		//string timeStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		cout << SctmTimeStep::Get().ElapsedTime() << "s" << "\t\t";

		string msg = "Time elapsed: ";
		cout.setf(std::ios::fixed);
		cout.precision(3);
		cout << msg << time << "s" << endl << endl;
		cout.unsetf(std::ios::fixed);
		cout.precision(6);
		return;
	}

	void SctmMessaging::PrintHeader(const char *header)
	{
		cout << "===> " << header << endl << endl;
		return;
	}

	void SctmMessaging::PrintFileError(const char *filename, const char *msg)
	{
		cout << endl << "XXXXX=>" << "cannot open or create file" << ' ' << filename << endl;
		cout << msg << endl;
		exit(1);
	}

	void SctmMessaging::PrintDirectoryError()
	{
		cout << endl << "XXXXX=>" << "no directory available" << endl;
		exit(1);
	}

	void SctmMessaging::PrintValue(double num)
	{
		cout << num << endl;
	}

	void SctmMessaging::PrintInvalidLineInParFile(const char *filename, int lineNum)
	{
		cout << endl << "XXXXX=>" << "Invalid line or value in parameter file: " << filename << " line:" << SctmConverter::IntToString(lineNum) << endl;
		exit(1);
	}

	SctmMessaging& SctmMessaging::Get()
	{
		static SctmMessaging instance;
		return instance;
	}

	void SctmMessaging::PrintInvalidParameterName(string& name)
	{
		cout << endl << "XXXXX=>" << "Invalid parameter name: " << name << endl;
		exit(1);
	}

	void SctmMessaging::PrintMessageLine(string line)
	{
		printLine(line);
	}






	SctmFileStream::SctmFileStream(string _filename, FileMode _mode)
	{
		this->fileName = _filename;
		fstream file;
		if (_mode == Write)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file)
			{
				//if the file doesn't exist, create it.
				file.open(this->fileName.c_str(), std::ios::out);
				if (!file)
					SctmMessaging::Get().PrintDirectoryError();
				else
					file.close();
			}
			else
			{
				file.close();//this is important because the existed file has been opened already.
				//if the file exist, truncate it
				file.open(this->fileName.c_str(), std::ios::out | std::ios::trunc);
				file.close();
			}
			return;
		}
		if (_mode == Append)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file)
			{
				//if the file doesn't exist, create it.
				file.open(this->fileName.c_str(), std::ios::out);
				if (!file)
					SctmMessaging::Get().PrintDirectoryError();
				else
					file.close();
			}
			else
			{
				file.close();//this is important because the existed file has been opened already.
			}
			return;
		}
		if (_mode == Read)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file.is_open())
				SctmMessaging::Get().PrintFileError(_filename.c_str());
			else
				file.close();
			return;
		}
	}

	void SctmFileStream::WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<double> &vec, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);

		tofile << title << endl;
		for (size_t ix = 0; ix != vec.size(); ++ix)
		{
			tofile << vec.at(ix);
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> vec4, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << '\t' << vec4.at(iv) << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> vec4, vector<double> vec5, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << '\t' << vec4.at(iv) << '\t' << vec5.at(iv) << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<double> &vec1, vector<double> &vec2, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<int> &vec1, vector<double> &vec2, vector<double> &vec3, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<int> &vec1, vector<int> &vec2, vector<double> &vec3, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << endl;
		}
		tofile.close();
	}

	void SctmFileStream::WriteVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<string> vec4, const char *title /*= "title not assigned"*/)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << title << endl;
		for (size_t iv = 0; iv != vec1.size(); ++iv)
		{
			tofile << vec1.at(iv) << '\t' << vec2.at(iv) << '\t' << vec3.at(iv) << '\t' << vec4.at(iv) << endl;
		}
		tofile.close();
	}




	void SctmFileStream::WriteLine(string &line)
	{
		std::ofstream tofile(this->fileName.c_str(), std::ios::app);
		tofile << line << '\n';
		tofile.close();
	}

	bool SctmFileStream::FileExisted(string _filename)
	{
		std::ifstream file(_filename.c_str());
		if (!file)
		{
			return false;
		}
		else
		{
			file.close();
			return true;
		}
	}

	void SctmFileStream::ReadVector(vector<int> &vec1, vector<double> &vec2, vector<double> &vec3)
	{
		std::ifstream infile(this->fileName.c_str(), std::ios::in);
		int vertID = 0;
		double val1 = 0;
		double val2 = 0;
		string title = "";
		std::getline(infile, title);

		vec1.clear(); vec2.clear(); vec3.clear();
		while (infile >> vertID >> val1 >> val2)
		{
			vec1.push_back(vertID);
			vec2.push_back(val1);
			vec3.push_back(val2);
		}
		infile.close();
		return;
	}

	void SctmFileStream::ReadVector(vector<double> &vec1, vector<double> &vec2, vector<double> &vec3, vector<double> &vec4)
	{
		std::ifstream infile(this->fileName.c_str(), std::ios::in);
		int stepNum = 0;
		double timestep = 0;
		double vg1 = 0;
		double vg2 = 0;
		double vg3 = 0;
		string title = "";
		string line = "";
		std::getline(infile, title);
		std::istringstream is;

		while (getline(infile, line))
		{	
			is.clear();
			is.str(line);
			is >> stepNum >> timestep >> vg1 >> vg2 >> vg3;
			vec1.push_back(timestep);
			vec2.push_back(vg1);
			vec3.push_back(vg2);
			vec4.push_back(vg3);
		}
		infile.close();
		return;
	}







	SctmTimeStep::SctmTimeStep()
	{
		this->temperature = SctmGlobalControl::Get().Temperature;
		this->currStepNumber = 0;
		this->currElapsedTime = 0;
		if (SctmGlobalControl::Get().SimTimeStepMode == "Internal")
		{
			generateTimeSequence();
		}
		else if (SctmGlobalControl::Get().SimTimeStepMode == "UserDefined")
		{
			readTimestep();
		}
		else
		{
			SCTM_ASSERT(SCTM_ERROR, 10054);
		}
	}

	void SctmTimeStep::GenerateNext()
	{
		currStepNumber += 1; // step number starts with 1
		//this->currTimeStep = getTimeStep_old();
		this->currTimeStep = getTimeStep();
		this->currElapsedTime = timeSequence.at(currStepNumber);
	}

	double SctmTimeStep::TimeStep() const
	{
		Normalization norm = Normalization(this->temperature);
		return norm.PushTime(currTimeStep);
	}

	double SctmTimeStep::ElapsedTime() const
	{
		//for output the time elapsed, normalization is not needed.
		//it should be normalized, but for history reasons, it is not normalized.
		return currElapsedTime;
	}

	int SctmTimeStep::StepNumber() const
	{
		return currStepNumber;
	}

	void SctmTimeStep::generateTimeSequence()
	{
		////////// time sequence parameter //////////
		//double startTimeDefault = 1e-12; // previously, 1e-12
		//if (SctmGlobalControl::Get().SimStartTime < startTimeDefault)
		//{
		//	startTimeDefault = SctmGlobalControl::Get().SimStartTime;
		//}
		double startTime = SctmGlobalControl::Get().SimStartTime;
		double endTime = SctmGlobalControl::Get().SimEndTime;
		int stepPerDecade = SctmGlobalControl::Get().SimStepsPerDecade;
		string stepScale = SctmGlobalControl::Get().SimTimeStepScale;
		double maxTimestep = SctmGlobalControl::Get().SimTimeStepMax;
		
		/////////////////////////////////////////////
		ParamBase *parBase = NULL;
		Normalization norm  = Normalization(this->temperature);
		
		//the major time steps
		double vgCellA = 0;
		double vgCellB = 0;
		double vgCellC = 0;
		string simStructure = SctmGlobalControl::Get().Structure;
		if (simStructure == "Triple" || simStructure == "TripleFull")
		{
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_voltage);
			vgCellA = dynamic_cast<Param<double> *>(parBase)->Value();
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_voltage);
			vgCellB = dynamic_cast<Param<double> *>(parBase)->Value();
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_voltage);
			vgCellC = dynamic_cast<Param<double> *>(parBase)->Value();
		}
		else if (SctmGlobalControl::Get().Structure == "Single")
		{
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_gate_voltage);
			vgCellA = dynamic_cast<Param<double> *>(parBase)->Value();
		}
		
		
		//deal with the start time
		//the time span start from 0 to time start only has one step
		timeSequence.push_back(0);
		VgSequenceCellA.push_back(vgCellA);
		VgSequenceCellB.push_back(vgCellB);
		VgSequenceCellC.push_back(vgCellC);
		
		timeSequence.push_back(startTime);
		VgSequenceCellA.push_back(vgCellA);
		VgSequenceCellB.push_back(vgCellB);
		VgSequenceCellC.push_back(vgCellC);

		//deal with the time span
		double startTime_exp = SctmMath::floor(SctmMath::log10(startTime));
		double endTime_exp = SctmMath::floor(SctmMath::log10(endTime));

		bool iterOver = false;
		double currTime_exp = 0;
		double spanStart = startTime;
		double spanEnd = 0;
		double spanStart_exp = startTime_exp;
		double spanEnd_exp = 0;
		while (!iterOver)
		{
			spanEnd_exp = spanStart_exp + 1;
			if (SctmMath::exp10(spanEnd_exp) > endTime || SctmMath::logically_equal(SctmMath::exp10(spanEnd_exp), endTime))
			{
				spanEnd = endTime;
				iterOver = true;
			}
			else
			{
				spanEnd = SctmMath::exp10(spanEnd_exp);
			}
			fillTimeStepInsideSpan(spanStart, spanEnd, vgCellA, vgCellB, vgCellC);
			spanStart = spanEnd;
			spanStart_exp = spanEnd_exp;
		}
	}

	double SctmTimeStep::getTimeStep()
	{
		double next = 0;
		next = timeSequence.at(currStepNumber) - timeSequence.at(currStepNumber - 1);
		return next;
	}

	bool SctmTimeStep::isEndTime(double time, double endTime)
	{
		double roundingError = 0.001;
		bool ret = false;
		ret = (SctmMath::abs((time - endTime) / endTime) < roundingError) || (time > endTime);
		return ret;
	}

	bool SctmTimeStep::End() const
	{
		static int totalEffectiveStep = timeSequence.size() - 1;
		return (currStepNumber == totalEffectiveStep);
	}

	SctmTimeStep& SctmTimeStep::Get()
	{
		static SctmTimeStep instance;
		return instance;
	}

	bool SctmTimeStep::IsMajorTime()
	{
		double eps = 1e-10;
		double diff = 0;
		bool ret = false;
		//currElapsedTime is in real value
		double currentTime = currElapsedTime;
		if (currStepNumber == 1)
		{
			return true;
		}
		diff = SctmMath::abs(SctmMath::log10(currentTime)) - SctmMath::abs(SctmMath::floor(SctmMath::log10(currentTime)));
		return SctmMath::abs(diff) < eps;
	}

	void SctmTimeStep::fillTimeStepInsideSpan(double starttime, double endtime, double vg1, double vg2, double vg3)
	{
		double timeStep = 0;
		double timeStep_exp = 0;
		double stepNum = SctmGlobalControl::Get().SimStepsPerDecade;
		double maxStep = SctmGlobalControl::Get().SimTimeStepMax;
		double currTime = starttime;
		bool iterOver = false;

		if (SctmGlobalControl::Get().SimTimeStepScale == "Linear")
		{
			timeStep = (endtime - starttime) / stepNum;
			if (!SctmMath::logically_equal(maxStep, 0) && timeStep > maxStep)
			{
				timeStep = maxStep;
			}

			currTime = starttime;
			iterOver = false;
			while (!iterOver)
			{
				currTime += timeStep;
				if (currTime > endtime || SctmMath::logically_equal(currTime, endtime))
				{
					currTime = endtime;
					iterOver = true;
				}
				timeSequence.push_back(currTime);
				VgSequenceCellA.push_back(vg1);
				VgSequenceCellB.push_back(vg2);
				VgSequenceCellC.push_back(vg3);
			}
		}
		else if (SctmGlobalControl::Get().SimTimeStepScale == "Exp10")
		{
			
			iterOver = false;
			while (!iterOver)
			{
				timeStep = currTime * //current time can't be 0
					(SctmMath::exp10((SctmMath::log10(endtime) - SctmMath::log10(starttime)) / stepNum) - 1);
				if (!SctmMath::logically_equal(maxStep, 0) && timeStep > maxStep)
				{
					timeStep = maxStep;
				}

				currTime += timeStep;
				if (currTime > endtime || isEndTime(currTime, endtime))
				{
					currTime = endtime;
					iterOver = true;
				}
				timeSequence.push_back(currTime);
				VgSequenceCellA.push_back(vg1);
				VgSequenceCellB.push_back(vg2);
				VgSequenceCellC.push_back(vg3);
			}
		}
	}

	double SctmTimeStep::VoltageCellA() const
	{
		Normalization norm = Normalization(this->temperature);
		double voltage = 0;
		voltage = VgSequenceCellA.at(this->currStepNumber);
		voltage = norm.PushPotential(voltage);
		return voltage;
	}

	double SctmTimeStep::VoltageCellB() const
	{
		Normalization norm = Normalization(this->temperature);
		double voltage = 0;
		voltage = VgSequenceCellB.at(this->currStepNumber);
		voltage = norm.PushPotential(voltage);
		return voltage;
	}

	double SctmTimeStep::VoltageCellC() const
	{
		Normalization norm = Normalization(this->temperature);
		double voltage = 0;
		voltage = VgSequenceCellC.at(this->currStepNumber);
		voltage = norm.PushPotential(voltage);
		return voltage;
	}

	bool SctmTimeStep::IsGateVoltageChanged()
	{
		bool isChanged_gate1 = true;
		bool isChanged_gate2 = true;
		bool isChanged_gate3 = true;

		if (currStepNumber == 0)
		{
			return true;
		}

		double lastGate1 = VgSequenceCellA.at(this->currStepNumber - 1);
		double lastGate2 = VgSequenceCellB.at(this->currStepNumber - 1);
		double lastGate3 = VgSequenceCellC.at(this->currStepNumber - 1);
		double currGate1 = VgSequenceCellA.at(this->currStepNumber);
		double currGate2 = VgSequenceCellB.at(this->currStepNumber);
		double currGate3 = VgSequenceCellC.at(this->currStepNumber);

		isChanged_gate1 = !SctmMath::logically_equal(lastGate1, currGate1);
		isChanged_gate2 = !SctmMath::logically_equal(lastGate2, currGate2);
		isChanged_gate3 = !SctmMath::logically_equal(lastGate3, currGate3);

		return isChanged_gate1 || isChanged_gate2 || isChanged_gate3;

	}

	void SctmTimeStep::readTimestep()
	{
		ParamBase *parBase = NULL;
		//the major time steps
		double vgCellA = 0;
		double vgCellB = 0;
		double vgCellC = 0;
		string simStructure = SctmGlobalControl::Get().Structure;
		if (simStructure == "Triple" || simStructure == "TripleFull")
		{
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate1_voltage);
			vgCellA = dynamic_cast<Param<double> *>(parBase)->Value();
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate2_voltage);
			vgCellB = dynamic_cast<Param<double> *>(parBase)->Value();
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::tc_gate3_voltage);
			vgCellC = dynamic_cast<Param<double> *>(parBase)->Value();
		}
		else if (SctmGlobalControl::Get().Structure == "Single")
		{
			parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_gate_voltage);
			vgCellA = dynamic_cast<Param<double> *>(parBase)->Value();
		}
		//set the initial time step information from parameter file
		timeSequence.push_back(0);
		VgSequenceCellA.push_back(vgCellA);
		VgSequenceCellB.push_back(vgCellB);
		VgSequenceCellC.push_back(vgCellC);

		//read and set the following time step from timestep file
		SctmData::Get().ReadTimestep(timeSequence, VgSequenceCellA, VgSequenceCellB, VgSequenceCellC);
	}









	SctmData::SctmData()
	{
		this->temperature = SctmGlobalControl::Get().Temperature;
		//directoryName = "E:\\PhD Study\\SimCTM\\SctmTest\\SolverPackTest\\";
		this->directoryName = SctmGlobalControl::Get().ProjectDirectory;
		this->pathSep = SctmEnv::Get().PathSep;
	}


	void SctmData::WriteElecDens(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Density" + pathSep + "eDensity" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization(this->temperature);
		vector<double> vecX;
		vector<double> vecY;
		vector<double> vecDen;

		FDVertex *currVert = NULL;
		for (size_t iVer = 0; iVer != vertices.size(); ++iVer)
		{
			currVert = vertices.at(iVer);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			vecDen.push_back(norm.PullDensity(currVert->Phys->GetPhysPrpty(PhysProperty::eDensity)));
		}

		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "electron density of time [" + numStr + "]"; 
		file.WriteVector(vecX, vecY, vecDen, title.c_str());
	}

	void SctmData::WritePotential(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Potential" + pathSep + "potential" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization(this->temperature);
		vector<double> vecX;
		vector<double> vecY;
		vector<double> vecPot;

		FDVertex *currVert = NULL;
		for (size_t iVer = 0; iVer != vertices.size(); ++iVer)
		{
			currVert = vertices.at(iVer);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			vecPot.push_back(norm.PullPotential(currVert->Phys->GetPhysPrpty(PhysProperty::ElectrostaticPotential)));
		}

		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "potential of time [" + numStr + "]"; 
		file.WriteVector(vecX, vecY, vecPot, title.c_str());
	}

	string SctmData::generateFileSuffix()
	{
		string ret = "";
		string step = SctmConverter::IntToString(SctmTimeStep::Get().StepNumber());

		ret = "_s" + step + ".txt";
		return ret;
	}

	void SctmData::WriteBandInfo(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Band" + pathSep + "band" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization(this->temperature);
		vector<double> vecX;
		vector<double> vecY;
		vector<double> vecCB;
		vector<double> vecVB;

		FDVertex *currVert = NULL;
		using MaterialDB::Mat;
		Mat::Name currMatName = Mat::ErrorMaterial;

		for (size_t iVer = 0; iVer != vertices.size(); ++iVer)
		{
			currVert = vertices.at(iVer);

			if (currVert->Phys->HasMultiPrpty(PhysProperty::ConductionBandEnergy))
			{
				vector<Mat::Name> mats = currVert->Phys->GetRelatedMaterialNames();
				for (size_t in = 0; in != mats.size(); ++in)
				{
					currMatName = mats.at(in);
					vecX.push_back(norm.PullLength(currVert->X));
					vecY.push_back(norm.PullLength(currVert->Y));
					vecCB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy, currMatName)));
					vecVB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy, currMatName)));
				}
			}
			else
			{
				vecX.push_back(norm.PullLength(currVert->X));
				vecY.push_back(norm.PullLength(currVert->Y));
				vecCB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy)));
				vecVB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy)));
			}
		}

		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "band structure of time [" + numStr + "] (x, y, Conduction band, Valence band)"; 
		file.WriteVector(vecX, vecY, vecCB, vecVB, title.c_str());
	}

	void SctmData::WriteTunnelCurrentFromSubs(FDDomain *domain, VertexMapDouble &currDensity)
	{
		fileName = directoryName + pathSep + "Current" + pathSep + "subsCurrent" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);
		
		int vertID = 0;
		vector<double> vecX;
		vector<double> vecY;
		vector<double> currDens;
		vector<string> direction;
		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		for (VertexMapDouble::iterator it = currDensity.begin(); it != currDensity.end(); ++it)
		{
			vertID = it->first;
			currVert = domain->GetVertex(vertID);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			currDens.push_back(norm.PullCurrDens(it->second));
			if (currVert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn)
			{
				direction.push_back("in");
			}
			else
			{
				direction.push_back("out");
			}
		}
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "substrate tunneling current density of time [" + numStr + "] (x, y, tunneling current density(in)/tunneling coefficient(out))";
		file.WriteVector(vecX, vecY, currDens, direction, title.c_str());
	}

	void SctmData::WriteElecCurrDens(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Current" + pathSep + "eCurrDens" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<double> vecX;
		vector<double> vecY;
		vector<double> eCurrDens;
		vector<double> eCurrDens_X;
		vector<double> eCurrDens_Y;
		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			eCurrDens.push_back(norm.PullCurrDens(currVert->Phys->GetPhysPrpty(PhysProperty::eCurrentDensity)));
			eCurrDens_X.push_back(norm.PullCurrDens(currVert->Phys->GetPhysPrpty(PhysProperty::eCurrentDensity_X)));
			eCurrDens_Y.push_back(norm.PullCurrDens(currVert->Phys->GetPhysPrpty(PhysProperty::eCurrentDensity_Y)));
		}
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "electron current density of time [" + numStr + "] (x, y, electron current density, eCurret in X, eCurrent in Y)"; 
		file.WriteVector(vecX, vecY, eCurrDens, eCurrDens_X, eCurrDens_Y, title.c_str());
	}

	void SctmData::WriteElecField(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "ElecField" + pathSep + "elecField" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<double> vecX;
		vector<double> vecY;
		vector<double> eElecFiled;
		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			eElecFiled.push_back(norm.PullElecField(currVert->Phys->GetPhysPrpty(PhysProperty::ElectricField)));
		}
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "electric field of time [" + numStr + "] (x, y, electric field)";
		file.WriteVector(vecX, vecY, eElecFiled, title.c_str());
	}

	void SctmData::WriteTotalElecDens(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "totalDensity.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		double area = 0;
		double freeElecDens = 0; // line density
		double trapElecDens = 0;
		Normalization norm = Normalization(this->temperature);
		FDVertex *vert = NULL;
		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			vert = vertices.at(iVert);
			area = vert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea);
			freeElecDens += area * vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
			trapElecDens += area * vert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
		}
		string timeStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string valStrFree = SctmConverter::DoubleToString(norm.PullLineDensity(freeElecDens));
		string valStrTrap = SctmConverter::DoubleToString(norm.PullLineDensity(trapElecDens));
		string valStrTotal = SctmConverter::DoubleToString(norm.PullLineDensity(freeElecDens + trapElecDens));
		string line = timeStr + "\t\t" + valStrFree + "\t\t" + valStrTrap + "\t\t" + valStrTotal;
		file.WriteLine(line);
	}

	void SctmData::WriteTunnelCoeff(FDDomain *domain, VertexMapDouble &inCurrDens, VertexMapDouble &outCurrCoeff)
	{
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "tunCoeff.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		int vertID = 0;
		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		VertexMapDouble::iterator it = inCurrDens.begin(); 
		string currDens = SctmConverter::DoubleToString(norm.PullCurrDens(it->second));
		it = outCurrCoeff.begin();
		string tunCoeff = SctmConverter::DoubleToString(norm.PullTunCoeff(it->second));
		
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string line = numStr + "\t\t" + currDens + "\t\t" + tunCoeff;
		file.WriteLine(line);
	}

	void SctmData::WriteTrappedInfo(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Trap" + pathSep + "trapOccupation" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<double> vecX;
		vector<double> vecY;
		vector<double> eTrappedDens;
		vector<double> occupation;
		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			eTrappedDens.push_back(norm.PullDensity(currVert->Trap->GetTrapPrpty(TrapProperty::eTrapped)));
			occupation.push_back(currVert->Trap->GetTrapPrpty(TrapProperty::eOccupation));
		}
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "occupation of electron trap [" + numStr + "] (x, y, trapped density, trap occupation)";
		file.WriteVector(vecX, vecY, eTrappedDens, occupation, title.c_str());
	}

	void SctmData::WriteFlatBandVoltageShift(FDDomain *domain)
	{
		double VfbShift = 0;
		Normalization norm = Normalization(this->temperature);

		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "VfbShift.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		VfbShift = SctmPhys::CalculateFlatbandShift_domain(domain);

		string timeStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string valStr = SctmConverter::DoubleToString(norm.PullPotential(VfbShift));
		string line = timeStr + "\t\t" + valStr;
		file.WriteLine(line);
	}

	void SctmData::WriteSubstrateResult(OneDimSubsSolver *subsSolver)
	{
		static bool firstRun = true;
		Normalization norm = Normalization(this->temperature);

		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "substrate.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		VertexMapDouble fermi_above_map;
		VertexMapDouble channel_potential_map;
		double fermi_above = 0;
		double channel_potential = 0;

		subsSolver->ReturnResult(fermi_above_map, channel_potential_map);

		if (firstRun)
		{
			string headline = "Time [s]\t\tChannel potential\t\tFermi energy above conduction band";
			file.WriteLine(headline);
			firstRun = false;
		}

		string time = "time = [" + SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime()) + "s]";
		file.WriteLine(time);
		int vertID = 0;
		string line = "";
		for (VertexMapDouble::iterator it = fermi_above_map.begin(); it != fermi_above_map.end(); ++it)
		{
			vertID = it->first;
			fermi_above = it->second;
			channel_potential = channel_potential_map[vertID];
			line = SctmConverter::DoubleToString(norm.PullPotential(channel_potential)) +
				"\t\t" + SctmConverter::DoubleToString(norm.PullPotential(fermi_above));
			file.WriteLine(line);
		}
	}

	SctmData& SctmData::Get()
	{
		static SctmData instance;
		return instance;
	}

	void SctmData::WriteTrapDensity(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "TrapInfo.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<double> vecX;
		vector<double> vecY;
		vector<double> trapDens;

		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;
		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			trapDens.push_back(norm.PullDensity(currVert->Trap->GetTrapPrpty(TrapProperty::eTrapDensity)));
		}

		string title = "trap distribution information";
		file.WriteVector(vecX, vecY, trapDens, title.c_str());
	}

	void SctmData::WriteTimerInfo(SctmTimer &timer)
	{
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "timing.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		double totalTime = timer.keywordTimer["Total"];
		double poissonTime = timer.keywordTimer["Poisson"];
		double ddTime = timer.keywordTimer["Transport"];

		string line = "";
		line = "Total time: " + SctmConverter::DoubleToString(totalTime) + "s";
		file.WriteLine(line);
		line = "Poisson time: " + SctmConverter::DoubleToString(poissonTime) + "s";
		file.WriteLine(line);
		line = "Drift-Diffusion time: " + SctmConverter::DoubleToString(ddTime) + "s";
		file.WriteLine(line);

	}

	void SctmData::WritePooleFrenkelInfo()
	{
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "simInfo.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		double time = SctmTimeStep::Get().ElapsedTime(); //elapsed time has been converted

		string line = "Energy decrease of Poole-Frenkel effect surpasses trap energy at " +
			SctmConverter::DoubleToString(time) + "s.";
		file.WriteLine(line);
	}

	void SctmData::WritePooleFrenkelDecrease(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "Trap" + pathSep + "PooleFrenkelDecrease" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<double> vecX;
		vector<double> vecY;
		vector<double> PF_decrease;

		Normalization norm = Normalization(this->temperature);
		FDVertex *currVert = NULL;

		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			PF_decrease.push_back(norm.PullEnergy(currVert->Trap->GetTrapPrpty(TrapProperty::eTrapEnergyDecreasePF)));
		}
		string numStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "energy decrease of Poole-Frenkel effect at [" + numStr + "] (x, y, energy decrease)";
		file.WriteVector(vecX, vecY, PF_decrease, title.c_str());
	}

	void SctmData::ReadSubsInfoFromFile(VertexMapDouble &fermiAboveMap, VertexMapDouble &channelPotMap)
	{
		fileName = directoryName + pathSep + "exchange" + pathSep + "substrate.in";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Read);
		vector<int> vertID;
		vector<double> pot;
		vector<double> fermiAbove;
		int id = 0;

		file.ReadVector(vertID, pot, fermiAbove);
		Normalization norm = Normalization(this->temperature);
		for (size_t iVert = 0; iVert != vertID.size(); ++iVert)
		{
			id = vertID.at(iVert);
			channelPotMap[id] = norm.PushPotential(pot.at(iVert));
			fermiAboveMap[id] = norm.PushPotential(fermiAbove.at(iVert));
		}
	}

	void SctmData::WriteVertexInfo(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + pathSep + "exchange" + pathSep + "subs_points.in";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);
		vector<int> vertID;
		vector<double> xCoords;
		vector<double> yCoords;

		Normalization norm = Normalization(this->temperature);
		FDVertex *vert = NULL;

		double cm_in_nm = SctmPhys::cm_in_nm;

		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			vert = vertices.at(iVert);
			vertID.push_back(vert->GetID());
			xCoords.push_back(norm.PullLength(vert->X) * cm_in_nm);
			yCoords.push_back(norm.PullLength(vert->Y) * cm_in_nm);
		}

		string title = "vertex ID\t\tx\t\ty\t(coordinates are in [nm])";
		file.WriteVector(vertID, xCoords, yCoords, title.c_str());
	}

	void SctmData::WriteVfbShiftEachInterface(FDDomain *domain)
	{
		fileName = directoryName + pathSep + "Substrate" + pathSep + "VfbInterface" + generateFileSuffix();
		SctmFileStream file_subs = SctmFileStream(fileName, SctmFileStream::Write);
		fileName = directoryName + pathSep + "exchange" + pathSep + "charge.in";
		SctmFileStream file_change = SctmFileStream(fileName, SctmFileStream::Write);



		Normalization norm = Normalization(this->temperature);
		vector<int> leftVertID;
		vector<int> rightVertID;
		vector<double> flatbandVoltageShift;

		static FDContact *subsContact = domain->GetContact("Channel");
		static vector<FDVertex *> &channelVerts = subsContact->GetContactVerts();

		FDVertex *lVert = NULL; //left vertex
		FDVertex *rVert = NULL; //right vertex
		double sum = 0;
		double vfbShift = 0;
		for (size_t iVert = 0; iVert < channelVerts.size()-1; iVert++)
		{
			//iterate the edge
			lVert = channelVerts.at(iVert);
			rVert = channelVerts.at(iVert + 1);

			sum = 0;
			vfbShift = SctmPhys::CalculateFlatbandShift_slice_for2D(lVert);
			if (lVert->WestVertex != NULL)
			{
				sum += vfbShift / 2;
			}
			else
			{
				sum += vfbShift;
			}

			vfbShift = SctmPhys::CalculateFlatbandShift_slice_for2D(rVert);
			if (rVert->EastVertex != NULL)
			{
				sum += vfbShift / 2;
			}
			else
			{
				sum += vfbShift;
			}

			leftVertID.push_back(lVert->GetID());
			rightVertID.push_back(rVert->GetID());
			flatbandVoltageShift.push_back(norm.PullPotential(sum));
		}

		string timeStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "flatband voltage shift of each interface segment at time=[" + timeStr + "] (left vertex ID, right vertex ID, voltage)";
		file_subs.WriteVector(leftVertID, rightVertID, flatbandVoltageShift, title.c_str());
		if (SctmGlobalControl::Get().CallPytaurus == "Initial" && SctmTimeStep::Get().TimeStep() == 1)
		{
			return;
		}
		file_change.WriteVector(leftVertID, rightVertID, flatbandVoltageShift, title.c_str());

	}

	void SctmData::WriteSubstrateFromInput()
	{
		fileName = directoryName + pathSep + "exchange" + pathSep + "substrate.in";
		SctmFileStream charge_file = SctmFileStream(fileName, SctmFileStream::Read);

		fileName = directoryName + pathSep + "Substrate" + pathSep + "substrate" + generateFileSuffix();
		SctmFileStream subs_file = SctmFileStream(fileName, SctmFileStream::Write);

		vector<int> vertID;
		vector<double> pot;
		vector<double> fermiAbove;

		charge_file.ReadVector(vertID, pot, fermiAbove);
		string timeStr = SctmConverter::DoubleToString(SctmTimeStep::Get().ElapsedTime());
		string title = "substrate result at time=[" + timeStr + "] (vertex ID, potential[V], fermi energy above cb edge [eV])";
		subs_file.WriteVector(vertID, pot, fermiAbove, title.c_str());
	}

	void SctmData::WriteTrappedDensRegionwise(FDDomain *domain)
	{
		//only for triple cell structures
		if (SctmGlobalControl::Get().Structure == "Single")
		{
			return;
		}

		static bool firstRun = true;
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "chargeRegionwise.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		if (firstRun)
		{
			string headline = "Time [s]\t\tTotal Charge\t\tMain Cell [ratio]\t\tOther Regions [ratio]";
			file.WriteLine(headline);
			firstRun = false;
		}

		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		
		double time = SctmTimeStep::Get().ElapsedTime(); //elapsed time, in real value
		double total = 0;
		double mainCell = 0;
		double outLateral = 0;

		double lineDens = 0; //line density (density * area)
		double area = 0;

		FDRegion *region = NULL;
		FDElement *elem = NULL;

		//Trap.Gate2
		region = domain->GetRegion("Trap.Gate2");
		vector<FDElement *> &elems = region->elements;
		for (size_t ielem = 0; ielem != elems.size(); ++ielem)
		{
			elem = elems.at(ielem);
			area = elem->Area / 4;
			lineDens = 0;
			lineDens += elem->SouthwestVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
			lineDens += elem->SoutheastVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
			lineDens += elem->NorthwestVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
			lineDens += elem->NortheastVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
			total += lineDens;
		}
		mainCell = total;

		static bool isLoad = false;
		string simStructure = SctmGlobalControl::Get().Structure;
		vector<string> outRegionsName;

		if (!isLoad && simStructure == "Triple")
		{
			outRegionsName.push_back("Trap.Gate1");
			outRegionsName.push_back("Trap.Iso2");
			outRegionsName.push_back("Trap.Iso3");
			outRegionsName.push_back("Trap.Gate3");
			isLoad = true;
		}
		else if (!isLoad && simStructure == "TripleFull")
		{
			outRegionsName.push_back("Trap.Iso1");
			outRegionsName.push_back("Trap.Gate1");
			outRegionsName.push_back("Trap.Iso2");
			outRegionsName.push_back("Trap.Iso3");
			outRegionsName.push_back("Trap.Gate3");
			outRegionsName.push_back("Trap.Iso4");
			isLoad = true;
		}

		string regName = "";
		for (size_t ir = 0; ir != outRegionsName.size(); ++ir)
		{
			regName = outRegionsName.at(ir);
			region = domain->GetRegion(regName);
			vector<FDElement *> &elems = region->elements;
			for (size_t ielem = 0; ielem != elems.size(); ++ielem)
			{
				elem = elems.at(ielem);
				area = elem->Area / 4;
				lineDens = 0;
				lineDens += elem->SouthwestVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
				lineDens += elem->SoutheastVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
				lineDens += elem->NorthwestVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
				lineDens += elem->NortheastVertex->Trap->GetTrapPrpty(TrapProperty::eTrapped) * area;
				outLateral += lineDens;
				total += lineDens;
			}
		}

		string line = SctmConverter::DoubleToString(time) + "\t\t" +
			SctmConverter::DoubleToString(norm.PullLineDensity(total)) + "\t\t" +
			(total == 0 ? "0" : SctmConverter::DoubleToString(mainCell / total)) + "\t\t" +
			(total == 0 ? "0" : SctmConverter::DoubleToString(outLateral / total));
		file.WriteLine(line);
	}

	void SctmData::ReadTimestep(vector<double> &timestep, vector<double> &vg1, vector<double> &vg2, vector<double> &vg3)
	{
		fileName = directoryName + pathSep + "timestep.in";
		SctmFileStream timestep_file = SctmFileStream(fileName, SctmFileStream::Read);

		if (!SctmFileStream::FileExisted(fileName))
		{
			SCTM_ASSERT(SCTM_ERROR, 10055);
		}
		timestep_file.ReadVector(timestep, vg1, vg2, vg3);
	}

	void SctmData::WriteTunnelOutDensity(FDDomain *domain, VertexMapDouble &tunToSubs, VertexMapDouble &tbToSubs, VertexMapDouble &tunToGate, VertexMapDouble &tbToGate)
	{
		static bool firstRun = true;
		fileName = directoryName + pathSep + "Miscellaneous" + pathSep + "tunnelOutDensity.txt";
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Append);

		if (firstRun)
		{
			string headline = "Time [s]\t\tDT-FN to subs\t\tTrap-to-band to subs\t\tDT-FN to gate\t\tTrap-to-band to Gate";
			file.WriteLine(headline);
			firstRun = false;
		}

		double lineDens_tunToSubs = 0; //line density
		double lineDens_tbToSubs = 0;
		double lineDens_tunToGate = 0;
		double lineDens_tbToGate = 0;

		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		FDVertex *vert = NULL;
		int vertID = 0;
		double controlLength = 0;
		double controlArea = 0;
		double lineDens = 0;
		double currDens = 0;
		double tunCoeff = 0;
		double eDens = 0;
		double tbFrequency = 0;
		double timestep = SctmTimeStep::Get().TimeStep();

		//DN-FN tunneling to substrate
		for (VertexMapDouble::iterator it = tunToSubs.begin(); it != tunToSubs.end(); ++it)
		{
			vertID = it->first;
			vert = domain->GetVertex(vertID);
			controlLength = (vert->EastLength + vert->WestLength) / 2;

			if (vert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelOut)
			{
				tunCoeff = -it->second;
				eDens = vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
				lineDens = eDens * tunCoeff * timestep * controlLength; //charge line density
				lineDens_tunToSubs += lineDens;
			}
			else if (vert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn)
			{
				currDens = it->second;
				lineDens = currDens * controlLength * timestep;
				lineDens_tunToSubs += -lineDens;
			}
		}
		//trap-to-band tunneling to substrate
		for (VertexMapDouble::iterator it = tbToSubs.begin(); it != tbToSubs.end(); ++it)
		{
			vertID = it->first;
			tunCoeff = it->second;
			vert = domain->GetVertex(vertID);

			controlArea = vert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea);
			eDens = vert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
			tbFrequency = vert->Trap->GetTrapPrpty(TrapProperty::eFrequency_T2B);

			lineDens = eDens * tunCoeff * controlArea * timestep * tbFrequency;
			lineDens_tbToSubs += lineDens;
		}
		//DT-FN tunneling to gate
		for (VertexMapDouble::iterator it = tunToGate.begin(); it != tunToGate.end(); ++it)
		{
			vertID = it->first;
			vert = domain->GetVertex(vertID);
			controlLength = (vert->EastLength + vert->WestLength) / 2;

			if (vert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelOut)
			{
				tunCoeff = -it->second;
				eDens = vert->Phys->GetPhysPrpty(PhysProperty::eDensity);
				lineDens = eDens * tunCoeff * timestep * controlLength;
				lineDens_tunToGate += lineDens;
			}
			else if (vert->BndCond.GetBCTunnelTag() == FDBoundary::eTunnelIn)
			{
				currDens = it->second;
				lineDens = currDens * controlLength * timestep;
				lineDens_tunToGate += -lineDens;
			}
		}
		//trap-to-band tunnel to gate
		for (VertexMapDouble::iterator it = tbToGate.begin(); it != tbToGate.end(); ++it)
		{
			vertID = it->first;
			tunCoeff = it->second;
			vert = domain->GetVertex(vertID);

			controlArea = vert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea);
			eDens = vert->Trap->GetTrapPrpty(TrapProperty::eTrapped);
			tbFrequency = vert->Trap->GetTrapPrpty(TrapProperty::eFrequency_T2B);

			lineDens = eDens * tunCoeff * controlArea * timestep * tbFrequency;
			lineDens_tbToGate += lineDens;
		}

		double time = SctmTimeStep::Get().ElapsedTime(); //elapsed time, in real value
		string line = SctmConverter::DoubleToString(time) + "\t\t" +
			SctmConverter::DoubleToString(norm.PullLineDensity(lineDens_tunToSubs)) + "\t\t" +
			SctmConverter::DoubleToString(norm.PullLineDensity(lineDens_tbToSubs)) + "\t\t" +
			SctmConverter::DoubleToString(norm.PullLineDensity(lineDens_tunToGate)) + "\t\t" +
			SctmConverter::DoubleToString(norm.PullLineDensity(lineDens_tbToGate)) + "\t\t";
		file.WriteLine(line);
	}





	string SctmConverter::IntToString(int num)
	{
		stringstream ss;
		ss << num;
		return ss.str();
	}

	string SctmConverter::DoubleToString(double num, bool useScientific /*= true*/, int numAfterPoionts /*= 3*/)
	{
		stringstream ss;
		ss << std::setprecision(numAfterPoionts) << num;
		return ss.str();
	}

	double SctmConverter::StringToDouble(const string &strVal)
	{
		std::istringstream is(strVal);
		double ret = 0;
		if (!(is >> ret))
		{
			throw BadParConversion("cannot convert to double(\"" + strVal + "\")");
		}
		return ret;
	}

	int SctmConverter::StringToInt(const string &strVal)
	{
		std::istringstream is(strVal);
		int ret = 0;
		if (!(is >> ret))
		{
			throw BadParConversion("cannot convert to int(\"" + strVal + "\")");
		}
		return ret;
	}

	bool SctmConverter::StringToBool(const string &strVal)
	{
		if (strVal == "True")
		{
			return true;
		}
		else if (strVal == "False")
		{
			return false;
		}
		else
		{
			throw BadParConversion("cannot convert to bool(\"" + strVal + "\")");
		}
	}






	void SctmGlobalControl::setGlobalCntrl_Directly()
	{
		Get().Temperature = 300;
		Get().GateVoltage = 18;

		Get().GateWorkFunction = 4.7; //for TaN gate

		Get().SimStartTime = 1e-15;
		Get().SimEndTime = 1;
		Get().SimStepsPerDecade = 10;

		//the length are in nm
		Get().XLength = 10;
		Get().YLengthTunnel = 4;
		Get().YLengthTrap = 6.5;
		Get().YLengthBlock = 15;
		Get().XGridNum = 5;
		Get().YGridNumTunnel = 50;
		Get().YGridNumTrap = 100;
		Get().YGridNumBlock = 50;

		Get().TunnelMaterial = Mat::SiO2;
		Get().TrapMaterial = Mat::Si3N4;
		Get().BlockMaterial = Mat::Al2O3;

		Get().UniformTrapDens = 6e19; // in [cm^-3]
		Get().SubstrateDoping = -1e17; // negative for P-type

		//lack of newly-added parameters
	}

	SctmGlobalControl::SctmGlobalControl()
	{
		//setGlobalCntrl_Directly();
		//setGloblaCntrl_FromParFile();
	}

	SctmGlobalControl& SctmGlobalControl::Get()
	{
		static SctmGlobalControl instance;
		return instance;
	}

	void SctmGlobalControl::setGlobalCntrl_FromParFile()
	{
		//set simulation parameters from file
		ParamBase *parBase = NULL;

		//Structure
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::structure);
		Get().Structure = dynamic_cast<Param<string> *>(parBase)->Value();

		//Solver
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::solver);
		Get().Solver = dynamic_cast<Param<string> *>(parBase)->Value();

		//Temperature
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::temperature);
		Get().Temperature = dynamic_cast<Param<double> *>(parBase)->Value();
		//SimStartTime
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_start);
		Get().SimStartTime = dynamic_cast<Param<double> *>(parBase)->Value();
		//SimEndTime
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_end);
		Get().SimEndTime = dynamic_cast<Param<double> *>(parBase)->Value();
		//SimStepsPerDecade
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_stepPerDecade);
		Get().SimStepsPerDecade = dynamic_cast<Param<int> *>(parBase)->Value();
		//SimTimeStepModel
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_stepMode);
		Get().SimTimeStepMode = dynamic_cast<Param<string> *>(parBase)->Value();
		//SimTimeStepScale
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_stepScale);
		Get().SimTimeStepScale = dynamic_cast<Param<string> *>(parBase)->Value();
		//SimTimeStepMax
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::time_stepMax);
		Get().SimTimeStepMax = dynamic_cast<Param<double> *>(parBase)->Value();

		//UniformTrapDens
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::trap_uniDensity);
		Get().UniformTrapDens = dynamic_cast<Param<double> *>(parBase)->Value();
		//SubstrateDoping
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::subs_type);
		string subsType = dynamic_cast<Param<string> *>(parBase)->Value();
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::subs_doping);
		double doping = dynamic_cast<Param<double> *>(parBase)->Value();
		Get().SubstrateDoping = doping * (subsType == "N" ? 1 : -1);
		//TrapDistribution
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::trap_distribution);
		Get().TrapDistribution = dynamic_cast<Param<string> *>(parBase)->Value();

		//TrapCaptureModel
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::trap_capture);
		Get().TrapCaptureModel = dynamic_cast<Param<string> *>(parBase)->Value();
		//PhysicsMFN
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::physics_mfn);
		bool mfn = dynamic_cast<Param<bool> *>(parBase)->Value();
		Get().PhysicsMFN = mfn;
		//PhysicsB2T
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::physics_b2t);
		bool b2t = dynamic_cast<Param<bool> *>(parBase)->Value();
		Get().PhysicsB2T = b2t;
		//PhyscisT2B
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::physics_t2b);
		bool t2b = dynamic_cast<Param<bool> *>(parBase)->Value();
		Get().PhysicsT2B = t2b;
		//PhysicsPFModel
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::physics_pf);
		Get().PhysicsPFModel = dynamic_cast<Param<string> *>(parBase)->Value();

		//for debug
		//TrapOccupation
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_trap_occupy);
		Get().TrapOccupation = dynamic_cast<Param<double> *>(parBase)->Value();
		//TrappedCell
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_trap_cell);
		Get().TrappedCell = dynamic_cast<Param<string> *>(parBase)->Value();
		//LateralTunneling
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_lateral_tunnel);
		Get().LateralTunneling = dynamic_cast<Param<bool> *>(parBase)->Value();
		//CallPytaurus
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_call_pytaurus);
		Get().CallPytaurus = dynamic_cast<Param<string> *>(parBase)->Value();
		//ClearCarrier
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_clear_carrier);
		Get().ClearCarrier = dynamic_cast<Param<bool> *>(parBase)->Value();

		//RetentionAfterPrgrm
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_rAfterP);
		Get().RetentionAfterPrgrm = dynamic_cast<Param<bool> *>(parBase)->Value();
		//RetentionEndTime
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::debug_rEndTime);
		Get().RetentionEndTime = dynamic_cast<Param<double> *>(parBase)->Value();

		//parameters below are used for the construction of SimpleONO domain
		//GateVoltage
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_gate_voltage);
		Get().GateVoltage = dynamic_cast<Param<double> *>(parBase)->Value();
		//GateWorkFunction
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_gate_workfunction);
		Get().GateWorkFunction = dynamic_cast<Param<double> *>(parBase)->Value();

		//Xlength
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_width_value);
		Get().XLength = dynamic_cast<Param<double> *>(parBase)->Value();
		//YLengthTunnel
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_tunnel_thick);
		Get().YLengthTunnel = dynamic_cast<Param<double> *>(parBase)->Value();
		//YLengthTrap
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_trap_thick);
		Get().YLengthTrap = dynamic_cast<Param<double> *>(parBase)->Value();
		//YLengthBlock
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_block_thick);
		Get().YLengthBlock = dynamic_cast<Param<double> *>(parBase)->Value();
		//XGridNum
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_width_grid);
		Get().XGridNum = dynamic_cast<Param<int> *>(parBase)->Value();
		//YGridNumTunnel
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_tunnel_grid);
		Get().YGridNumTunnel = dynamic_cast<Param<int> *>(parBase)->Value();
		//YGridNumTrap
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_trap_grid);
		Get().YGridNumTrap = dynamic_cast<Param<int> *>(parBase)->Value();
		//YGridnumBlock
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_block_grid);
		Get().YGridNumBlock = dynamic_cast<Param<int> *>(parBase)->Value();

		//TunnelMaterial
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_tunnel_material);
		Get().TunnelMaterial = Mat::Parse(dynamic_cast<Param<string> *>(parBase)->Value());
		//TrapMaterial
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_trap_material);
		Get().TrapMaterial = Mat::Parse(dynamic_cast<Param<string> *>(parBase)->Value());
		//BlockMaterial
		parBase = SctmParameterParser::Get().GetPar(SctmParameterParser::sc_block_material);
		Get().BlockMaterial = Mat::Parse(dynamic_cast<Param<string> *>(parBase)->Value());
	}

	void SctmGlobalControl::SetGlobalControl(string defaultParPath, string prjpath)
	{
		SctmGlobalControl::Get().ProjectDirectory = prjpath;
		SctmGlobalControl::Get().DefaulParFile = defaultParPath;
		SctmGlobalControl::Get().UserParFile = prjpath + SctmEnv::Get().PathSep + "user.param";

		setGlobalCntrl_FromParFile();
	}






	SctmParameterParser::SctmParameterParser()
	{
		//string prjPath = SctmGlobalControl::Get().ProjectDirectory;
		defaultParFile = SctmGlobalControl::Get().DefaulParFile;
		userParFile = SctmGlobalControl::Get().UserParFile;

		if (!SctmFileStream::FileExisted(defaultParFile))
		{
			SctmMessaging::Get().PrintFileError(defaultParFile.c_str(), "The default parameter file is missing.");
			SCTM_ASSERT(SCTM_ERROR, 10037);
		}
		ReadParFile(defaultParFile, defaultParMap);
		if (SctmFileStream::FileExisted(userParFile))
		{
			ReadParFile(userParFile, userParMap);
		}
		else
		{
			SctmMessaging::Get().PrintFileError(userParFile.c_str(), "The user parameter file is missing, use default.");
		}
	}

	bool SctmParameterParser::isCommentOrSpaceLine(string &line)
	{
		bool isSpace = (line.find_first_not_of(" \t") == string::npos);
		bool isComment = false;
		if (isSpace)
		{
			return true;
		}
		else
		{
			isComment = trimSpace(line).at(0) == '#';
		}
		return isComment;
	}

	string SctmParameterParser::getParToken(string &line)
	{
		size_t posBegin = line.find_first_not_of(" \t");
		if (posBegin == string::npos)
		{
			//a space line
			return "";
		}
		size_t posEnd = line.find_first_of(":");
		if (posBegin == posEnd)
		{
			//no valid token before :
			return "";
		}
		string effStr = line.substr(posBegin, posEnd - posBegin);
		return trimSpace(effStr);
	}

	string SctmParameterParser::trimSpace(string &line)
	{
		size_t posBegin = line.find_first_not_of(" \t");
		if (posBegin == string::npos)
		{
			//a space line
			return "";
		}
		size_t posEnd = line.find_last_not_of(" \t");
		return line.substr(posBegin, posEnd - posBegin + 1);
	}

	string SctmParameterParser::trimComment(string &line)
	{
		string ret;
		size_t pos = line.find_first_of("#");
		ret = line.substr(0, pos);
		return ret;
	}

	string SctmParameterParser::getParValStr(string &line)
	{
		int posBegin = line.find_first_of(":");
		if (posBegin == line.length() -1)
		{
			//no valid value string after :
			return "";
		}
		string valStr = line.substr(posBegin + 1);
		return trimSpace(valStr);
	}

	bool SctmParameterParser::isValid(string &line)
	{
		size_t pos = line.find_first_of(":");
		return pos != string::npos;
	}

	void SctmParameterParser::parseParValue(std::map<SctmParameterParser::ParName, ParamBase*> &mapToSet, string &name, string &valStr)
	{
		using MaterialDB::Mat;
		static Mat::Name currMat;
		double valDouble = 0;
		bool valBool = false;
		int valInt = 0;

		if (name == "structure")
		{
			ParName pName = ParName::structure;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}
		if (name == "solver")
		{
			ParName pName = ParName::solver;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}
		if (name == "temperature")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::temperature, valDouble);
			mapToSet[ParName::temperature] = par;
			return;
		}
		if (name == "time.start")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::time_start, valDouble);
			mapToSet[ParName::time_start] = par;
			return;
		}
		if (name == "time.end")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::time_end, valDouble);
			mapToSet[ParName::time_end] = par;
			return;
		}
		if (name == "time.stepPerDecade")
		{
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(ParName::time_stepPerDecade, valInt);
			mapToSet[ParName::time_stepPerDecade] = par;
			return;
		}
		if (name == "time.stepMode")
		{
			Param<string> *par = new Param<string>(ParName::time_stepMode, valStr);
			mapToSet[ParName::time_stepMode] = par;
			return;
		}
		if (name == "time.stepScale")
		{
			Param<string> *par = new Param<string>(ParName::time_stepScale, valStr);
			mapToSet[ParName::time_stepScale] = par;
			return;
		}
		if (name == "time.stepMax")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::time_stepMax, valDouble);
			mapToSet[ParName::time_stepMax] = par;
			return;
		}
		if (name == "sc.gate.voltage")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_gate_voltage, valDouble);
			mapToSet[ParName::sc_gate_voltage] = par;
			return;
		}
		if (name == "sc.gate.workfunction")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_gate_workfunction, valDouble);
			mapToSet[ParName::sc_gate_workfunction] = par;
			return;
		}
		if (name == "subs.type")
		{
			Param<string> *par = new Param<string>(ParName::subs_type, valStr);
			mapToSet[ParName::subs_type] = par;
			return;
		}
		if (name == "subs.doping")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::subs_doping, valDouble);
			mapToSet[ParName::subs_doping] = par;
			return;
		}
		if (name == "trap.uniDensity")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::trap_uniDensity, valDouble);
			mapToSet[ParName::trap_uniDensity] = par;
			return;
		}
		if (name == "trap.distribution")
		{
			Param<string> *par = new Param<string>(ParName::trap_distribution, valStr);
			mapToSet[ParName::trap_distribution] = par;
			return;
		}
		if (name == "trap.capture")
		{
			Param<string> *par = new Param<string>(ParName::trap_capture, valStr);
			mapToSet[ParName::trap_capture] = par;
			return;
		}
		if (name == "physics.mfn")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::physics_mfn, valBool);
			mapToSet[ParName::physics_mfn] = par;
			return;
		}
		if (name == "physics.b2t")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::physics_b2t, valBool);
			mapToSet[ParName::physics_b2t] = par;
			return;
		}
		if (name == "physics.t2b")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::physics_t2b, valBool);
			mapToSet[ParName::physics_t2b] = par;
			return;
		}
		if (name == "physics.pf")
		{
			Param<string> *par = new Param<string>(ParName::physics_pf, valStr);
			mapToSet[ParName::physics_pf] = par;
			return;
		}
		if (name == "debug.lateral.tunnel")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::debug_lateral_tunnel, valBool);
			mapToSet[ParName::debug_lateral_tunnel] = par;
			return;
		}
		if (name == "debug.trap.occupy")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::debug_trap_occupy, valDouble);
			mapToSet[ParName::debug_trap_occupy] = par;
			return;
		}
		if (name == "debug.trap.cell")
		{
			Param<string> *par = new Param<string>(ParName::debug_trap_cell, valStr);
			mapToSet[ParName::debug_trap_cell] = par;
			return;
		}
		if (name == "debug.call.pytaurus")
		{
			Param<string> *par = new Param<string>(ParName::debug_call_pytaurus, valStr);
			mapToSet[ParName::debug_call_pytaurus] = par;
			return;
		}
		if (name == "debug.clear.carrier")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::debug_clear_carrier, valBool);
			mapToSet[ParName::debug_clear_carrier] = par;
			return;
		}
		if (name == "debug.rAfterP")
		{
			valBool = SctmConverter::StringToBool(valStr);
			Param<bool> *par = new Param<bool>(ParName::debug_rAfterP, valBool);
			mapToSet[ParName::debug_rAfterP] = par;
			return;
		}
		if (name == "debug.rEndTime")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::debug_rEndTime, valDouble);
			mapToSet[ParName::debug_rEndTime] = par;
			return;
		}
		//parameters for simulation structure
		//single cell
		if (name == "sc.width.value")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_width_value, valDouble);
			mapToSet[ParName::sc_width_value] = par;
			return;
		}
		if (name == "sc.width.grid")
		{
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(ParName::sc_width_grid, valInt);
			mapToSet[ParName::sc_width_grid] = par;
			return;
		}
		if (name == "sc.tunnel.thick")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_tunnel_thick, valDouble);
			mapToSet[ParName::sc_tunnel_thick] = par;
			return;
		}
		if (name == "sc.tunnel.grid")
		{
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(ParName::sc_tunnel_grid, valInt);
			mapToSet[ParName::sc_tunnel_grid] = par;
			return;
		}
		if (name == "sc.tunnel.material")
		{
			Param<string> *par = new Param<string>(ParName::sc_tunnel_material, valStr);
			mapToSet[ParName::sc_tunnel_material] = par;
			return;
		}
		if (name == "sc.trap.thick")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_trap_thick, valDouble);
			mapToSet[ParName::sc_trap_thick] = par;
			return;
		}
		if (name == "sc.trap.grid")
		{
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(ParName::sc_trap_grid, valInt);
			mapToSet[ParName::sc_trap_grid] = par;
			return;
		}
		if (name == "sc.trap.material")
		{
			Param<string> *par = new Param<string>(ParName::sc_trap_material, valStr);
			mapToSet[ParName::sc_trap_material] = par;
			return;
		}
		if (name == "sc.block.thick")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(ParName::sc_block_thick, valDouble);
			mapToSet[ParName::sc_block_thick] = par;
			return;
		}
		if (name == "sc.block.grid")
		{
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(ParName::sc_block_grid, valInt);
			mapToSet[ParName::sc_block_grid] = par;
			return;
		}
		if (name == "sc.block.material")
		{
			Param<string> *par = new Param<string>(ParName::sc_block_material, valStr);
			mapToSet[ParName::sc_block_material] = par;
			return;
		}
		//triple cell
		if (name == "tc.gate1.voltage")
		{
			ParName pName = ParName::tc_gate1_voltage;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate1.workfunction")
		{
			ParName pName = ParName::tc_gate1_workfunction;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate1.width")
		{
			ParName pName = ParName::tc_gate1_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate1.width.grid")
		{
			ParName pName = ParName::tc_gate1_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate2.voltage")
		{
			ParName pName = ParName::tc_gate2_voltage;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate2.workfunction")
		{
			ParName pName = ParName::tc_gate2_workfunction;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate2.width")
		{
			ParName pName = ParName::tc_gate2_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate2.width.grid")
		{
			ParName pName = ParName::tc_gate2_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate3.voltage")
		{
			ParName pName = ParName::tc_gate3_voltage;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate3.workfunction")
		{
			ParName pName = ParName::tc_gate3_workfunction;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate3.width")
		{
			ParName pName = ParName::tc_gate3_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.gate3.width.grid")
		{
			ParName pName = ParName::tc_gate3_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso.material")
		{
			ParName pName = ParName::tc_iso_material;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso.thick")
		{
			ParName pName = ParName::tc_iso_thick;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso.thick.grid")
		{
			ParName pName = ParName::tc_iso_thick_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso1.width")
		{
			ParName pName = ParName::tc_iso1_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso1.width.grid")
		{
			ParName pName = ParName::tc_iso1_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso2.width")
		{
			ParName pName = ParName::tc_iso2_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso2.width.grid")
		{
			ParName pName = ParName::tc_iso2_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso3.width")
		{
			ParName pName = ParName::tc_iso3_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso3.width.grid")
		{
			ParName pName = ParName::tc_iso3_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso4.width")
		{
			ParName pName = ParName::tc_iso4_width;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.iso4.width.grid")
		{
			ParName pName = ParName::tc_iso4_width_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.tunnel.thick")
		{
			ParName pName = ParName::tc_tunnel_thick;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.tunnel.thick.grid")
		{
			ParName pName = ParName::tc_tunnel_thick_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.tunnel.material")
		{
			ParName pName = ParName::tc_tunnel_material;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.trap.thick")
		{
			ParName pName = ParName::tc_trap_thick;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.trap.thick.grid")
		{
			ParName pName = ParName::tc_trap_thick_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.trap.material")
		{
			ParName pName = ParName::tc_trap_material;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.block.thick")
		{
			ParName pName = ParName::tc_block_thick;
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = new Param<double>(pName, valDouble);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.block.thick.grid")
		{
			ParName pName = ParName::tc_block_thick_grid;
			valInt = SctmConverter::StringToInt(valStr);
			Param<int> *par = new Param<int>(pName, valInt);
			mapToSet[pName] = par;
			return;
		}
		if (name == "tc.block.material")
		{
			ParName pName = ParName::tc_block_material;
			Param<string> *par = new Param<string>(pName, valStr);
			mapToSet[pName] = par;
			return;
		}

		//parameters for material properties
		if (name == "material")
		{
			currMat = MaterialDB::Mat::Parse(valStr);
			return;
		}
		if (name == "bandgap")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_bandgap, valDouble);
				mapToSet[ParName::Si_bandgap] = par;
				break;
			case MaterialDB::Mat::SiO2:
				par = new Param<double>(ParName::SiO2_bandgap, valDouble);
				mapToSet[ParName::SiO2_bandgap] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_bandgap, valDouble);
				mapToSet[ParName::Si3N4_bandgap] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_bandgap, valDouble);
				mapToSet[ParName::HfO2_bandgap] = par;
				break;
			case MaterialDB::Mat::Al2O3:
				par = new Param<double>(ParName::Al2O3_bandgap, valDouble);
				mapToSet[ParName::Al2O3_bandgap] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "dielectricConstant")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_dielectricConstant, valDouble);
				mapToSet[ParName::Si_dielectricConstant] = par;
				break;
			case MaterialDB::Mat::SiO2:
				par = new Param<double>(ParName::SiO2_dielectricConstant, valDouble);
				mapToSet[ParName::SiO2_dielectricConstant] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_dielectricConstant, valDouble);
				mapToSet[ParName::Si3N4_dielectricConstant] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_dielectricConstant, valDouble);
				mapToSet[ParName::HfO2_dielectricConstant] = par;
				break;
			case MaterialDB::Mat::Al2O3:
				par = new Param<double>(ParName::Al2O3_dielectricConstant, valDouble);
				mapToSet[ParName::Al2O3_dielectricConstant] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "electronAffinity")
		{
			Param<double> *par = NULL;
			valDouble = SctmConverter::StringToDouble(valStr);
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_electronAffinity, valDouble);
				mapToSet[ParName::Si_electronAffinity] = par;
				break;
			case MaterialDB::Mat::SiO2:
				par = new Param<double>(ParName::SiO2_electronAffinity, valDouble);
				mapToSet[ParName::SiO2_electronAffinity] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_electronAffinity, valDouble);
				mapToSet[ParName::Si3N4_electronAffinity] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_electronAffinity, valDouble);
				mapToSet[ParName::HfO2_electronAffinity] = par;
				break;
			case MaterialDB::Mat::Al2O3:
				par = new Param<double>(ParName::Al2O3_electronAffinity, valDouble);
				mapToSet[ParName::Al2O3_electronAffinity] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eMass")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_eMass, valDouble);
				mapToSet[ParName::Si_eMass] = par;
				break;
			case MaterialDB::Mat::SiO2:
				par = new Param<double>(ParName::SiO2_eMass, valDouble);
				mapToSet[ParName::SiO2_eMass] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eMass, valDouble);
				mapToSet[ParName::Si3N4_eMass] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eMass, valDouble);
				mapToSet[ParName::HfO2_eMass] = par;
				break;
			case MaterialDB::Mat::Al2O3:
				par = new Param<double>(ParName::Al2O3_eMass, valDouble);
				mapToSet[ParName::Al2O3_eMass] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eDOSMass")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_eDOSMass, valDouble);
				mapToSet[ParName::Si_eDOSMass] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eDOSMass, valDouble);
				mapToSet[ParName::Si3N4_eDOSMass] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eDOSMass, valDouble);
				mapToSet[ParName::HfO2_eDOSMass] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "hMass")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_hMass, valDouble);
				mapToSet[ParName::Si_hMass] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eMobility")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Silicon:
				par = new Param<double>(ParName::Si_eMobility, valDouble);
				mapToSet[ParName::Si_eMobility] = par;
				break;
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eMobility, valDouble);
				mapToSet[ParName::Si3N4_eMobility] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eMobility, valDouble);
				mapToSet[ParName::HfO2_eMobility] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eTrapEnergy")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eTrapEnergy, valDouble);
				mapToSet[ParName::Si3N4_eTrapEnergy] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eTrapEnergy, valDouble);
				mapToSet[ParName::HfO2_eTrapEnergy] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eXsection")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eXsection, valDouble);
				mapToSet[ParName::Si3N4_eXsection] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eXsection, valDouble);
				mapToSet[ParName::HfO2_eXsection] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eFrequencyT2B")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eFrequencyT2B, valDouble);
				mapToSet[ParName::Si3N4_eFrequencyT2B] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eFrequencyT2B, valDouble);
				mapToSet[ParName::HfO2_eFrequencyT2B] = par;
				break;
			default:
				break;
			}
			return;
		}
		if (name == "eFrequencyPF")
		{
			valDouble = SctmConverter::StringToDouble(valStr);
			Param<double> *par = NULL;
			switch (currMat)
			{
			case MaterialDB::Mat::Si3N4:
				par = new Param<double>(ParName::Si3N4_eFrequencyPF, valDouble);
				mapToSet[ParName::Si3N4_eFrequencyPF] = par;
				break;
			case MaterialDB::Mat::HfO2:
				par = new Param<double>(ParName::HfO2_eFrequencyPF, valDouble);
				mapToSet[ParName::HfO2_eFrequencyPF] = par;
				break;
			default:
				break;
			}
			return;
		}
		//process the parameters used in Pytaurus
		if (name == "tc.drain.voltage" || name == "tc.drain.voltage.read" || name == "tc.gate.voltage.pass" ||
			name == "tc.gate.voltage.read")
		{
			return;
		}
		//fall to match the above names
		SctmMessaging::Get().PrintInvalidParameterName(name);
		SCTM_ASSERT(SCTM_ERROR, 10047);
	}

	void SctmParameterParser::ReadParFile(string &file, std::map<ParName, ParamBase*> &mapToSet)
	{
		std::ifstream infile(file.c_str());
		string aLine;
		string effline; // remove inline comment
		string parToken;
		string parValue;
		int lineCnt = 0;

		while (std::getline(infile, aLine))
		{
			lineCnt += 1;
			if (isCommentOrSpaceLine(aLine))
			{
				continue;;
			}
			if (!isValid(aLine))
			{
				SctmMessaging::Get().PrintInvalidLineInParFile(file.c_str(), lineCnt);
				SCTM_ASSERT(SCTM_ERROR, 10034);
			}
			effline = trimComment(aLine);
			parToken = getParToken(effline);
			parValue = getParValStr(effline);

			try
			{
				parseParValue(mapToSet, parToken, parValue);
			}
			catch (BadParConversion)
			{
				SctmMessaging::Get().PrintInvalidLineInParFile(file.c_str(), lineCnt);
				SCTM_ASSERT(SCTM_ERROR, 10034);
			}
		}
	}

	SctmParameterParser& SctmParameterParser::Get()
	{
		static SctmParameterParser parser;
		return parser;
	}

	ParamBase * SctmParameterParser::GetPar(ParName name)
	{
		if (userParMap.find(name) != userParMap.end())
		{
			return userParMap[name];
		}
		if (defaultParMap.find(name) != defaultParMap.end())
		{
			return defaultParMap[name];
		}
		else
		{
			//message out
			SCTM_ASSERT(SCTM_ERROR, 10035);
		}
		return NULL;
	}





	SctmEnv& SctmEnv::Get()
	{
		static SctmEnv env;
		return env;
	}

	SctmEnv::SctmEnv()
	{
		if (SCTM_ENV == "Windows")
		{
			DebugPrjPath = "E:\\PhD Study\\SimCTM\\SctmTest\\SolverPackTest";
			DefaultParamPath = "E:\\MyCode\\SimCTM\\SimCTM\\default.param";
			//DefaultParamPath = "E:\\PhD Study\\SimCTM\\default.param";
			ClearPrjPyPath = "E:\\\"PhD Study\"\\SimCTM\\SctmPy\\DeleteData.py";
			PathSep = "\\";
		}
		else if (SCTM_ENV == "Linux")
		{
			DebugPrjPath = "/home/lunzhy/SimCTM/debug";
			DefaultParamPath = "/home/lunzhy/SimCTM/default.param";
			PytaurusPath = "/home/lunzhy/SimCTM/code/Pytaurus/bin/pyt.py";
			PathSep = "/";
		}
	}

	bool SctmEnv::IsLinux()
	{
		return SCTM_ENV == "Linux";
	}

	bool SctmEnv::IsWindows()
	{
		return SCTM_ENV == "Windows";
	}





	void SctmPyCaller::PyClean(string folderPath)
	{
		string command = SctmEnv::Get().PytaurusPath + " clean " + folderPath;
		std::system(command.c_str());
	}

	void SctmPyCaller::PySolve()
	{
		Normalization norm = Normalization(SctmGlobalControl::Get().Temperature);
		string callPytaurusMode = SctmGlobalControl::Get().CallPytaurus;
		string command = "";
		//call pytaurus every time the gate voltage changed
		if (!SctmTimeStep::Get().IsGateVoltageChanged())
		{
			if (callPytaurusMode == "Never")
			{
				return;
			}
			if (callPytaurusMode == "Initial")
			{
				if (SctmTimeStep::Get().StepNumber() != 1) //i.e. first step
				{
					return;
				}
			}
			else if (callPytaurusMode == "EveryStep")
			{
				//do nothing
			}
			else if (callPytaurusMode == "Major")
			{
				if (!SctmTimeStep::Get().IsMajorTime())
				{
					return;
				}
			}
			else
			{
				SCTM_ASSERT(SCTM_ERROR, 10053);
			}
			command = SctmEnv::Get().PytaurusPath + " solve " + SctmGlobalControl::Get().ProjectDirectory;
		}
		else
		{
			command = SctmEnv::Get().PytaurusPath + " solve " + SctmGlobalControl::Get().ProjectDirectory
				+ " -vg1=" + SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellA()))
				+ " -vg2=" + SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellB()))
				+ " -vg3=" + SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellC()));
		}
		
		string header = "Solve channel potential using Pytaurus. Vgate1 = " + 
			SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellA()))
			+ " Vgate2 = " + SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellB()))
			+ " Vgate3 = " + SctmConverter::DoubleToString(norm.PullPotential(SctmTimeStep::Get().VoltageCellC()));
		SctmMessaging::Get().PrintHeader(header.c_str());
		
		std::system(command.c_str());
	}

	void SctmPyCaller::PyPrepare(string folderPath)
	{
		string command = SctmEnv::Get().PytaurusPath + " prepare " + folderPath;
		std::system(command.c_str());
	}

	void SctmPyCaller::PyBuildStructure()
	{
		string callPytaurusMode = SctmGlobalControl::Get().CallPytaurus;
		if (callPytaurusMode == "Never")
		{
			return;
		}

		SctmMessaging::Get().PrintHeader("Build structure using Pytaurus.");
		string command = SctmEnv::Get().PytaurusPath + " structure " +
			SctmGlobalControl::Get().ProjectDirectory;
		std::system(command.c_str());
	}

	void SctmPyCaller::PyParseAvgVfb()
	{
		string command = SctmEnv::Get().PytaurusPath + " parsevth " + SctmGlobalControl::Get().ProjectDirectory;
		std::system(command.c_str());
	}





}
