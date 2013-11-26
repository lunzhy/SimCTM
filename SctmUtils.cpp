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
#include "Material.h"
#include "Normalization.h"
#include <iomanip>
#include <cstring>
#include <fstream>
#include "SctmPhys.h"

using std::cout;
using std::endl;
using std::fstream;
using std::stringstream;
using SctmPhys::PhysProperty;

namespace SctmUtils
{
	SctmMessaging UtilsMsg = SctmMessaging();
	SctmTimer UtilsTimer = SctmTimer();
	SctmDebug UtilsDebug = SctmDebug();
	SctmTimeStep UtilsTimeStep = SctmTimeStep();
	SctmData UtilsData = SctmData();

	void SctmTimer::Start()
	{
		start_time = clock();
		set_time = 0;
		end_time = 0;
		return;
	}

	void SctmTimer::Set()
	{
		set_time = clock();
		return;
	}

	void SctmTimer::End()
	{
		end_time = clock();
		return;
	}

	double SctmTimer::SinceLastSet()
	{
		double time = 0;
		clock_t current_time = clock();
		//TODO: is this conversion correct?
		time = (double)((current_time - set_time) / (clock_t)clockPerSecond);
		return time;
	}

	double SctmTimer::TotalTime()
	{
		double time = 0;
		time = (double)((end_time - start_time ) / (clock_t)clockPerSecond);
		return time;
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
			msg = "[SctmPhys.cpp] Non-existed physical property.";
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
			msg = "[DDSolver.cpp] Errors occurred when passing the boundary current density from tunneling layer to ddsolver.";
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
		Normalization norm = Normalization();
		for (size_t iVert = 0; iVert != domain->vertices.size(); ++iVert)
		{
			currVert = domain->GetVertex(iVert);
			PrintValue(currVert->GetID());
			cout << " -- ";
			//PrintValue(currVert->IsAtContact());
			//PrintValue(currVert->IsAtBoundary(FDBoundary::eCurrentDensity));
			//PrintValue(currVert->BndCond.Valid(FDBoundary::eCurrentDensity));
			//if (currVert->BndCond.Valid(FDBoundary::eCurrentDensity)) { PrintBCType(currVert->BndCond); }
			PrintValue(currVert->IsAtBoundary(FDBoundary::Potential));
			PrintValue(currVert->BndCond.Valid(FDBoundary::Potential));
			if (currVert->IsAtBoundary(FDBoundary::Potential))
			{
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
			//if (currVert->BndCond.Valid(FDBoundary::eDensity)) { PrintDirectionVector(currVert->BndCond.GetBCNormVector(FDBoundary::eDensity)); }
			//PrintValue(currVert->EastLength);
			//PrintValue(currVert->WestLength);
			//PrintValue(currVert->SouthLength);
			//PrintValue(currVert->NorthLength);
			//PrintValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->GetID());
			//PrintValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->GetID());
			//PrintValue(currVert->SouthVertex==NULL ? -1 : currVert->SouthVertex->GetID());
			//PrintValue(currVert->NorthVertex==NULL ? -1 : currVert->NorthVertex->GetID());
			//cout << " -- ";
			//PrintValue(norm.PullLength(currVert->EastLength));
			//PrintValue(norm.PullLength(currVert->WestLength));
			//PrintValue(norm.PullLength(currVert->SouthLength));
			//PrintValue(norm.PullLength(currVert->NorthLength));
			//cout << " -- ";
			//PrintValue(currVert->NorthwestElem==NULL ? -1 : currVert->NorthwestElem->GetInternalID());
			//PrintValue(currVert->NortheastElem==NULL ? -1 : currVert->NortheastElem->GetInternalID());
			//PrintValue(currVert->SouthwestElem==NULL ? -1 : currVert->SouthwestElem->GetInternalID());
			//PrintValue(currVert->SoutheastElem==NULL ? -1 : currVert->SoutheastElem->GetInternalID());
			cout << " -- ";
			PrintValue(currVert->Phys->GetPhysPrpty(PhysProperty::DensityControlArea));
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
		cout << dv;
	}

	void SctmDebug::WritePoisson(FDDomain *domain)
	{
		UtilsData.WritePotential(domain->GetVertices());
	}

	void SctmDebug::WriteBandInfo(FDDomain *domain)
	{
		UtilsData.WriteBandInfo(domain->GetVertices());
	}

	void SctmDebug::WriteDensity(FDDomain *domain)
	{
		UtilsData.WriteElecDens(domain->GetDDVerts());
	}



	void SctmMessaging::printLine(string &line)
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
		string timeStr = ConvertToString::Double(UtilsTimeStep.ElapsedTime());
		cout << timeStr << "s" << "\t\t";

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

	void SctmMessaging::PrintFileError(const char *filename)
	{
		cout << endl << "XXXXX=>" << "cannot open or create file" << ' ' << filename << endl;
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
					UtilsMsg.PrintDirectoryError();
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
		}
		if (_mode == Read)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file.is_open())
				UtilsMsg.PrintFileError(_filename.c_str());
			else
				file.close();
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



	SctmTimeStep::SctmTimeStep()
	{
		this->currStepNumber = 0;
		this->currElapsedTime = 0;
	}

	void SctmTimeStep::GenerateNext()
	{
		currStepNumber += 1; // step number starts with 1
		this->currTimeStep = getTimeStep();
		this->currElapsedTime += this->currTimeStep;
	}

	double SctmTimeStep::TimeStep() const
	{
		return currTimeStep;
	}

	double SctmTimeStep::ElapsedTime() const
	{
		Normalization norm = Normalization();
		return norm.PullTime(currElapsedTime);
	}

	double SctmTimeStep::getTimeStep()
	{
		//TODO: currently, constant time step is used in the simulation
		Normalization norm = Normalization();
		double next = 0; // in [s]

		while(false)
		{
			if (currStepNumber <= 10)
			{
				next = 1e-15;
				break;
			}
			if (currStepNumber <= 19)
			{
				next = 1e-14;
				break;
			}
			if (currStepNumber <= 28 )
			{
				next = 1e-13;
				break;
			}
			if (currStepNumber <= 37)
			{
				next = 1e-12;
				break;
			}
			if (currStepNumber <= 46)
			{
				next = 1e-11;
				break;
			}
			break;
		}

		while(false)
		{
			if (currStepNumber <= 10)
			{
				next = 1e-15;
				break;
			}
			if (currStepNumber <= 13)
			{
				next = 3e-14;
				break;
			}
			if (currStepNumber <= 16 )
			{
				next = 3e-13;
				break;
			}
			if (currStepNumber <= 19)
			{
				next = 3e-12;
				break;
			}
			if (currStepNumber <= 22)
			{
				next = 3e-11;
				break;
			}
			break;
		}
		
		next = 1e-5;
		return norm.PushTime(next);
	}

	int SctmTimeStep::StepNumber() const
	{
		return currStepNumber;
	}

	SctmData::SctmData()
	{
		directoryName = "E:\\PhD Study\\SimCTM\\SctmTest\\SolverPackTest\\";
	}


	void SctmData::WriteElecDens(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + "eDensity" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization();
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

		string numStr = ConvertToString::Double(UtilsTimeStep.ElapsedTime());
		string title = "electron density of time [" + numStr + "]"; 
		file.WriteVector(vecX, vecY, vecDen, title.c_str());
	}

	void SctmData::WritePotential(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + "potential" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization();
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

		string numStr = ConvertToString::Double(UtilsTimeStep.ElapsedTime());
		string title = "potential of time [" + numStr + "]"; 
		file.WriteVector(vecX, vecY, vecPot, title.c_str());
	}

	string SctmData::generateFileSuffix()
	{
		string ret = "";
		string step = ConvertToString::Int(UtilsTimeStep.StepNumber());

		ret = "_s" + step + ".txt";
		return ret;
	}

	void SctmData::WriteBandInfo(vector<FDVertex *> &vertices)
	{
		fileName = directoryName + "band" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);

		Normalization norm = Normalization();
		vector<double> vecX;
		vector<double> vecY;
		vector<double> vecCB;
		vector<double> vecVB;

		FDVertex *currVert = NULL;
		for (size_t iVer = 0; iVer != vertices.size(); ++iVer)
		{
			currVert = vertices.at(iVer);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			vecCB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ConductionBandEnergy)));
			vecVB.push_back(norm.PullEnergy(currVert->Phys->GetPhysPrpty(PhysProperty::ValenceBandEnergy)));
		}

		string numStr = ConvertToString::Double(UtilsTimeStep.ElapsedTime());
		string title = "band structure of time [" + numStr + "] (x, y, Conduction band, Valence band)"; 
		file.WriteVector(vecX, vecY, vecCB, vecVB, title.c_str());
	}

	void SctmData::WriteTunnelCurrentFromSubs(FDDomain *domain, VertexMapDouble &currDensity)
	{
		fileName = directoryName + "subsCurrent" + generateFileSuffix();
		SctmFileStream file = SctmFileStream(fileName, SctmFileStream::Write);
		
		int vertID = 0;
		vector<double> vecX;
		vector<double> vecY;
		vector<double> currDens;
		Normalization norm = Normalization();
		FDVertex *currVert = NULL;

		for (VertexMapDouble::iterator it = currDensity.begin(); it != currDensity.end(); ++it)
		{
			vertID = it->first;
			currVert = domain->GetVertex(vertID);
			vecX.push_back(norm.PullLength(currVert->X));
			vecY.push_back(norm.PullLength(currVert->Y));
			currDens.push_back(norm.PullCurrDens(it->second));
		}
		string numStr = ConvertToString::Double(UtilsTimeStep.ElapsedTime());
		string title = "band structure of time [" + numStr + "] (x, y, tunneling current density)";
		file.WriteVector(vecX, vecY, currDens, title.c_str());
	}


	string ConvertToString::Int(int num)
	{
		stringstream ss;
		ss << num;
		return ss.str();
	}

	string ConvertToString::Double(double num, bool useScientific /*= true*/, int numAfterPoionts /*= 3*/)
	{
		stringstream ss;
		ss << num;
		return ss.str();
	}


}