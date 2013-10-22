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

using std::cout;
using std::endl;
using std::fstream;

namespace SctmUtils
{
	SctmMessaging UtilsMsg = SctmMessaging();
	SctmTimer UtilsTimer = SctmTimer();
	SctmDebug UtilsDebug = SctmDebug();
	SctmTimeStep UtilsTimeStep = SctmTimeStep();

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
		default:
			msg = "Untracked error";
		}
		PrintErrorInfo(msg);
	}

	void SctmDebug::PrintDomainDetails(FDDomain &domain)
	{
		if (!this->enable)
			return;
		using namespace MaterialDB;
		FDVertex *currVert = NULL;
		Normalization norm = Normalization();
		for (size_t iVert = 0; iVert != domain.vertices.size(); ++iVert)
		{
			currVert = domain.GetVertex(iVert);
			PrintValue(currVert->GetID());
			cout << " -- ";
			PrintValue(currVert->IsAtContact());
			PrintValue(currVert->IsAtBoundary(FDBoundary::eCurrentDensity));
			PrintValue(currVert->BndCond.Valid(FDBoundary::eCurrentDensity));
			if (currVert->BndCond.Valid(FDBoundary::eCurrentDensity)) { PrintBCType(currVert->BndCond); }
			//PrintValue(currVert->EastLength);
			//PrintValue(currVert->WestLength);
			//PrintValue(currVert->SouthLength);
			//PrintValue(currVert->NorthLength);
			//PrintValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->GetID());
			//PrintValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->GetID());
			//PrintValue(currVert->SouthVertex==NULL ? -1 : currVert->SouthVertex->GetID());
			//PrintValue(currVert->NorthVertex==NULL ? -1 : currVert->NorthVertex->GetID());
			cout << " -- ";
			PrintValue(norm.PullLength(currVert->EastLength));
			PrintValue(norm.PullLength(currVert->WestLength));
			PrintValue(norm.PullLength(currVert->SouthLength));
			PrintValue(norm.PullLength(currVert->NorthLength));
			cout << " -- ";
			//PrintValue(currVert->NorthwestElem==NULL ? -1 : currVert->NorthwestElem->GetInternalID());
			//PrintValue(currVert->NortheastElem==NULL ? -1 : currVert->NortheastElem->GetInternalID());
			//PrintValue(currVert->SouthwestElem==NULL ? -1 : currVert->SouthwestElem->GetInternalID());
			//PrintValue(currVert->SoutheastElem==NULL ? -1 : currVert->SoutheastElem->GetInternalID());
			cout << " -- ";
			//PrintValue(currVert->Phys.GetPhysPrpty(PhysProperty::eMobility));
			//PrintValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->SouthVertex==NULL ? -1 : currVert->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			//PrintValue(currVert->NorthVertex==NULL ? -1 : currVert->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			cout << " -- ";
			//PrintValue(currVert->NorthwestElem==NULL ? -1 : GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			//PrintValue(currVert->NortheastElem==NULL ? -1 : GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_Bandgap));;
			//PrintValue(currVert->SouthwestElem==NULL ? -1 : GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			//PrintValue(currVert->SoutheastElem==NULL ? -1 : GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_Bandgap));
			cout << endl;
		}
	}

	void SctmDebug::PrintBCType(FDBoundary &bc)
	{
		if (!this->enable)
			return;
		string typestring;
		switch (bc.GetBCType(FDBoundary::eCurrentDensity))
		{
		case FDBoundary::BC_Dirichlet:
			typestring = "Dirichlet";
			break;
		case FDBoundary::BC_Neumann:
			typestring = "Neumann";
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


	SctmFileOperator::SctmFileOperator(string _filename, FileMode _mode)
	{
		this->fileName = _filename;
		//if the file doesn't exist, create it.
		fstream file;
		if (_mode == Write)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file)
				file.open(this->fileName.c_str(), std::ios::out);
			else
			{
				file.close();//this is important because the existed file has been opened already.
				file.open(this->fileName.c_str(), std::ios::out | std::ios::trunc);
				file.close();
			}
		}
		if (_mode = Read)
		{
			file.open(this->fileName.c_str(), std::ios::in);
			if (!file.is_open())
				UtilsMsg.PrintFileError(_filename.c_str());
		}
	}

	void SctmFileOperator::Write2DVectorForOrigin(vector<double> &vecX, vector<double> &vecY, vector<vector<double>> &vector2D, const char *title)
	{
		fstream tofile;
		tofile.open(this->fileName.c_str(), std::ios::app);
		if (!tofile) 
			UtilsMsg.PrintFileError(this->fileName.c_str());
		
		tofile << title << endl;
		for (size_t ix = 0; ix != vecX.size(); ++ix)
		{
			for (size_t iy = 0; iy != vecY.size(); ++iy)
			{
				tofile << vecX.at(ix) << '\t' << vecY.at(iy) << '\t' << vector2D.at(ix).at(iy) << endl;
			}
		}
		tofile.close();
	}

	void SctmFileOperator::WriteVector(vector<double> &vec, const char *title)
	{
		fstream tofile;
		tofile.open(this->fileName.c_str(), std::ios::app);
		if (!tofile) 
			UtilsMsg.PrintFileError(this->fileName.c_str());

		tofile << title << endl;
		for (size_t ix = 0; ix != vec.size(); ++ix)
		{
			tofile << vec.at(ix) <<  '\t';
		}
		tofile << endl;
		tofile.close();
	}

	void SctmFileOperator::ReadTunnelParameter(vector<double> &cbedges, vector<double> &elecfields)
	{
		cbedges.clear(); elecfields.clear();
		std::ifstream file(this->fileName.c_str());

		double val = 0;
		while (!file.eof())
		{
			file >> val; elecfields.push_back(val);
			file >> val; cbedges.push_back(val);
		}
		file.close();
	}

	void SctmFileOperator::WriteDDResult(vector<FDVertex *> &vertices, const char *title)
	{
		fstream tofile;
		tofile.open(this->fileName.c_str(), std::ios::app);
		
		Normalization norm = Normalization();
		FDVertex *currVert = NULL;
		tofile << title << endl;
		for (size_t iVert = 0; iVert != vertices.size(); ++iVert)
		{
			currVert = vertices.at(iVert);
			tofile << norm.PullLength(currVert->X) << '\t' << norm.PullLength(currVert->Y) << '\t' << 
				norm.PullDensity(currVert->Phys.GetPhysPrpty(PhysProperty::eDensity)) << endl;
		}
		tofile.close();
	}

	SctmTimeStep::SctmTimeStep()
	{
		Normalization norm = Normalization();
		this->timeNormFactor = norm.timeFactor;
	}

	double SctmTimeStep::NextTimeStep()
	{
		//TODO: currently, constant time step is used in the simulation
		double timestep = 1e-11; // in [s]
		return timestep / timeNormFactor;
	}

}