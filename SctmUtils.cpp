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
#include <iomanip>
#include <cstring>

using std::cout;
using std::endl;

namespace SctmUtils
{
	SctmMessaging UtilsMsg = SctmMessaging();
	SctmTimer UtilsTimer = SctmTimer();
	SctmDebug UtilsDebug = SctmDebug();

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
		for (size_t iVert = 0; iVert != domain.vertices.size(); ++iVert)
		{
			currVert = domain.GetVertex(iVert);
			printValue(currVert->GetInternalID()); cout << " -- ";
			printValue(currVert->IsAtContact());
			printValue(currVert->IsAtBoundary());
			if (currVert->BndCond.Valid()) { printBCType(currVert->BndCond); }
			printValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->GetInternalID());
			printValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->GetInternalID());
			printValue(currVert->SouthVertex==NULL ? -1 : currVert->SouthVertex->GetInternalID());
			printValue(currVert->NorthVertex==NULL ? -1 : currVert->NorthVertex->GetInternalID()); cout << " -- ";
			printValue(currVert->EastLength);
			printValue(currVert->WestLength);
			printValue(currVert->SouthLength);
			printValue(currVert->NorthLength); cout << " -- ";
			printValue(currVert->NorthwestElem==NULL ? -1 : currVert->NorthwestElem->GetInternalID());
			printValue(currVert->NortheastElem==NULL ? -1 : currVert->NortheastElem->GetInternalID());
			printValue(currVert->SouthwestElem==NULL ? -1 : currVert->SouthwestElem->GetInternalID());
			printValue(currVert->SoutheastElem==NULL ? -1 : currVert->SoutheastElem->GetInternalID()); cout << " -- ";
			printValue(currVert->EastVertex==NULL ? -1 : currVert->EastVertex->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			printValue(currVert->WestVertex==NULL ? -1 : currVert->WestVertex->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			printValue(currVert->SouthVertex==NULL ? -1 : currVert->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity));
			printValue(currVert->NorthVertex==NULL ? -1 : currVert->Phys.GetPhysPrpty(PhysProperty::ElectronAffinity)); cout << " -- ";
			printValue(currVert->NorthwestElem==NULL ? -1 : GetMatPrpty(currVert->NorthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			printValue(currVert->NortheastElem==NULL ? -1 : GetMatPrpty(currVert->NortheastElem->Region->Mat, MatProperty::Mat_Bandgap));;
			printValue(currVert->SouthwestElem==NULL ? -1 : GetMatPrpty(currVert->SouthwestElem->Region->Mat, MatProperty::Mat_Bandgap));
			printValue(currVert->SoutheastElem==NULL ? -1 : GetMatPrpty(currVert->SoutheastElem->Region->Mat, MatProperty::Mat_Bandgap));
			cout << endl;
		}
	}

	void SctmDebug::printBCType(FDBoundary &bctype)
	{
		if (!this->enable)
			return;
		string typestring;
		switch (bctype.GetBCType(FDBoundary::Potential))
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
		printValue(typestring);
	}

	void SctmDebug::PrintSparseMatrix(Eigen::SparseMatrix<double> &matrix)
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
		cout << msg << time << "s" << endl;
		cout.unsetf(std::ios::fixed);
		cout.precision(6);
		return;
	}

	void SctmMessaging::PrintHeader(const char *header)
	{
		cout << "===> " << header << endl << endl;
		return;
	}
}

