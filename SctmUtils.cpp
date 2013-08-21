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

#include <iostream>
using std::cout;
using std::endl;

namespace SctmUtils
{
	SctmMessaging Msg = SctmMessaging();

	void SctmTimer::Start()
	{
		start_time = clock();
		return;
	}

	void SctmTimer::Reset()
	{
		start_time = 0;
		end_time = 0;
		duration = 0;
		Start();
		return;
	}

	void SctmTimer::End()
	{
		end_time = clock();
		duration = (end_time - start_time) / clockPerSecond;
		return;
	}

	double SctmTimer::Duration()
	{
		End();
		return duration;
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
			msg	= "[DomainDetails.cpp] Could not find the boundary condition name.";
			break;
		default:
			msg = "Untracked error";
		}
		PrintErrorInfo(msg);
	}

	void SctmDebug::PrintDomainDetails(FDDomain &domain)
	{
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
		cout << matrix << endl;
	}

	void SctmDebug::PrintVector(const std::vector<double> &vec)
	{
		for (size_t iVec = 0; iVec != vec.size(); ++iVec)
		{
			printValue(vec.at(iVec));
		}
		cout << endl;
	}

	void SctmMessaging::printMessage(string msg)
	{
		std::cout << msg << std::endl;
	}

	void SctmMessaging::PrintHeader(string header)
	{
		string msg = "----------------------------" + header + "----------------------------";
		printMessage(msg);
	}
}

