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
#include <iostream>
using std::cout;
using std::endl;

namespace SctmUtils
{
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

	void SctmDebug::Message(string msg)
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
			msg = "Non-existed physical property.";
			break;
		case 10002:
			msg = "Non-existed material property.";
			break;
		case 10003:
			msg = "Not rectangular element in constructing elements.";
			break;
		case 10004:
			msg = "Error in calculating vertex-related physical value using material-base property.";
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
		Message(msg);
	}
}

