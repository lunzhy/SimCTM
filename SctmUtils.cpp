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
		case 1:
			msg = "Non-existed physical property.";
			break;
		case 2:
			msg = "Non-existed material property.";
			break;
		case 3:
			msg = "Not rectangular element in constructing elements.";
			break;
		case 4:
			msg = "Error in calculating vertex-related physical value using material-base property.";
			break;
		case 5:
			msg = "[MatrixSolver.cpp] The size of right-hand side vector and solution vector are not equal in matrix solver.";
			break;
		case 6:
			msg = "[MatrixSolver.cpp] The solver of SparseLU fails to solver matrix equation.";
			break;
		default:
			msg = "Untracked error";
		}
		Message(msg);
	}
}

