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
			msg = "Non-existed physical property";
		case 2:
			msg = "Non-existed material property";
		case 3:
			msg = "Not rectangular element in constructing elements";
		default:
			msg = "Untracked error";
		}
		Message(msg);
	}
}

