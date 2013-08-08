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
#pragma once

#define DEBUG
#include <string>
#include <ctime>

#ifdef DEBUG
	#define SCTM_ASSERT(cond, err_code) if (!(cond)) { SctmUtils::ErrorCodeParser(err_code); }
#else
	#define SCTM_ASSERT(cond, err_code)
#endif // DEBUG


using std::string;
namespace SctmUtils
{
	void Message(string msg);
	void ErrorCodeParser(int err_code);

	class SctmTimer
	{
		static const int clockPerSecond = CLOCKS_PER_SEC;
	public:
		SctmTimer():start_time(0), end_time(0), duration(0)
		{
			Start();
		}
		void Start();
		void Reset();
		void End();
		double Duration();

	protected:
		double duration;
		std::clock_t start_time;
		std::clock_t end_time;
	};
}