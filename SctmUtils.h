/**
* @file SctmUtils.h contains the utilities used in the simulation
* @brief
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
#include <stdexcept>
#include <string>
#include <ctime>

#define SCTM_CHECK(cond, err_code) if (!(cond)) { throw std::runtime_error(SctmUtils::ErrorCodeParser(err_code)); }

namespace SctmUtils
{
	std::string ErrorCodeParser(int err_code)
	{
		std::string msg;
		switch (err_code)
		{
			
		}
		return msg;
	}

	class SctmTimer
	{
		static const int clockPerSecond = CLOCKS_PER_SEC;
	public:
		SctmTimer():start_time(0), end_time(0), duration(0)
		{
			Start();
		}
		void Start()
		{
			start_time = clock();
			return;
		}
		void Reset()
		{
			start_time = 0;
			end_time = 0;
			duration = 0;
			Start();
			return;
		}
		void End()
		{
			end_time = clock();
			duration = (end_time - start_time) / clockPerSecond;
			return;
		}
		double Duration()
		{
			End();
			return duration;
		}
	protected:
		double duration;
		std::clock_t start_time;
		std::clock_t end_time;

	};
}