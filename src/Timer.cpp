#include "DebugMacros.h"
#include "Timer.h"

Timer::Timer(const std::string& phrase)
                : _i{std::chrono::high_resolution_clock::now()}, _phrase{phrase}
{
}

Timer::~Timer()
{
	auto _f = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(_f - _i);
	DEBUG(_phrase << " completed in " << duration.count() << " ms." << std::endl);
}
