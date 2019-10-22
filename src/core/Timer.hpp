#ifndef CORE_TIMER_HPP
#define CORE_TIMER_HPP

#include <chrono>
#include <string>

class Timer
{
public:
	Timer(const std::string& _phrase);
	~Timer();
private:
	std::chrono::high_resolution_clock::time_point _i;
	std::string _phrase;
};

#endif // CORE_TIMER_HPP
