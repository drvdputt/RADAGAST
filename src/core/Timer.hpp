#ifndef CORE_TIMER_H_
#define CORE_TIMER_H_

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

#endif /* CORE_TIMER_H_ */
