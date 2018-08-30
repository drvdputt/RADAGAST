/** This file contains a main which will run all tests. Or something. */

#include "LineProfile.h"
#include "Testing.h"

class Test
{
public:
	Test(std::function<void()> run, std::string name) : _run{run}, _name{name} {};

	void run() const
	{
		std::cout << "-----------------------------------------------------------------"
		             "---------------\n";
		std::cout << "running " << _name << std::endl;
		_run();
		std::cout << "-----------------------------------------------------------------"
		             "---------------"
		          << std::endl;
	}

private:
	std::function<void()> _run;
	std::string _name;
};

int main()
{
	std::vector<Test> tv = {{test_addToSpectrum, "LineProfile addToSpectrum"},
			     {test_integrateSpectrum, "LineProfile integrateSpectrum"},
			     {Testing::testChemistry, "chemistry"},
			     {Testing::testACollapse, "A matrix collapse"},
			     {Testing::testFromFilesvsHardCoded, "HydrogenFromFiles vs HydrogenHardcoded"}
	};
	for (const Test& t : tv)
		t.run();
}
