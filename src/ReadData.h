#ifndef _READDATA_H_
#define _READDATA_H_

#include <vector>
#include <string>

namespace ReadData
{

void recombinationContinuum(std::string file, std::vector<double>& fileFrequencyv,
		std::vector<double>& fileThresholdv, std::vector<double>& fileTemperaturev,
		std::vector<std::vector<double>>& fileGammaDaggervv);

}
/* namespace ReadData */

#endif /* _READDATA_H_ */
