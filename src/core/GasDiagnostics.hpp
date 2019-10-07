#ifndef CORE_GASDIAGNOSTICS_H_
#define CORE_GASDIAGNOSTICS_H_

#include "Array.hpp"

#include <map>
#include <vector>

/** This class should be used to store any extra information that is not put into the standard
    GasState, to avoid making it too bloaty. Maybe it will make sense to design it so that not
    all the members have to be set. The amount of detail can maybe be set using some sort of
    configuration file. */
class GasDiagnostics
{
public:
	/** Configuration */
	bool saveLevelPopulations() const { return _saveLevelPopulations; }

	/** The photoelectric heating contribution by each grain population */
	Array photoelectricHeating() const { return _photoelectricHeating; }

	/** Other heating rates, stored as a map with keys */
	const std::map<std::string, double> otherHeating() const { return _heatingm; }
	const std::map<std::string, double> cooling() const { return _coolingm; }

	/** Names of the reactions in the chemical network used */
	std::vector<std::string> reactionNames() const { return _reactionNamev; }
	/** Rates of the reactions in the chemical network used */
	Array reactionRates() const { return _reactionRatev; }

	/** Population density of the included H levels */
	Array hPopulations() const { return _hPopulationv; }
	/** Population density of the included H2 levels */
	Array h2Populations() const { return _h2Populationv; }

	/** Get the value of a diagnostic that was added using putUserValue */
	double userValue(const std::string& key) { return _userValuem.at(key); }

	void setPhotoelectricHeating(const Array& peh) { _photoelectricHeating = peh; }
	void setReactionNames(const std::vector<std::string>& rn) { _reactionNamev = rn; }
	void setReactionRates(const Array& rr) { _reactionRatev = rr; }
	void setHPopulations(const Array& hp) { _hPopulationv = hp; }
	void setH2Populations(const Array& h2p) { _h2Populationv = h2p; }
	void setHeating(const std::string& key, double value);
	void setCooling(const std::string& key, double value);
	void setUserValue(const std::string& key, double value);

private:
	bool _saveLevelPopulations{true};

	Array _photoelectricHeating;

	std::map<std::string, double> _heatingm;
	std::map<std::string, double> _coolingm;

	std::vector<std::string> _reactionNamev;
	Array _reactionRatev;

	Array _hPopulationv;
	Array _h2Populationv;

	std::map<std::string, double> _userValuem;
};

#endif /* CORE_GASDIAGNOSTICS_H_ */
