#ifndef CORE_SPECIESMODELMANAGER_HPP
#define CORE_SPECIESMODELMANAGER_HPP

#include "H2Data.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"

class HData;
class H2Data;

/** This class is responsible for loading and storing the right data for the specialized species
    models (H and H2), and (sometimes polymorphically) constructing such models while giving
    them access to the necessary data through a pointer. */
class SpeciesModelManager
{
public:
	/** Create an instance of SpeciesModelManager with specific settings for H and H2.
	    Depending the given options, different data sets might be loaded, and the behaviour
	    of the make functions will be influenced. The options are undocumented for now, as
	    the are subject to change. Look at the 'if's in the implementation. */
	SpeciesModelManager(const std::string& hOption, const std::string& h2Option);

	/** Create a new HModel, which is given a pointer to the H data stored here. */
	std::unique_ptr<HModel> makeHModel() const;

	/** Create a new H2Model polymorphically. Depending on the heuristics given, either a
	    BigH2Model or a SmallH2Model might be created. A BigH2Model will get access to the
	    H2 data stored here, while a SmallH2Model ideally does not need it. */
	std::unique_ptr<H2Model> makeH2Model(/* heuristic parameters */) const;

private:
	std::unique_ptr<HData> _hData;
	std::unique_ptr<H2Data> _h2Data;
};

#endif // CORE_SPECIESMODELMANAGER_HPP
