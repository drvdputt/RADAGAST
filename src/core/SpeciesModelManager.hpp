#ifndef CORE_SPECIESMODELMANAGER_HPP
#define CORE_SPECIESMODELMANAGER_HPP

class HData;

class SpeciesModelManager
{
public:
	SpeciesModelManager(/* h settings, h2 settings */);
	HModel makeHModel() const;
	H2Model makeH2Model(/* heuristic parameters */) const;

private:
	unique_ptr<HData> _hData;
	unique_ptr<H2Data> _h2Data;
}

#endif // CORE_SPECIESMODELMANAGER_HPP
