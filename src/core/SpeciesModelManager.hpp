#ifndef CORE_SPECIESMODELMANAGER_HPP
#define CORE_SPECIESMODELMANAGER_HPP

#include "H2Data.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"
#include "LookupTable.hpp"
#include "SimpleColumnFile.hpp"

namespace GasModule
{
    class HData;
    class H2Data;

    /** This class is responsible for loading and storing the right data for the specialized
        species models (H and H2), and (sometimes polymorphically) constructing such models while
        giving them access to the necessary data through a pointer. */
    class SpeciesModelManager
    {
    public:
        /** Create an instance of SpeciesModelManager. The necessary data sets are prepared
            according to the options in Options.hpp. */
        SpeciesModelManager();

        /** Create a new HModel, which is given a pointer to the H data stored here. A pointer
            to the radiation field is requested, to allow the H model to precalculate some
            quantities. */
        std::unique_ptr<HModel> makeHModel(const Spectrum* meanIntensityv) const;

        /** Create a new H2Model polymorphically. Depending on the heuristics given, either a
            BigH2Model or a SmallH2Model might be created. A BigH2Model will get access to the
            H2 data stored here, while a SmallH2Model ideally does not need it. A pointer to the
            radiation field is requested, to allow the H2 model to precalculate some
            quantities. */
        std::unique_ptr<H2Model> makeH2Model(const Spectrum* meanIntensityv) const;

    private:
        std::unique_ptr<HData> _hData;
        std::unique_ptr<H2Data> _h2Data;
        LookupTable _h2LTECool{"dat/h2/lte_cooling.dat", 2, 380};
    };
}
#endif  // CORE_SPECIESMODELMANAGER_HPP
