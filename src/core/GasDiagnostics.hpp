#ifndef CORE_GASDIAGNOSTICS_HPP
#define CORE_GASDIAGNOSTICS_HPP

#include <map>
#include <valarray>
#include <vector>

namespace GasModule
{
    /** This class should be used to store any extra information that is not put into the standard
        GasState, to avoid making it too bloaty. Maybe it will make sense to design it so that not
        all the members have to be set. The amount of detail can maybe be set using some sort of
        configuration file. */
    class GasDiagnostics
    {
    public:
        /** Configuration */
        bool saveLevelPopulations() const { return _saveLevelPopulations; }

        /** Invidivual photoelectric heating contributions. TODO: make this work per population
            and/or per size. */
        const std::valarray<double>& photoelectricHeating() const { return _photoelectricHeating; }

        /** All heating rates (including total photoelectric), stored as a map with keys */
        const std::map<std::string, double>& heating() const { return _heatingm; }

        /** All cooling rates, stored as a map with keys */
        const std::map<std::string, double>& cooling() const { return _coolingm; }

        /** Names of the reactions in the chemical network used */
        const std::vector<std::string>& reactionNames() const { return _reactionNamev; }

        /** Rates of the reactions in the chemical network used */
        const std::valarray<double>& reactionRates() const { return _reactionRatev; }

        /** Population density of the included H levels */
        const std::valarray<double>& hPopulations() const { return _hPopulationv; }

        /** Population density of the included H2 levels */
        const std::valarray<double>& h2Populations() const { return _h2Populationv; }

        /** Things that fall outside of a specific category */
        const std::map<std::string, double>& userValues() const { return _userValuem; }

        void setPhotoelectricHeating(const std::valarray<double>& peh) { _photoelectricHeating = peh; }
        void setReactionNames(const std::vector<std::string>& rn) { _reactionNamev = rn; }
        void setReactionRates(const std::valarray<double>& rr) { _reactionRatev = rr; }
        void setHPopulations(const std::valarray<double>& hp) { _hPopulationv = hp; }
        void setH2Populations(const std::valarray<double>& h2p) { _h2Populationv = h2p; }
        void setHeating(const std::string& key, double value);
        void setCooling(const std::string& key, double value);
        void setUserValue(const std::string& key, double value);

    private:
        bool _saveLevelPopulations{true};

        std::valarray<double> _photoelectricHeating;

        std::map<std::string, double> _heatingm;
        std::map<std::string, double> _coolingm;

        std::vector<std::string> _reactionNamev;
        std::valarray<double> _reactionRatev;

        std::valarray<double> _hPopulationv;
        std::valarray<double> _h2Populationv;

        std::map<std::string, double> _userValuem;
    };
}
#endif  // CORE_GASDIAGNOSTICS_HPP
