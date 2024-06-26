#ifndef CORE_GASDIAGNOSTICS_HPP
#define CORE_GASDIAGNOSTICS_HPP

#include <map>
#include <string>
#include <valarray>
#include <vector>

namespace RADAGAST
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

        /** All heating rates (including total photoelectric), stored as a map with keys */
        const std::map<std::string, double>& heating() const { return _heatingm; }

        /** All cooling rates, stored as a map with keys */
        const std::map<std::string, double>& cooling() const { return _coolingm; }

        /** Names and rates of the reactions in the chemical network */
        const std::map<std::string, double> reactionInfov() const { return _reactionInfov; }

        /** Population density of the included H levels */
        const std::valarray<double>& hPopulations() const { return _hPopulationv; }

        /** Population density of the included H2 levels */
        const std::valarray<double>& h2Populations() const { return _h2Populationv; }

        /** Things that fall outside of a specific category */
        const std::map<std::string, double>& userValues() const { return _userValuem; }

        void setReactionInfo(const std::vector<std::string>& reactionNames, const std::vector<double>& reactionRates);
        void setHPopulations(const std::valarray<double>& hp) { _hPopulationv = hp; }
        void setH2Populations(const std::valarray<double>& h2p) { _h2Populationv = h2p; }

        /** Add heating entry with the label "heat:name" */
        void setHeating(const std::string& name, double value);
        /** Add cooling entry with the label "cool:name" */
        void setCooling(const std::string& name, double value);
        void setUserValue(const std::string& key, double value);

    private:
        bool _saveLevelPopulations{true};

        std::map<std::string, double> _heatingm;
        std::map<std::string, double> _coolingm;
        std::map<std::string, double> _reactionInfov;

        std::valarray<double> _hPopulationv;
        std::valarray<double> _h2Populationv;

        std::map<std::string, double> _userValuem;
    };
}
#endif  // CORE_GASDIAGNOSTICS_HPP
