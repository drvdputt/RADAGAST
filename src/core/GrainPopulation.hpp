#ifndef CORE_GRAINPOPULATION_HPP
#define CORE_GRAINPOPULATION_HPP

#include "Array.hpp"
#include "GrainInterface.hpp"
#include "LookupTable.hpp"
#include <memory>
#include <vector>

namespace GasModule
{
    class GrainH2Formation;
    class GrainPhotoelectricData;

    /** This class contains the constants for a grain population. Given that the grain
        properties for a number of cells are the same, it can be used by multiple threads
        simultaneously, with the exception of the setDensityv() function. The latter should only
        be used in serial code. */
    class GrainPopulation
    {
    public:
        /** Default constructor and destructor, implemented as default in the cpp file (deals with
            forward declaration to unique_ptr types). */
        GrainPopulation();
        ~GrainPopulation();

        /** The grain type can be one of the types listed in the enum above. Then, an array of
            sizes, number densities, and one of temperatures needs to be specified (consider even
            using grain temperature distributions (for each grain size!) in the future). The
            temperatures are important for the H2 formation rate on the surfaces of the grains. The
            absorption efficiency also needs to be given. It needs to be given for each grain size,
            and for each point of the frequency grid of the input radiation field. It needs to
            extend to low enough frequencies, if reasonable grain temperatures are wanted.
            Otherwise, the grey body emission of the grains will be wrong.*/
        GrainPopulation(GasModule::GrainTypeLabel type, const Array& sizev, const Array& densityv,
                        const Array& temperaturev, const Array& frequencyv, const std::vector<Array>& qAbsvv);

        /** Undelete the move constructor (needed to be able to put these objects into a vector
            std::vector). Due to the unique_ptr members, it's better to have '= default' in the cpp
            file (so we can just forward declare GrainH2Formation and GrainPhotoelectricData
            here). */
        GrainPopulation(GrainPopulation&&);

        void setDensityv(const Array& densityv);

        /** @name Trivial getters. */
        /**@{*/
        const std::vector<Array>& qAbsvv() const { return _qAbsvv; }
        const std::valarray<double>& initialTemperaturev() const { return _temperaturev; }
        const std::valarray<double>& densityv() const { return _densityv; }
        const std::valarray<double>& sizev() const { return _sizev; }
        const Array& qAbsv(int m) const { return _qAbsvv[m]; }
        double density(int m) const { return _densityv[m]; }
        double size(int m) const { return _sizev[m]; }
        size_t numSizes() const;
        /**@}*/

        /** If not a nullptr, then the object can be used to calculate the h2 formation rate. */
        const GrainH2Formation* h2formationData() const { return _h2formation.get(); }

        /** If not a nullptr, then the object can be used to calculate the photoelectric effect
            (charge distribution and heating rate). */
        const GrainPhotoelectricData* photoelectricData() const { return _photoelectricData.get(); }

        /** The total thermal emission of a single grain of size _sizev[m], for a given
            temperature. Uses a lookup table, which is created at construction. [erg s-1] */
        double totalThermalEmission(int m, double T) const;

        /** Some checks for inconsistencies. Program will abort if this fails. */
        void test() const;

    private:
        Array _sizev;
        Array _densityv;
        Array _temperaturev;
        std::vector<Array> _qAbsvv;
        std::unique_ptr<GrainH2Formation> _h2formation;
        std::unique_ptr<GrainPhotoelectricData> _photoelectricData;
        // Total grey body emission for a set of temperatures (first argument of
        // LookupTable::evaluate is the size index)
        LookupTable _greybodyCoolingCurves;
    };
}  // namespace GasModule

#endif  // CORE_GRAINPOPULATION_HPP
