#ifndef CORE_GRAININTERFACE_HPP
#define CORE_GRAININTERFACE_HPP

#include <functional>
#include <memory>
#include <valarray>
#include <vector>

namespace GasModule
{
    class GrainPopulation;

    /** List of grain types which have built-in values for H2 formation and the photoelectric
        effect. For anything else, the label 'OTHER' can be used, but then no contribution to the
        H2 formation or photoelectric heating will be made by this grain population. */
    enum class GrainTypeLabel { CAR, SIL, OTHER };

    /** Class that a client code should use to pass the grain properties. This class is a wrapper
        around a list of grain models (GrainPopulation objects). */
    class GrainInterface
    {
    public:
        /** Default constructor, equivalent to no grains at all. */
        GrainInterface();
        ~GrainInterface();

        /** Clients can use this to add populations to this GrainInterface. It is a wrapper around
            the constructor of GrainPopulation. This way we avoid having to pull in the
            dependencies of GrainPopulation. The frequency list should be the same as iFrequencyv
            of the GasInterface (TODO: make this more elegant), and indicates for which frequencies
            the given qAbs is tabulated. */
        void addPopulation(GasModule::GrainTypeLabel type, const std::valarray<double>& sizev,
                           const std::valarray<double>& densityv, const std::valarray<double>& temperaturev,
                           const std::valarray<double>& frequencyv, const std::vector<std::valarray<double>>& qAbsvv);

        /** Change the number densities of the grains of an existing population (at index p) in
            this GrainInterface. */
        void changePopulationDensityv(int p, const std::valarray<double>& densityv);

        /** Undelete the move constructor, return this object from functions etc. This is necessary
            because the presence of a unique_ptr. Very annoying if you don't know about this
            behaviour, especially since the compiler errors can be quite cryptic. */
        GrainInterface(GrainInterface&&);

        /** The number of grain populations. */
        size_t numPopulations() const;

        /** Get a pointer to the whole population vector. Return nullptr if no grains are
            present. */
        const std::vector<GrainPopulation>* populationv() const;

        /** A quick test to see if all population objects in the vector have reasonable
            contents. */
        void test() const;

    private:
        std::unique_ptr<std::vector<GrainPopulation>> _populationv;
    };
} /* namespace GasModule */

#endif  // CORE_GRAININTERFACE_HPP
