#ifndef CORE_GRAININTERFACE_HPP
#define CORE_GRAININTERFACE_HPP

#include <functional>
#include <memory>
#include <valarray>
#include <vector>

class GrainPopulation;

namespace GasModule
{
    /** Class that a client code should use to pass the grain properties. This class is a wrapper
        around a list of grain models (GrainPopulation objects). */
    class GrainInterface
    {
    public:
        /** Constructor which takes a vector of predefined populations. Please pass the vector
            using a unique pointer and std::move. Ownership over the vector of populations will
            then be transferred to this class. */
        GrainInterface(std::unique_ptr<std::vector<GrainPopulation>> populationvToMove);

        /** Default constructor, equivalent to no grains at all. */
        GrainInterface();
        ~GrainInterface();

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
