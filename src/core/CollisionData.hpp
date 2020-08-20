#ifndef CORE_COLLISIONDATA_HPP
#define CORE_COLLISIONDATA_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "Error.hpp"
#include <array>
#include <map>
#include <vector>

namespace RADAGAST
{
    /** Class which contains collision coefficient data. This type of data is typically tabulated
        per transition (a transition is just an abstract pair of indices here), and has a list of
        temperatures It stores the data in a contiguous way (Using an Eigen::Matrix), and maps
        pairs of level indices to entries (columns in this matrix). It also contains a list of
        temperatures, one for every row of the matrix. I chose to store the data in a
        'column-per-transition' format because Eigen works using a column based format by default.
        That way, the coefficient at different temperatures for the same transition is stored
        contiguously. */
    class CollisionData
    {
    public:
        // CONSTRUCTION AND SETUP

        /** Construct a CollisionData object. Not ready for use yet. Call the 'prepare' function,
            and then insert data for the number of transitions and temperatures you specified. */
        CollisionData();

        /** will use the given temperature grid. The second argument should equal the total number
            of transitions that will be inserted into this object. */
        void prepare(const Array& temperaturev, size_t numTransitions);

        /** Use this function to put all the data into the map. Arguments: initial level index,
            final level index, Array containing q(T) for this transition. Will throw an error if
            the number of transitions inserted exceeds the size that was reserved at
            construction. */
        void insertDataForTransition(const Array& qForEachTv, int i, int f);

        /** Performs a couple of checks on the data. To be called when the user is done putting in
            the data using the insert function. Will throw an error if there is still a column of
            zeros, or if there is a negative value in the table. Will also throw an error if
            something's off with the map or the transition list. */
        void check() const;

        /** Get a list of all the transitions (initial, final) that were inserted. */
        const std::vector<std::array<int, 2>>& transitionv() const { return _transitionv; }

        /** Interpolate all transitions at the same temperature. Returns an Array containing the
            collision coefficients for each transition in transitionv(). */
        Array qv(double T) const;

    private:
        /** Stores the temperature associated with each temperature index. [K]. */
        Array _temperaturev;

        /** The main data. Basically a copy of what was in the input file. First index =
            temperature index, second index = transition index. [cm3 s-1]. */
        EMatrix _qvv;

        /** A map from {initial level index, final level index} to transition index. */
        std::map<std::array<int, 2>, int> _transitionToColm;

        /** Also store a list of all the transitions. This makes the whole thing more
            bidirectional.*/
        std::vector<std::array<int, 2>> _transitionv;
    };
}
#endif  // CORE_COLLISIONDATA_HPP
