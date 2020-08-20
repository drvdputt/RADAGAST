#include "CollisionData.hpp"
#include "Error.hpp"
#include "TemplatedUtils.hpp"

namespace RADAGAST
{
    CollisionData::CollisionData() = default;

    void CollisionData::prepare(const Array& temperaturev, size_t numTransitions)
    {
        _temperaturev = temperaturev;
        _qvv = EMatrix::Zero(temperaturev.size(), numTransitions);

        // This is usually not necessary, since in recommended use, prepare is only called on a
        // freshly constructed instance. But just in case.
        _transitionToColm.clear();
        _transitionv.clear();
    }

    void CollisionData::insertDataForTransition(const Array& qForEachTv, int i, int f)

    {
        // Index for new transition == previous + 1 == current size of map
        int transitionIndex = _transitionToColm.size();
        if (transitionIndex >= _qvv.cols())
            Error::runtime("CollisionData object is full! Can't insert more transitions.");

        // The number of rows is already set. An error should appear when qForEachTv is of the wrong size.
        // Creating the map does not transfer ownership. The assigment is a copy operation.
        _qvv.col(transitionIndex) = Eigen::Map<const EVector>(&qForEachTv[0], qForEachTv.size());
        _transitionToColm.insert({{i, f}, transitionIndex});
        _transitionv.emplace_back(std::array<int, 2>{i, f});
    }

    void CollisionData::check() const
    {
        Error::equalCheck<size_t>("Transition map size and data matrix number of columns", _transitionToColm.size(),
                                  _qvv.cols());
        Error::equalCheck("Transition map and transition list size", _transitionToColm.size(), _transitionv.size());

        for (int c = 0; c < _qvv.cols(); c++)
        {
            if ((_qvv.col(c).array() == 0).all())
                Error::runtime("Zero column in q-data matrix. Some data was probably not "
                               "filled in");
            if ((_qvv.col(c).array() < 0).any()) Error::runtime("Negative value in collision data!");
        }
    }

    Array CollisionData::qv(double T) const
    {
        // indices between which to interpolate
        int left, right;
        std::tie(left, right) = TemplatedUtils::saneIndexPair(T, _temperaturev);

        int numTransitions = _transitionv.size();
        Array result(numTransitions);
        for (int transitionIndex = 0; transitionIndex < numTransitions; transitionIndex++)
        {
            // extrapolates linearly based on last two points, so need to clip at 0
            result[transitionIndex] = std::max(
                0., TemplatedUtils::interpolateLinear(T, _temperaturev[left], _temperaturev[right],
                                                      _qvv(left, transitionIndex), _qvv(right, transitionIndex)));
        }
        return result;
    }
}
