#include "GrainPopulation.hpp"
#include "Error.hpp"
#include "GrainH2Formation.hpp"
#include "GrainPhotoelectricData.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "Testing.hpp"

namespace GasModule
{
    GrainPopulation::GrainPopulation() = default;

    GrainPopulation::~GrainPopulation() = default;

    GrainPopulation::GrainPopulation(GasModule::GrainTypeLabel type, const Array& sizev, const Array& densityv,
                                     const Array& temperaturev, const Array& frequencyv,
                                     const std::vector<Array>& qAbsvv)
        : _sizev{sizev}, _densityv{densityv}, _temperaturev{temperaturev}, _qAbsvv{qAbsvv}
    {
        Error::equalCheck("sizev.size() and densityv.size()", _sizev.size(), _densityv.size());
        Error::equalCheck("sizev.size() and temperaturev.size()", _sizev.size(), _temperaturev.size());
        Error::equalCheck("sizev.size() and qAbsvv.size()", _sizev.size(), _qAbsvv.size());

        if (type == GasModule::GrainTypeLabel::CAR)
        {
            _h2formation = std::make_unique<GrainH2Formation>(GrainH2FormationData::carSurface,
                                                              GrainH2FormationData::grainHeatingPerH2Formed_car);
            _photoelectricData = std::make_unique<GrainPhotoelectricData>(true);
        }
        else if (type == GasModule::GrainTypeLabel::SIL)
        {
            _h2formation = std::make_unique<GrainH2Formation>(GrainH2FormationData::silSurface,
                                                              GrainH2FormationData::grainHeatingPerH2Formed_sil);
            _photoelectricData = std::make_unique<GrainPhotoelectricData>(false);
        }
        // else, no data is available and these will be nullptr

        // Optimization: pre-calculate the modified greybody cooling curve here
        Array Tv = Testing::generateGeometricGridv(100, Options::grainsolution_minGrainTemp,
                                                   Options::grainsolution_maxGrainTemp);
        Table<2> bbEmission(_sizev.size(), Tv.size());
        Array greybodyIntegrandv(frequencyv.size());

        for (size_t m = 0; m < _sizev.size(); m++)
        {
            double crossSection = Constant::PI * _sizev[m] * _sizev[m];
            for (size_t t = 0; t < Tv.size(); t++)
            {
                for (size_t i = 0; i < frequencyv.size(); i++)
                    greybodyIntegrandv[i] = qAbsvv[m][i] * SpecialFunctions::planck(frequencyv[i], Tv[t]);

                bbEmission(m, t) = Constant::FPI * crossSection
                                   * TemplatedUtils::integrate<double, Array, Array>(frequencyv, greybodyIntegrandv);
            }
        }
        _greybodyCoolingCurves = LookupTable(Tv, bbEmission);
    }

    GrainPopulation::GrainPopulation(GrainPopulation&&) = default;

    void GrainPopulation::setDensityv(const Array& densityv)
    {
        Error::equalCheck("Size of new and old list of densities", densityv.size(), _densityv.size());
        _densityv = densityv;
    }

    size_t GrainPopulation::numSizes() const { return _sizev.size(); }

    void GrainPopulation::test() const
    {
        if (TemplatedUtils::contains(0., _sizev)) Error::runtime("Grain of size 0 not allowed!");
    }

    double GrainPopulation::totalThermalEmission(int m, double T) const
    {
        return _greybodyCoolingCurves.evaluate(m, T);
    }
}  // namespace GasModule
