#include "GrainPopulation.hpp"
#include "Error.hpp"
#include "GrainH2Formation.hpp"
#include "GrainPhotoelectricData.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "Testing.hpp"

GrainPopulation::GrainPopulation() = default;

GrainPopulation::~GrainPopulation() = default;

GrainPopulation::GrainPopulation(GrainTypeLabel type, const Array& sizev, const Array& densityv,
                                 const Array& temperaturev, const Array& frequencyv, const std::vector<Array>& qAbsvv)
    : _sizev{sizev}, _densityv{densityv}, _temperaturev{temperaturev}, _qAbsvv{qAbsvv}
{
    Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
    Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(), temperaturev.size());
    Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());

    if (type == GrainTypeLabel::CAR)
    {
        _h2formation = std::make_unique<GrainH2Formation>(GrainH2FormationData::carSurface,
                                                          GrainH2FormationData::grainHeatingPerH2Formed_car);
        _photoelectricData = std::make_unique<GrainPhotoelectricData>(true);
    }
    else if (type == GrainTypeLabel::SIL)
    {
        _h2formation = std::make_unique<GrainH2Formation>(GrainH2FormationData::silSurface,
                                                          GrainH2FormationData::grainHeatingPerH2Formed_sil);
        _photoelectricData = std::make_unique<GrainPhotoelectricData>(false);
    }
    // else, no data is available and these will be nullptr

    // Optimization: pre-calculate the modified blackbody cooling curve here
    Array Tv =
        Testing::generateGeometricGridv(100, Options::grainsolution_minGrainTemp, Options::grainsolution_maxGrainTemp);
    Table<2> bbEmission(numSizes(), Tv.size());
    for (size_t m = 0; m < numSizes(); m++)
    {
        double crossSection = Constant::PI * _sizev[m] * _sizev[m];
        for (int t = 0; t < Tv.size(); t++)
        {
            Array blackbodyIntegrandv(frequencyv.size());
            for (size_t i = 0; i < frequencyv.size(); i++)
                blackbodyIntegrandv[i] = qAbsvv[m][i] * SpecialFunctions::planck(frequencyv[i], Tv[t]);

            bbEmission(m, t) =
                crossSection * TemplatedUtils::integrate<double, Array, Array>(frequencyv, blackbodyIntegrandv);
        }
    }
    _blackbodyCoolingCurves = LookupTable(Tv, bbEmission);
}

GrainPopulation::GrainPopulation(GrainPopulation&&) = default;

void GrainPopulation::test() const
{
    if (TemplatedUtils::contains(0., _sizev)) Error::runtime("Grain of size 0 not allowed!");
}

double GrainPopulation::totalThermalEmission(int m, double T) const
{
    return _blackbodyCoolingCurves.evaluate(m, T);
}
