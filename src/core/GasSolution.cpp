#include "GasSolution.hpp"
#include "CollisionParameters.hpp"
#include "DebugMacros.hpp"
#include "FreeBound.hpp"
#include "FreeFree.hpp"
#include "GasDiagnostics.hpp"
#include "GrainPopulation.hpp"
#include "Ionization.hpp"
#include "Options.hpp"

namespace RADAGAST
{
    GasSolution::GasSolution(const RADAGAST::GrainInterface* gri, const Spectrum* meanIntensity, double fshield,
                             const SpeciesIndex* speciesIndex, std::unique_ptr<HModel> hModel,
                             std::unique_ptr<H2Model> h2Model, const FreeBound* freeBound, const FreeFree* freeFree)
        : _sv(speciesIndex), _hSolution(std::move(hModel)),
          _h2Solution(std::move(h2Model)), _meanIntensity{meanIntensity}, _fshield{fshield},
          _freeBound{freeBound}, _freeFree{freeFree}, _ionHeatPerH{Ionization::heatingPerH(*meanIntensity)}
    {
        int numPop = gri->numPopulations();
        if (numPop)
        {
            _grainSolutionv.reserve(numPop);
            for (int i = 0; i < numPop; i++) _grainSolutionv.emplace_back(&gri->populationv()->at(i), meanIntensity);
        }
    }

    void GasSolution::makeZero()
    {
        _sv.setDensities(EVector::Zero(_sv.size()));
        solveLevels();
    }

    void GasSolution::solveLevels()
    {
        CollisionParameters cp = {_t, _sv, _h2Solution->orthoPara()};
        _hSolution->solve(nH(), cp);
        _h2Solution->solve(nH2(), cp, _kGrainH2FormationRateCoeff * _sv.nH());
    }

    void GasSolution::solveGrains()
    {
        for (auto& g : _grainSolutionv) g.recalculate(_t, _sv);

        _kGrainH2FormationRateCoeff = 0.;
        for (const auto& g : _grainSolutionv) _kGrainH2FormationRateCoeff += g.surfaceH2FormationRateCoeff(_t);
    }

    Array GasSolution::emissivityv(const Array& eFrequencyv) const
    {
        Array lineEmv(eFrequencyv.size());
        lineEmv += _hSolution->emissivityv(eFrequencyv);
        lineEmv += _h2Solution->emissivityv(eFrequencyv);

        Array contEmCoeffv(eFrequencyv.size());
        _freeBound->addEmissionCoefficientv(_t, eFrequencyv, contEmCoeffv);
        _freeFree->addEmissionCoefficientv(_t, eFrequencyv, contEmCoeffv);

        return lineEmv + (np() * ne() / Constant::FPI) * contEmCoeffv;
    }

    Array GasSolution::opacityv(const Array& oFrequencyv) const
    {
        size_t numFreq = oFrequencyv.size();
        Array lineOpv(numFreq);
        lineOpv += _hSolution->opacityv(oFrequencyv);
        lineOpv += _h2Solution->opacityv(oFrequencyv);

        Array contOpv(numFreq);
        _freeFree->addOpacityCoefficientv(_t, oFrequencyv, contOpv);
        contOpv *= np() * ne();

        Array totalOpv(numFreq);
        for (size_t i = 0; i < numFreq; i++)
        {
            // evaluate ionization cross section at given frequencies
            double ionizOp_iFreq = nH() * Ionization::crossSection(oFrequencyv[i]);
            totalOpv[i] = ionizOp_iFreq + contOpv[i] + lineOpv[i];
        }
        return totalOpv;
    }

    double GasSolution::cooling() const
    {
        double freefreeCool = _freeFree->cooling(np() * ne(), _t);
        double hRecCool = Ionization::cooling(nH(), np(), ne(), _t);
        double grainCool = grainCooling();
        DEBUG("Cooling contributions: FF" << freefreeCool << " FB " << hRecCool << " Grcol " << grainCool << '\n');
        return freefreeCool + hRecCool + grainCool;
    }

    double GasSolution::heating()
    {
        double hLine = _hSolution->netHeating();
        double hPhotoIonHeat = nH() * _ionHeatPerH;

        // scale these H2 quantities with the shielding factor for now. A more physically
        // motivated approach might be implemented later at a lower level, inside the H2 model
        double h2Line = _h2Solution->netHeating();
        double dissHeat = _h2Solution->dissociationHeating();

        double grainHeat = grainHeating();
        DEBUG("Heating contributions: Hln " << hLine << " H2ln " << h2Line << " Hphot " << hPhotoIonHeat << " H2diss "
                                            << dissHeat << " Grphot " << grainHeat << '\n');
        return hLine + h2Line + hPhotoIonHeat + dissHeat + grainHeat;
    }

    double GasSolution::grainHeating()
    {
        // This assumes that solveGrains() was already called
        double grainPhotoelectricHeating = 0.;
        for (auto& g : _grainSolutionv) grainPhotoelectricHeating += g.photoelectricGasHeating();
        return grainPhotoelectricHeating;
    }

    double GasSolution::grainCooling() const
    {
        if (!Options::cooling_gasGrainCollisions) return 0.;

        // This assumes that solveGrains was already called
        double gasGrainCooling = 0.;
        for (const auto& g : _grainSolutionv) gasGrainCooling += g.collisionalGasCooling();
        return gasGrainCooling;
    }

    void GasSolution::fillDiagnostics(GasDiagnostics* gd)

    {
        if (!gd) Error::runtime("GasDiagnostics is nullptr!");

        double h2form = kGrainH2FormationRateCoeff();
        double h2dissoc = _h2Solution->dissociationRate();

        double hphotoion = Ionization::photoRateCoeff(*_meanIntensity);
        double hcolion = Ionization::collisionalRateCoeff(_t);
        double hrec = Ionization::recombinationRateCoeff(_t);

        gd->setReactionInfo({"H2 form", "H2 diss", "H photoion", "H colion", "H rec"},
                            {h2form, h2dissoc, hphotoion, hcolion, hrec});

        if (gd->saveLevelPopulations())
        {
            EVector hlv = _hSolution->levelSolution()->nv();
            gd->setHPopulations(Array(hlv.data(), hlv.size()));
            if (_h2Solution->hasLevels())
            {
                EVector h2lv = _h2Solution->levelSolution()->nv();
                gd->setH2Populations(Array(h2lv.data(), h2lv.size()));
            }
        }
        double netHline = _hSolution->netHeating();
        double netH2line = _h2Solution->netHeating();

        gd->setHeating("H photoion", nH() * _ionHeatPerH);
        gd->setCooling("H rec", Ionization::cooling(nH(), np(), ne(), _t));
        gd->setHeating("H dxc", netHline);
        gd->setCooling("H exc", -netHline);

        gd->setHeating("H2 dxc", netH2line);
        gd->setCooling("H2 exc", -netH2line);
        gd->setHeating("H2 diss", _h2Solution->dissociationHeating());
        gd->setCooling("freefree", _freeFree->cooling(np() * ne(), _t));

        // I need this per grain size. Doing this thing for now.
        double grainPhotoHeat = grainHeating();
        double grainCollCool = grainCooling();
        gd->setHeating("grain photo", grainPhotoHeat);
        gd->setCooling("grain col", grainCollCool);

        // Other things will be written to the 'user values' of gasDiagnostics, somewhere in
        // these calls
        _h2Solution->extraDiagnostics(*gd);

        // write grain temperatures
        for (size_t i = 0; i < _grainSolutionv.size(); ++i)
        {
            std::string label = "pop" + std::to_string(i);
            _grainSolutionv[i].extraDiagnostics(*gd, label);
        }
    }

    void GasSolution::setGasState(RADAGAST::GasState& g) const
    {
        g.setMembers(_t, {_sv.data(), _sv.size()}, _hSolution->n2s());
    }

    double GasSolution::kDissH2Levels() const { return _h2Solution->dissociationRate(); }

    double GasSolution::kGrainH2FormationRateCoeff() const { return _kGrainH2FormationRateCoeff; }
}
