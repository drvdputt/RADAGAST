#include "GasSolution.hpp"
#include "FreeBound.hpp"
#include "FreeFree.hpp"
#include "GasDiagnostics.hpp"
#include "GasStruct.hpp"
#include "GrainPhotoelectricEffect.hpp"
#include "Ionization.hpp"
#include "Options.hpp"

GasSolution::GasSolution(const GasModule::GrainInterface& gri,
                         const Spectrum& specificIntensity, const SpeciesIndex& speciesIndex,
                         std::unique_ptr<HModel> hModel, std::unique_ptr<H2Model> h2Model,
                         const FreeBound& freeBound, const FreeFree& freeFree)
                : _specificIntensity{specificIntensity}, _sv(speciesIndex),
                  _hSolution(std::move(hModel)),
                  _h2Solution(std::move(h2Model)), _freeBound{freeBound}, _freeFree{freeFree}
{
	_grainSolutionv.reserve(gri.numPopulations());
	for (const auto& p : *gri.populationv())
		_grainSolutionv.emplace_back(p);
}

void GasSolution::makeZero()
{
	_sv.setDensities(EVector::Zero(_sv.size()));
	solveLevels();
}

void GasSolution::solveLevels(double formH2)
{
	GasStruct gas = {_t, _sv};
	_hSolution->solve(nH(), gas, _specificIntensity);
	_h2Solution->solve(nH2(), gas, _specificIntensity, formH2);
}

void GasSolution::solveGrains()
{
	GrainPhotoelectricCalculator::Environment env(
	                _specificIntensity, _t, ne(), np(), {-1, 1, 0, 0},
	                {ne(), np(), nH(), nH2()},
	                {Constant::ELECTRONMASS, Constant::PROTONMASS, Constant::HMASS_CGS,
	                 2 * Constant::HMASS_CGS});

	for (auto& g : _grainSolutionv)
	{
		// I think charge distributions are best recalculated before the new
		// temperatures (they affect the collisional processes)
		g.recalculateChargeDistributions(env);
		g.recalculateTemperatures(env);
	}
}

Array GasSolution::emissivityv(const Array& eFrequencyv) const
{
	Array lineEmv(eFrequencyv.size());
	lineEmv += _hSolution->emissivityv(eFrequencyv);
	lineEmv += _h2Solution->emissivityv(eFrequencyv);

	Array contEmCoeffv(eFrequencyv.size());
	_freeBound.addEmissionCoefficientv(_t, eFrequencyv, contEmCoeffv);
	_freeFree.addEmissionCoefficientv(_t, eFrequencyv, contEmCoeffv);

	return lineEmv + (np() * ne() / Constant::FPI) * contEmCoeffv;
}

Array GasSolution::opacityv(const Array& oFrequencyv) const
{
	size_t numFreq = oFrequencyv.size();
	Array lineOpv(numFreq);
	lineOpv += _hSolution->opacityv(oFrequencyv);
	lineOpv += _h2Solution->opacityv(oFrequencyv);

	Array contOpv(numFreq);
	_freeFree.addOpacityCoefficientv(_t, oFrequencyv, contOpv);
	contOpv *= np() * ne();

	Array totalOpv(numFreq);
	for (size_t i = 0; i < numFreq; i++)
	{
		// TODO: this should actually be the average over the cross section for
		// this frequency bin
		double ionizOp_iFreq = nH() * Ionization::crossSection(oFrequencyv[i]);
		totalOpv[i] = ionizOp_iFreq + contOpv[i] + lineOpv[i];
	}
	return totalOpv;
}

double GasSolution::cooling() const
{
	double freefreeCool = _freeFree.cooling(np() * ne(), _t);
	double hRecCool = Ionization::cooling(nH(), np(), ne(), _t);
	return freefreeCool + hRecCool;
}

double GasSolution::heating() const
{
	// double freefreeHeat = _freeFree->heating(np(s) * ne(s), s.T, s.specificIntensity);
	// TODO: decide whether to keep the above, as it is negligible in any case I can imagine
	double lineHeat = _hSolution->netHeating() + _h2Solution->netHeating();
	double hPhotoIonHeat = Ionization::heating(np(), ne(), _t, _specificIntensity);
	double dissHeat = _h2Solution->dissociationHeating(_specificIntensity);
	double grainHeat = grainHeating();
	return lineHeat + hPhotoIonHeat + dissHeat + grainHeat;
}

double GasSolution::grainHeating(double* photoHeat, double* collCool) const
{
	double grainPhotoelectricHeating = 0;
	double gasGrainCooling = 0;

	// Specify the environment parameters
	GrainPhotoelectricCalculator::Environment env(
	                _specificIntensity, _t, ne(), np(), {-1, 1, 0, 0},
	                {ne(), np(), nH(), nH2()},
	                {Constant::ELECTRONMASS, Constant::PROTONMASS, Constant::HMASS_CGS,
	                 2 * Constant::HMASS_CGS});

	for (const auto& g : _grainSolutionv)
	{
		// This assumes that solveGrains was already called
		grainPhotoelectricHeating += g.photoelectricGasHeating(env);
		if (Options::cooling_gasGrainCollisions)
			gasGrainCooling += g.collisionalGasCooling(env);
	}

	if (photoHeat != nullptr)
		*photoHeat = grainPhotoelectricHeating;
	if (collCool != nullptr)
		*collCool = gasGrainCooling;

	return grainPhotoelectricHeating - gasGrainCooling;
}

void GasSolution::fillDiagnostics(GasDiagnostics* gd) const

{
	if (!gd)
		Error::runtime("GasDiagnostics is nullptr!");

	double h2form = kGrainH2FormationRateCoeff();
	double h2dissoc = _h2Solution->dissociationRate(_specificIntensity);

	double hphotoion = Ionization::photoRateCoeff(_specificIntensity);
	double hcolion = Ionization::collisionalRateCoeff(_t);
	double hrec = Ionization::recombinationRateCoeff(_t);

	gd->setReactionNames({"h2form", "h2dissoc", "hphotoion", "hcolion", "hrec"});
	gd->setReactionRates({h2form, h2dissoc, hphotoion, hcolion, hrec});

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

	gd->setHeating("H ion", Ionization::heating(np(), ne(), _t, _specificIntensity));
	gd->setCooling("Hrec", Ionization::cooling(nH(), np(), ne(), _t));
	gd->setHeating("H deexc", netHline);
	gd->setCooling("H exc", -netHline);

	gd->setHeating("H2 deexc", netH2line);
	gd->setCooling("H2 exc", -netH2line);
	gd->setHeating("H2 dissoc", _h2Solution->dissociationHeating(_specificIntensity));
	// gd->setHeating("freefree", x);
	gd->setCooling("freefree", _freeFree.cooling(np() * ne(), _t));

	// I need this per grain size. Doing this thing for now.
	double grainPhotoHeat = 0;
	double grainCollCool = 0;
	double totalGrainHeat = grainHeating(&grainPhotoHeat, &grainCollCool);
	gd->setPhotoelectricHeating(Array({totalGrainHeat}));
	gd->setHeating("total grainphoto", grainPhotoHeat);
	gd->setCooling("grain collisions", grainCollCool);
}

GasModule::GasState GasSolution::makeGasState(const Array& oFrequencyv,
                                              const Array& eFrequencyv) const
{
	Array emv, opv;
	if (eFrequencyv.size() > 2)
		emv = emissivityv(eFrequencyv);
	if (oFrequencyv.size() > 2)
	{
		opv = opacityv(oFrequencyv);
	}
	Array densityv(_sv.data(), _sv.size());
	return {emv, opv, _t, densityv};
}

double GasSolution::kDissH2Levels() const
{
	return _h2Solution->dissociationRate(_specificIntensity);
}

double GasSolution::kGrainH2FormationRateCoeff() const
{
	double total = 0;
	for (const auto& g : _grainSolutionv)
		total += g.surfaceH2FormationRateCoeff(_t);
	return total;
}
