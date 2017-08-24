#include "GasInterfaceImpl.h"
#include "Constants.h"

#include "ChemicalNetwork.h"
#include "ChemistrySolver.h"
#include "DebugMacros.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "GrainProperties.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "Testing.h"

using namespace std;

GasInterfaceImpl::GasInterfaceImpl(unique_ptr<NLevel> atomModel, bool molecular,
                                   const Array& frequencyv)
                : _frequencyv(frequencyv), _atomicLevels(std::move(atomModel)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
	if (molecular)
	{
		_molecularLevels = make_unique<H2Levels>(make_shared<H2FromFiles>(), frequencyv);
		_chemSolver = make_unique<ChemistrySolver>(move(make_unique<ChemicalNetwork>()));
	}

	ine = ChemicalNetwork::speciesIndexm.at("e-");
	inp = ChemicalNetwork::speciesIndexm.at("H+");
	inH = ChemicalNetwork::speciesIndexm.at("H");
	inH2 = ChemicalNetwork::speciesIndexm.at("H2");
}

GasInterfaceImpl::~GasInterfaceImpl() = default;

void GasInterfaceImpl::solveInitialGuess(GasModule::GasState& gs, double n, double T,
                                         const GasModule::GrainInfo& gi) const
{
	Array isrfGuess(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];
		isrfGuess[iFreq] = SpecialFunctions::planck(freq, T);
	}
	solveBalance(gs, n, T, isrfGuess, gi);
}

void GasInterfaceImpl::solveBalance(GasModule::GasState& gs, double n, double Tinit,
                                    const Array& specificIntensityv,
                                    const GasModule::GrainInfo& gi) const
{
#ifndef SILENT
	double isrf = TemplatedUtils::integrate<double>(_frequencyv, specificIntensityv);

	DEBUG("Solving balance under isrf of " << isrf << " erg / s / cm2 / sr = ");
	DEBUG(isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);
#endif
	Solution s;

	if (n > 0)
	{
		const double Tmax = 1000000.;
		const double Tmin = 10;
		const double logTmax = log10(Tmax);
		const double logTmin = log10(Tmin);
		const double logTtolerance = 1.e-4;

		double logTinit = log10(Tinit);

		/* Lambda function that will be used by the search algorithm. The state of the
		   system will be updated every time the algorithm calls this function. The return
		   value indicates whether the temperature should increase (there is net heating so
		   we need a higher temperature leading to more cooling) or decrease (there is net
		   cooling so we need a lower temperature leading to less cooling). */
		int counter = 0;
		function<int(double)> evaluateThermalBalance = [&](double logT) -> int {
			counter++;
			s = calculateDensities(n, pow(10., logT), specificIntensityv);
			double netPowerIn = heating(s, gi) - cooling(s);
			DEBUG("Cycle " << counter << ": logT = " << logT
			               << "; netHeating = " << netPowerIn << endl
			               << endl);
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(
		                evaluateThermalBalance, logTinit, logTtolerance, logTmax, logTmin);

		// Evaluate the densities for one last time, using the final temperature.
		s = calculateDensities(n, pow(10., logTfinal), specificIntensityv);
	}
	else
	{
		s = calculateDensities(0, 0, specificIntensityv);
	}

	const Array& emv = emissivityv(s);
	const Array& opv = opacityv(s);
	const Array& scv = scatteringOpacityv(s);

#ifdef SANITY
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		if (specificIntensityv[iFreq] < 0 || emv[iFreq] < 0 || opv[iFreq] < 0 ||
		    scv[iFreq] < 0)
		{
			Error::runtime("GasModule: negative value in one of the optical "
			               "properties");
		}
	}
#endif
	// Put the relevant data into the gas state
	gs = GasModule::GasState(specificIntensityv, emv, opv, scv, s.T, f(s));
}

GasInterfaceImpl::Solution
GasInterfaceImpl::calculateDensities(double ntotal, double T, const Array& specificIntensityv) const
{
	Solution s;
	s.T = T;
	s.specificIntensityv = specificIntensityv;

	/* Lambda function, because it is only needed in this scope. The [&] passes the current
	   scope by reference, so the lambda function can modify s. */
	auto solveLevelBalances = [&]() {
		double nH = s.abundancev[inH];
		double ne = s.abundancev[ine];
		double np = s.abundancev[inp];
		double nH2 = s.abundancev[inH2];
		s.HSolution = _atomicLevels->solveBalance(nH, ne, np, T, specificIntensityv);
		if (_molecularLevels)
			s.H2Solution = _molecularLevels->solveBalance(nH2, ne, np, T,
			                                              specificIntensityv);
	};

	if (ntotal > 0)
	{
		DEBUG("Calculating state for T = " << T << "K" << endl);

		// Initial guess for the chemistry. Rather important for getting good convergence.
		double iniNH2 = ntotal / 4;
		double iniAtomAndIon = ntotal / 2;
		double guessF = Ionization::solveBalance(iniAtomAndIon, T, _frequencyv,
		                                         specificIntensityv);
		s.abundancev = EVector(4);
		s.abundancev << guessF * iniAtomAndIon, guessF * iniAtomAndIon,
		                (1 - guessF) * iniAtomAndIon, iniNH2;
		/* Note that the total density of H nuclei is 0 * ne + 1 * np + 1 * nH / 2 + 2 * nH2
		   / 4 = 0 + n / 2 + 2n / 2 = ntotal */

		// Solve the levels for the first time using the initial guess
		solveLevelBalances();

		int counter = 0;
		bool stopCriterion = false;
		while (!stopCriterion)
		{
			EVector previousAbundancev = s.abundancev;

			// When including H2
			if (_molecularLevels)
			{
				// Calculate fixed rate coefficients
				double kFromH2Levels =
				                _molecularLevels->dissociationRate(s.H2Solution);
				EVector reactionRates = _chemSolver->chemicalNetwork()->rateCoeffv(
				                T, _frequencyv, specificIntensityv, kFromH2Levels);

				// Solve chemistry network
				s.abundancev = _chemSolver->solveBalance(reactionRates,
				                                         s.abundancev);

				/* TODO: Add effect of grain charging to chemical network. I think
				   it might be possible to do this by imposing a conservation
				   equation for the number of electrons: ne + nH + nH2 = (ne + nH +
				   nH2)_0 + <Cg>*ng The average grain charge <Gg> should be updated
				   together with the rates I guess?  Another option would be to
				   include the grain charge rates into the network as extra
				   reactions. The production vector would be (1 0 0 0) while the
				   reactant vector would be zero (the grains don't disappear when
				   they lose an electron) Grain recombination / charge exchange
				   reaction could also be added. I need to think about whether the
				   'disappearing' particles will cause problems when couples with
				   conservation equations. */
			}
			// When ignoring H2
			else
			{
				// Just solve the ionization balance in the nebular approximation
				double f = Ionization::solveBalance(ntotal, T, _frequencyv,
				                                    specificIntensityv);
				DEBUG("Ionized fraction = " << f << endl);

				// Neutral fraction
				s.abundancev[inH] = ntotal * (1 - f);
				// Ionized fraction
				s.abundancev[inp] = ntotal * f;
				// Electron density is simply equal to proton density
				s.abundancev[ine] = s.abundancev[inp];
				s.abundancev[inH2] = 0;
			}
			solveLevelBalances();

			// An abundance has converged if it changes by less than 1%
			EArray changev = s.abundancev - previousAbundancev;
			// Or if it is negligible compared to the norm (or sum maybe?)
			double norm = s.abundancev.norm();
			Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
			                changev.abs() <= 0.01 * previousAbundancev.array() ||
			                s.abundancev.array() < 1.e-99 * norm;
			counter++;
			DEBUG("Chemistry: " << counter << endl
			                    << s.abundancev << endl
			                    << "convergence: \n"
			                    << convergedv << endl);

			// Currently, the implementation without molecules does not need iteration.
			stopCriterion = !_molecularLevels || convergedv.all();
		}
	}
	else
	{
		s.HSolution = _atomicLevels->solveBalance(0, 0, 0, T, specificIntensityv);
		if (_molecularLevels)
			s.H2Solution = _molecularLevels->solveBalance(0, 0, 0, T,
			                                              specificIntensityv);
		s.abundancev = EVector::Zero(ChemicalNetwork::speciesIndexm.size());
	}
	return s;
}

Array GasInterfaceImpl::emissivityv(const Solution& s) const
{
	Array lineEmv = _atomicLevels->emissivityv(s.HSolution);
	if (_molecularLevels)
		lineEmv += _molecularLevels->emissivityv(s.H2Solution);

	Array contEmCoeffv(_frequencyv.size());
	_freeBound->addEmissionCoefficientv(s.T, contEmCoeffv);
	_freeFree->addEmissionCoefficientv(s.T, contEmCoeffv);

	return lineEmv + (np_ne(s) / Constant::FPI) * contEmCoeffv;
}

Array GasInterfaceImpl::opacityv(const Solution& s) const
{
	Array lineOp = _atomicLevels->opacityv(s.HSolution);
	if (_molecularLevels)
		lineOp += _molecularLevels->emissivityv(s.H2Solution);

	Array contOpCoeffv(_frequencyv.size());
	_freeFree->addOpacityCoefficientv(s.T, contOpCoeffv);

	double npne = np_ne(s);
	double nH0 = nAtomic(s);
	Array totalOp(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double ionizOp_iFreq = nH0 * Ionization::crossSection(_frequencyv[iFreq]);
		totalOp[iFreq] = ionizOp_iFreq + npne * contOpCoeffv[iFreq] + lineOp[iFreq];
#ifdef SANITY
		if (totalOp[iFreq] < 0)
		{
			cout << "Negative opacity!";
		}
#endif
	}
	return totalOp;
}

Array GasInterfaceImpl::scatteringOpacityv(const Solution& s) const
{
	return Array(s.specificIntensityv.size());
}

double GasInterfaceImpl::cooling(const Solution& s) const
{
	double lineCool = lineCooling(s);
	double contCool = continuumCooling(s);
	DEBUG("cooling: line / cont = " << lineCool << " / " << contCool << endl);
	return lineCool + contCool;
}

double GasInterfaceImpl::heating(const Solution& s) const
{
	double lineHeat = lineHeating(s);
	double contHeat = continuumHeating(s);
	DEBUG("heating: line / cont = " << lineHeat << " / " << contHeat << endl);
	return lineHeat + contHeat;
}

double GasInterfaceImpl::heating(const Solution& s, const GasModule::GrainInfo& g) const
{
	double grainPhotoelectricHeating{0};

	for (const auto& pop : g.populationv())
	{
		auto type = pop._type;
		if (GrainProperties::heatingAvailable(type))
		{
			// Choose the correct parameters for the photoelectric heating calculation
			GrainPhotoelectricEffect gpe(type);

			// Specify the environment parameters
			double ne = s.abundancev[ine];
			double np = s.abundancev[inp];
			GrainPhotoelectricEffect::Environment env(
			                _frequencyv, s.specificIntensityv, s.T, ne, np, {-1, 1},
			                {ne, np}, {Constant::ELECTRONMASS, Constant::PROTONMASS});

			grainPhotoelectricHeating += gpe.heatingRate(env, pop);
		}
	}
	return heating(s) + grainPhotoelectricHeating;
}

double GasInterfaceImpl::lineCooling(const Solution& s) const
{
	double result = _atomicLevels->cooling(s.HSolution);
	if (_molecularLevels)
		result += _molecularLevels->cooling(s.H2Solution);
	return result;
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	double result = _atomicLevels->heating(s.HSolution);
	if (_molecularLevels)
		result += _molecularLevels->heating(s.H2Solution);
	return result;
}

double GasInterfaceImpl::continuumCooling(const Solution& s) const
{
	return _freeFree->cooling(np_ne(s), s.T) +
	       Ionization::cooling(s.abundancev(inH), s.abundancev(inp), s.abundancev(ine), s.T);
}

double GasInterfaceImpl::continuumHeating(const Solution& s) const
{
	double result = _freeFree->heating(np_ne(s), s.T, s.specificIntensityv);
	result += Ionization::heating(s.abundancev(inp), s.abundancev(ine), s.T, _frequencyv,
	                              s.specificIntensityv);
	if (_molecularLevels)
		result += _molecularLevels->dissociationHeating(s.H2Solution);
	return result;
}
