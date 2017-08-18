#ifndef GASMODULE_GIT_SRC_GRAINSPECIES_H_
#define GASMODULE_GIT_SRC_GRAINSPECIES_H_

/** This class groups together some functions that provide grain properties that are needed
    throughout the code. This way, we can provide extensibility, while keeping the photoelectric
    heating and the grain surface H2 formation rate implementations manageable. This class should
    probably be 'static', as in, all the members being static functions. */
class GrainSpecies
{
public:
	GrainSpecies();

	// See 2013-Röllig et al. table C.1.
	typedef struct SurfaceInteractionParameters
	{

		SurfaceInteractionParameters() = default;
		SurfaceInteractionParameters(double EH2, double Es, double EHp, double EHc,
		                             double aSqrt, double nuH2, double nuHc, double F)
		                : _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt},
		                  _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
		{
		}
		const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0},
		                _f{0};
	} SurfaceInteractionParameters;

	static virtual SurfaceInteractionParameters sfcInteractionPar()
	{
		if (t == GrainType::CAR)
			return {520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10};
		else /* if  (t == GrainType::SILICATE) */
			return {320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10};
	}
};

#endif /* GASMODULE_GIT_SRC_GRAINSPECIES_H_ */
