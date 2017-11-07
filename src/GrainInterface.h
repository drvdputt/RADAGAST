#ifndef GASMODULE_GIT_SRC_GRAININTERFACE_H_
#define GASMODULE_GIT_SRC_GRAININTERFACE_H_

#include <functional>
#include <memory>
#include <valarray>
#include <vector>

class GrainType;

namespace GasModule
{

// NEEDED FOR H2 FORMATION //
/** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al. table
    C.1. */
typedef struct SfcInteractionPar
{
	SfcInteractionPar() = default;
	SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt, double nuH2,
	                  double nuHc, double F);
	/* This boolean is set to false when the default constructor is called, signifiying that the
	   created object does not contain any useful information (meaning that the graintype for
	   which it was contructed is not supported, and that the graintype given in the
	   corresponding static function should be skipped for the H2 formation rate calculation. */
	bool _valid{false};
	const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
} SfcInteractionPar;

/** Function signature for custom ionization potential. */
typedef std::function<double(double a, int Z)> IonizationPotentialf;

/** Function signature for custom photoelectric yield. */
typedef std::function<double(double a, int Z, double hnu)> PhotoelectricYieldf;

/** Function signature for custom autionization threshold. */
typedef std::function<double(double a)> AutoIonizationThresholdf;

/** Function signature for custom sticking coefficient. */
typedef std::function<double(double a, int Z, double z_i)> StickingCoefficientf;

/** List of grain types which have built-in values for these properties. */
enum class GrainTypeLabel
{
	CAR,
	SIL
};

/** Class that a client code should use to pass the grain properties in a cell. This class contains
    a list of grain models, for each of which a set of properties must be given. These properties,
    are listed in the nested 'Population' class. Another class, closely related to this one but not
    to be included by the client code, is GrainTypeProperties. It contains a bunch of functions that
    return parameters pertaining to the different choices of GrainType listed above. Those functions
    are in a separate file because they have no use in the public interface. */
class GrainInterface
{
public:
	/** The properties that need to be given per grain model that needs to be included. The
	    grain type can be one of the types listed in the enum above. Then, an array of sizes,
	    number densities, and one of temperatures needs to be specified (consider even using
	    grain temperature distributions (for each grain size!) in the future). The temperatures
	    are important for the H2 formation rate on the surfaces of the grains. The absorption
	    efficiency also needs to be given (need to review its exact definition). It represents
	    the number of photons absorbed from an ambient radiation field at a certain wavelength.
	    It needs to be given for each grain size, and for each point of the frequency grid. */
	class Population
	{
	public:
		/** For CAR or SIL grain type. The built-in values for some of the quantities for h2
		    formation and the photoelectric effect will be used. Remember to use move
		    construction, because of the unique_ptr member. TODO: make this more elegant or
		    something. Not using a move constructor will not compile. This is something I
		    should ask peter about. */
		Population(GrainTypeLabel type, const std::valarray<double>& sizev,
		           const std::valarray<double>& densityv,
		           const std::valarray<double>& temperaturev,
		           const std::vector<std::valarray<double>>& qAbsvv);

		/** For custom grain types. Some lambda functions need to be provided. */
		Population(const std::valarray<double>& sizev,
		           const std::valarray<double>& densityv,
		           const std::valarray<double>& temperaturev,
		           const std::vector<std::valarray<double>>& qAbsvv,
		           const SfcInteractionPar& sfcInteractionPar, bool heatingAvailable,
		           double workFunction, IonizationPotentialf ionizationpotentialf,
		           PhotoelectricYieldf photoElectricYieldf,
		           AutoIonizationThresholdf autoIonizationThresholdf,
		           StickingCoefficientf stickingCoefficientf);

		// Undelete the move constructor
		Population(Population&&);
		~Population();

		const GrainType* type() const { return _type.get(); }

		const std::valarray<double>&_sizev, _densityv, _temperaturev;
		const std::vector<std::valarray<double>>& _qAbsvv;

	private:
		/** Properties which should be provided either by calling the big constructor, or by
		    functions of the builtin grain types. In the latter case, the following pointer
		    will be initialized, and the correct functions will be called using
		    polymorphism. A builtin grain type instance will be assigned to this pointer
		    based on the "type" argument in the top constructor. */
		std::unique_ptr<GrainType> _type;
	};

	/** Constructor which takes a vector of predefined populations (by pointer). */
	GrainInterface(const std::vector<Population>* populationv);

	/** Default constructor, equivalent to no grains at all. */
	GrainInterface();

	~GrainInterface();

	int numPopulations() const;
	const Population* population(int i) const;
	const std::vector<Population>* populationv() const;

private:
	// Default constructor is empty vector
	const std::vector<Population>* _populationv;
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININTERFACE_H_ */
