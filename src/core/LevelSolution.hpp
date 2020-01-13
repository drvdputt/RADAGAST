#ifndef CORE_LEVELSOLUTION_HPP
#define CORE_LEVELSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

class LevelCoefficients;

/** Stores non-constant quantities that have to do with a level system, and claculates derived
    quantities */
class LevelSolution
{
public:
    /** Pass a pointer to the correct LevelCoefficients object, for access to the constants.
	    Make sure that the arguments passed to setCvv and setNv are dimensionally compatible
	    with the given LevelCoefficients. A temperature and density are also needed to
	    calculate the line width and normalization. Once these have been set externally,
	    using the setters, the other functions can be called safely. */
    LevelSolution(const LevelCoefficients* lc) : _levelCoefficients{lc} {};

    /** Set new temperature */
    void setT(double t) { _t = t; }
    double t() const { return _t; }

    /** Update the collision coefficients */
    void setCvv(const EMatrix& cvv)
    {
        _isCvvSet = true;
        _cvv = cvv;
    }

    /** Update the level populations */
    void setNv(const EVector& nv)
    {
        _isNvSet = true;
        _nv = nv;
    }
    const EVector& nv() const { return _nv; }

    /** The spectrum emitted by the line transitions, expressed as the emission coefficient
	    j_nu f * (erg/cm3/s/Hz). */
    Array emissivityv(const Array& eFrequencyv) const;

    /** The line opacity alpha_nu, equivalent to kappaRho for dust [cm-1]. */
    Array opacityv(const Array& oFrequencyv) const;

    /** Net heating due to (de-)excitation [erg s-1 cm-3] */
    double netHeating() const;

    /** True if nv has been filled in */
    bool isNvSet() const { return _isNvSet; }

    /** True if cvv has been filled in */
    bool isCvvSet() const { return _isCvvSet; }

private:
    const LevelCoefficients* _levelCoefficients;
    double _t{0.};
    EVector _nv;
    EMatrix _cvv;
    bool _isNvSet{false};
    bool _isCvvSet{false};
};

#endif  // CORE_LEVELSOLUTION_HPP
