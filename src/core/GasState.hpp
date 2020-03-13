#ifndef CORE_GASSTATE_HPP
#define CORE_GASSTATE_HPP

#include <string>
#include <valarray>

namespace GasModule
{
    class GasState;
    class GasInterface;
    class GasInterfaceImpl;
    // forward declare friend function
    namespace Testing
    {
        void writeGasState(const std::string&, const GasModule::GasInterface&, const GasModule::GasState&);
    }

    /** This class is used to store a minimal set of results at the end of the equilibrium
        calculation, since it is far to memory-consuming to store every detail. A client is
        supposed to store these objects (one for each cell), and can choose to pass these to
        certain functions of the @c GasInterface to obtain the desired physical properties.
        Functions taking a GasState as an argument will recalculate the requested quantities, if
        necessary. */
    class GasState
    {
        friend class GasInterfaceImpl;
        friend void Testing::writeGasState(const std::string&, const GasModule::GasInterface&,
                                           const GasModule::GasState&);

    public:
        /** Set all the members of the gas state to new values. This function should not be used by
            clients, besides for testing purposes. */
        void setMembers(double t, const std::valarray<double>& nv, double n2s)
        {
            _t = t;
            _nv = nv;
            _n2s = n2s;
        }

        /** Returns the gas temperature in K. */
        double temperature() const { return _t; }

        /** Return density at a certain index. Abstract for now, but maybe a solution will be
            available in GasInterface later. */
        double density(int i) const { return _nv[i]; }

    private:
        double _t{0.};
        std::valarray<double> _nv{4};
        double _n2s{0.};
    };
}

#endif  // CORE_GASSTATE_HPP
