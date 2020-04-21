#ifndef CORE_CHEMISTRY_HPP
#define CORE_CHEMISTRY_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "SpeciesIndex.hpp"
#include <map>
#include <vector>

namespace GasModule
{
    /** This class can store a list of reactions, and solve the chemical equilibrium based on the
        coefficients of the reactions and their reaction rates. */
    class Chemistry
    {
    public:
        virtual ~Chemistry() = default;

        /** Set up the internal SpeciesIndex. If one was already present, it is reset. */
        void registerSpecies(const std::vector<std::string>& namev);

        /** Add another species. Be sure to call prepareCoefficients when done adding species
            and/or reactions. */
        void addSpecies(const std::string& name);

        /** Read-only acces to the SpeciesIndex object, so the result vector can be interpreted,
            and optionally a SpeciesVector can be made. */
        const SpeciesIndex& speciesIndex() const { return _speciesIndex; }

        /** Function to provide a clear syntax for adding reactions in the setup of the chemical
            network. Each reaction is given a number, and the reaction is added to the reaction
            index using the given name as a key. The last, optional argument allows the user to
            change the density power of the reactants (by default, it is equal to the
            stoichiometric coefficient). Subclasses typically implement a constructor which
            calls addReaction and prepareCoefficients. Make sure to call prepareCoefficients
            when done adding species and/or reactions.*/
        void addReaction(const std::string& reactionName, const std::vector<std::string>& reactantNamev,
                         const std::vector<int>& reactantStoichv, const std::vector<std::string>& productNamev,
                         const std::vector<double>& productStoichv, const std::vector<int>& reactantPowerv);

        /** Version without the optional argument. */
        void addReaction(const std::string& reactionName, const std::vector<std::string>& reactantNamev,
                         const std::vector<int>& reactantStoichv, const std::vector<std::string>& productNamev,
                         const std::vector<double>& productStoichv);

        /** This should be called after reactions have been added */
        void prepareCoefficients();

        /** Look up the index of a reaction, based on the name that was given in to addReaction. */
        int reactionIndex(const std::string& reactionName) const;

        /** The number of reactions. This determines the size of the reaction rate vector */
        int numReactions() const { return _reactionv.size(); }

        /** The number of species. This determines the size of the species density vector */
        int numSpecies() const { return _numSpecies; }

        /** Solves the chemical network given a certain rate coefficient vector (indexed on the
            reactions). An initial value n0v can be given. A vector containing updated densities is
            returned. */
        EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v, double maxTime = -1) const;

        /** Evaluate the rate of change for each species [cm-3 s-1]. A vector for the total
            reaction speeds is also given (rateCoeffv * density product). It will be used as a
            workspace, and can also be used to diagnose the speed of each reaction [s-1]. It needs
            to be of size _numReactions. The output for Fv and the current nv are passed by
            pointer, because their memory is owned by two gsl_vector instances. An
            Eigen::Map<EVector> is used in the implementation (such an object can not simply be
            passed as an argument, as there are lots of weird semantics involved). */
        void evaluateFv(double* FvOutput, const double* nv, const EVector& rateCoeffv, EVector& kv) const;

        /** Evaluate the Jacobian of Fv [s-1]. Every column j is the derivative of Fv towards n_j.
            For optimization (and for diagnostics), a matrix to put the jacobian of the reaction
            speed vector is passed. It needs to be of size (_numReactions, _numSpecies). */
        void evaluateJvv(double* JvvOutputRowMajor, const double* nv, const EVector& rateCoeffv, EMatrix& Jkvv) const;

    private:
        /** Solve the chemistry by evolving the system until equilibrium. */
        EVector solveTimeDep(const EVector& rateCoeffv, const EVector& n0v, double maxTime) const;

        /** Calculate the density factor needed to calculate the speed of reaction r. Formula:
            Product_i q n_i ^ Rir, where n_i are the elements of nv, and Rir = _rStoichvv(i, r). */
        double reactionSpeed(const double* nv, const EVector& rateCoeffv, size_t r) const;

        /** Calculate the derivative of the density product for reaction r with respect to the
            density j. Formula: (Rjr - 1) * n_j^{Rjr - 1} * Product_{i != j} n_i ^ Rir */
        void reactionSpeedJacobian(EMatrix& Jkvv, const double* nv, const EVector& rateCoeffv) const;

        // Keep track of index for each species name.
        SpeciesIndex _speciesIndex;

        class Reaction
        {
        public:
            /** See documentation of Chemistry::addReaction(). Create a new reaction with custom
                exponents for the reactant densities. The last argument needs to have the same
                size as rNamev and rCoeffv. Exponents are integer since it's faster and
                currently we don't need real ones. */
            Reaction(const std::vector<std::string>& rNamev, const std::vector<int>& rCoeffv,
                     const std::vector<std::string>& pNamev, const std::vector<double>& pCoeffv,
                     const std::vector<int>& rPowerv);

            const std::vector<std::string> _rNamev, _pNamev;
            const std::vector<int> _rCoeffv;
            const std::vector<double> _pCoeffv;
            const std::vector<int> _rPowerv;
        };

        std::vector<Reaction> _reactionv;
        std::map<std::string, int> _reactionIndexm;

        // Filled in by prepareCoefficients(), indexed on (species, reaction)
        EMatrix_int _rStoichvv;  // reactant coefficients
        EMatrix_int _rPowervv;   // reactant density power in reaction speed
        EMatrix _netStoichvv;    // product - reactant coefficients
        size_t _numSpecies;
        int _numReactions;
    };
}
#endif  // CORE_CHEMISTRY_HPP
