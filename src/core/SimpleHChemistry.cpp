#include "SimpleHChemistry.hpp"
#include "DebugMacros.hpp"
#include "Ionization.hpp"

namespace RADAGAST
{
    SimpleHChemistry::SimpleHChemistry()
    {
        registerSpecies({"e-", "H+", "H", "H2"});

        // CONSERVATION EQUATIONS
        // Protons
        // addConserved({"H+", "H", "H2"}, {1, 1, 2});

        // Electrons
        // addConserved({"e-", "H", "H2"}, {1, 1, 2});

        // REACTIONS
        // Photoionization
        // H + gamma -> e- + H+
        addReaction("H photoionization", {"H"}, {1}, {"e-", "H+"}, {1, 1});

        // Ionization by collision with electron
        // H + e- -> H+ + 2e-
        addReaction("H collisional ionization", {"H", "e-"}, {1, 1}, {"H+", "e-"}, {1, 2});

        // Radiative recombination
        // e- + H+ -> H + gamma
        addReaction("H radiative recombination", {"e-", "H+"}, {1, 1}, {"H"}, {1});

        // Dissociation after excitation
        // H2 -> H + H
        addReaction("H2 dissociation", {"H2"}, {1}, {"H"}, {2});

        // H2 formation on grain surfaces grain + 2H -> H2. Need to overide the power of nH by 1
        // (see last argument), otherwise the rate will scale as nH * nH.
        addReaction("H2 formation", {"H"}, {2}, {"H2"}, {1}, {1});

        prepareCoefficients();
    }

    EVector SimpleHChemistry::rateCoeffv(double T, const Spectrum& meanIntensity, double kDissFromH2Levels,
                                         double kH2FormationGrain) const
    {
        EVector k(numReactions());
        k(reactionIndex("H photoionization")) = Ionization::photoRateCoeff(meanIntensity);
        k(reactionIndex("H collisional ionization")) = Ionization::collisionalRateCoeff(T);
        k(reactionIndex("H radiative recombination")) = Ionization::recombinationRateCoeff(T);
        k(reactionIndex("H2 dissociation")) = kDissFromH2Levels;
        // The rate of this reaction is twice the H2 formation rate. See comment in constructor
        // implementation.
        k(reactionIndex("H2 formation")) = kH2FormationGrain;
        DEBUG("0 photo-ion; 1 coll-ion; 2 rad-rec; 3 dissoc; 4 H2form\n" << k << '\n');
        return k;
    }
}
