#include "SimpleHChemistry.hpp"
#include "DebugMacros.hpp"
#include "Ionization.hpp"

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

    /* H2 formation on grain surfaces H + H -> H2 Need to put 1 here for the H
	   stoichiometry, because the rate scales with nH instead of nH^2 (see rateCoeffv()). By
	   then using twice (not the word 'half' in the label) the reaction rate, we end up with
	   an equal tempo of H2 formation, but one that scales only linearly with nH. */

    // TODO: find a way around the 0.5 coefficient, so that I can use integer coefficients
    addReaction("half H2 formation", {"H"}, {1}, {"H2"}, {.5});

    prepareCoefficients();
}

EVector SimpleHChemistry::rateCoeffv(double T, const Spectrum& specificIntensity, double kDissFromH2Levels,
                                     double kH2FormationGrain) const
{
    EVector k(numReactions());
    k(reactionIndex("H photoionization")) = Ionization::photoRateCoeff(specificIntensity);
    k(reactionIndex("H collisional ionization")) = Ionization::collisionalRateCoeff(T);
    k(reactionIndex("H radiative recombination")) = Ionization::recombinationRateCoeff(T);
    k(reactionIndex("H2 dissociation")) = kDissFromH2Levels;
    /* The rate of this reaction is twice the H2 formation rate. See comment in constructor
	   implementation. */
    k(reactionIndex("half H2 formation")) = 2 * kH2FormationGrain;

    DEBUG("0 photo-ion; 1 coll-ion; 2 rad-rec; 3 dissoc; 4 H2form\n" << k << '\n');

    return k;
}
