// Includes
#include <ctime>
#include <iostream>
#include <type_traits>
#include "Execute.hpp"
#include "RNG.h"
#include "Oscillation.h"

// MAIN
int main()
{
    // What example are we using? The class and the distance function.
    using TheModel = WhittledAway::Oscillation;

    // Make sure TheModel is derived from WhittledAway::Model
    static_assert(std::is_base_of<WhittledAway::Model, TheModel>::value,
                    "Examples must derive from WhittledAway::Model.");

    const auto& dist_func = TheModel::parameter_distance;

    // Create random number generators
    // The first one is used to generate reference points
    // and the second is used for the Nested Sampling
    unsigned long seed0 = 0;
    unsigned long seed1 = time(0);
    InfoNest::RNG rng0(seed0);
    InfoNest::RNG rng1(seed1);

    // Define run parameters
    constexpr double depth         = 30.0;
    constexpr size_t num_reps      = 1000;
    constexpr size_t num_particles = 10;
    constexpr size_t mcmc_steps    = 1000;

    // Do the run.
    InfoNest::execute<TheModel>(rng0, rng1, depth, num_reps, num_particles,
                                mcmc_steps, dist_func,
                                InfoNest::Mode::conditional_entropy,
                                10000);

    return 0;
}

