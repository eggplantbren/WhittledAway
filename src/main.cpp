// Includes
#include <iostream>
#include <ctime>
#include "Execute.hpp"
#include "RNG.h"
#include "MyModel.h"

// MAIN
int main()
{
    // What example are we using? The class and the distance function.
    typedef WhittledAway::MyModel TheExample;
    const auto& dist_func = TheExample::parameter_distance;

    // Create random number generators
    // The first one is used to generate reference points
    // and the second is used for the Nested Sampling
    unsigned long seed0 = 0;
    unsigned long seed1 = time(0);
    InfoNest::RNG rng0(seed0);
    InfoNest::RNG rng1(seed1);

    // Define run parameters
    constexpr double depth         = 20.0;
    constexpr size_t num_reps      = 1000;
    constexpr size_t num_particles = 10;
    constexpr size_t mcmc_steps    = 1000;

    // Do the run.
    InfoNest::execute<TheExample>(rng0, rng1, depth, num_reps, num_particles,
                                  mcmc_steps, dist_func,
                                  InfoNest::Mode::conditional_entropy);

    return 0;
}

