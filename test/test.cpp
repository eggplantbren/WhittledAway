// Includes
#include <iostream>
#include "Models.h"
#include "RNG.h"

using namespace WhittledAway;

// MAIN
int main()
{
    // What example are we using? The class and the distance function.
    using TheModel = Oscillation;

    // Make sure TheModel is derived from WhittledAway::Model
    static_assert(std::is_base_of<Model, TheModel>::value,
                    "Examples must derive from WhittledAway::Model.");

    // Create random number generator
    InfoNest::RNG rng(1);

    // Make a particle
    TheModel o1;
    o1.generate(rng);

    // Print the particle
    o1.print(std::cout);
    std::cout << std::endl;

    // Do MCMC
    for(int i=0; i<10000000; ++i)
    {
        // Make a proposal
        auto o2 = o1;
        double logH = o2.perturb(rng);

        // Accept?
        if(rng.rand() <= exp(logH))
            o1 = o2;

        if((i+1)%100 == 0)
        {
            std::cerr << (i+1) << std::endl;
            o1.print(std::cout);
            std::cout << std::endl;
        }
}


    return 0;
}

