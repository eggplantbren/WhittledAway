#ifndef WhittledAway_AR1_h
#define WhittledAway_AR1_h

// Includes
#include <armadillo>
#include <Eigen/Dense>
#include <ostream>
#include <vector>
#include "Model.h"
#include "RNG.h"

namespace WhittledAway
{

class AR1 : public Model
{
    private:

        // Mean level, beta, L
        double mu, beta, L;

        // Calculate log likelihood
        void calculate_C();
        void calculate_logl();

        // Helper for perturb
        double perturb_parameters(InfoNest::RNG& rng);

    public:

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Override
        void generate_data(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

    public:
        // A few options to use for `distance`
        static double parameter_distance
                            (const AR1& s1, const AR1& s2);

};

} // namespace WhittledAway

#endif

