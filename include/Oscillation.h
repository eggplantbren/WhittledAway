#ifndef WhittledAway_Oscillation_h
#define WhittledAway_Oscillation_h

// Includes
#include <armadillo>
#include <Eigen/Dense>
#include <ostream>
#include <vector>
#include "Model.h"
#include "RNG.h"

namespace WhittledAway
{

class Oscillation : public Model
{
    private:

        // Noise sd
        static constexpr double sigma = 1.0;

    private:

        // Amplitude, log-period, and quality
        double A, log10_period, quality;

        // Calculate log likelihood
        void calculate_C();
        void calculate_logl();

        // Helper for perturb
        double perturb_parameters(InfoNest::RNG& rng);

    public:

        // Constructor needs to make C of the right size
        Oscillation();

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

    public:
        // A few options to use for `distance`
        static double parameter_distance
                            (const Oscillation& s1, const Oscillation& s2);

};

} // namespace WhittledAway

#endif

