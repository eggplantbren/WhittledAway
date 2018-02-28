#ifndef WhittledAway_MyModel_h
#define WhittledAway_MyModel_h

// Includes
#include <ostream>
#include <vector>
#include "RNG.h"

namespace WhittledAway
{

class MyModel
{
    private:

        // Number of data points and so on
        static constexpr size_t N = 101;

        // Noise sd
        static constexpr double sigma = 1.0;

    private:

        // Amplitude, log-period, and phase
        double A, log10_period, phi;

        // Model curve
        std::vector<double> t;
        std::vector<double> mu;

        // The data
        std::vector<double> y;

        // Log likelihood, i.e. ln p(y | params).
        // Only needed for one of the proposals
        double logl;

        // Calculate log likelihood
        void calculate_logl(bool whittle=false);

        // Compute mu from the current parameter values
        void calculate_mu();

        // Helper for perturb
        double perturb_parameters(InfoNest::RNG& rng);

    public:
        // Do-nothing constructor
        MyModel();

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Metropolis proposal
        double perturb(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

    public:
        // A few options to use for `distance`
        static double parameter_distance
                            (const MyModel& s1, const MyModel& s2);

};

} // namespace WhittledAway

#endif

