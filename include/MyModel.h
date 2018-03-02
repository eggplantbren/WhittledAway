#ifndef WhittledAway_MyModel_h
#define WhittledAway_MyModel_h

// Includes
#include <Eigen/Dense>
#include <ostream>
#include <vector>
#include "RNG.h"

namespace WhittledAway
{

class MyModel
{
    private:

        // Number of data points
        static constexpr size_t N = 100;

        // Noise sd
        static constexpr double sigma = 1.0;

    private:

        // Amplitude, log-period, and quality
        double A, log10_period, quality;

        // Data
        Eigen::VectorXd y;

        // Covariance matrix and its Cholesky decomposition
        Eigen::MatrixXd C;
        Eigen::LLT<Eigen::MatrixXd> L;
        Eigen::MatrixXd Lmat;

        // Log likelihood, i.e. ln p(y | params).
        // Only needed for one of the proposals
        double logl;

        // Calculate log likelihood
        void calculate_C();
        void generate_data(InfoNest::RNG& rng);
        void calculate_logl(bool whittle=false);

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

