#ifndef WhittledAway_Model_h
#define WhittledAway_Model_h

#include <armadillo>
#include <cstdlib>
#include <Eigen/Dense>
#include "RNG.h"

namespace WhittledAway
{

/*
* A base class for different examples.
*/
class Model
{
    protected:

        // Number of data points
        static constexpr size_t N = 100;

        // Whittle likelihood?
        static constexpr bool whittle = false;

        // Log likelihood, i.e. ln p(y | params).
        // Only needed for one of the proposals
        double logl;

        // Data in an Eigen vector and an Armadillo vector
        Eigen::VectorXd y;
        arma::cx_vec y_fft;

        // Covariance matrix and its Cholesky decomposition
        Eigen::MatrixXd C;
        Eigen::LLT<Eigen::MatrixXd> L;
        Eigen::MatrixXd Lmat;

        void generate_data(InfoNest::RNG& rng);

        // Member functions required
        virtual void calculate_C() = 0;
        virtual void calculate_logl() = 0;
        virtual double perturb_parameters(InfoNest::RNG& rng) = 0;

        // Generate from the distribution
        virtual void generate(InfoNest::RNG& rng) = 0;

        // Printing to stream
        virtual void print(std::ostream& out) const = 0;

    public:

        // Constructor just sets up vectors and matrices to the
        // appropriate size
        Model();

        // Perturber
        double perturb(InfoNest::RNG& rng);
};


} // namespace WhittledAway

#endif

