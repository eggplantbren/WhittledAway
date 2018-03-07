#ifndef WhittledAway_WhiteNoise_h
#define WhittledAway_WhiteNoise_h

// Includes
#include <armadillo>
#include <Eigen/Dense>
#include <ostream>
#include <vector>
#include "Model.h"
#include "RNG.h"

namespace WhittledAway
{

class WhiteNoise : public Model
{
    private:

        // Noise sd
        double sigma;

        void generate_data(InfoNest::RNG& rng);
        void calculate_C();
        void calculate_logl();

        // Helper for perturb
        double perturb_parameters(InfoNest::RNG& rng);

    public:

        // Generate from the distribution
        void generate(InfoNest::RNG& rng);

        // Printing to stream
        void print(std::ostream& out) const;

    public:
        // A few options to use for `distance`
        static double parameter_distance
                            (const WhiteNoise& s1, const WhiteNoise& s2);

};

} // namespace WhittledAway

#endif

