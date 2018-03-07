#include "Model.h"

namespace WhittledAway
{

Model::Model()
:y(N)
,y_fft(N)
,C(N, N)
{

}



double Model::perturb(InfoNest::RNG& rng)
{
    double logH = -logl;
    logH += perturb_parameters(rng);
    logH += logl;
    return logH;
}

void Model::generate_data(InfoNest::RNG& rng)
{
    calculate_C();

    Eigen::VectorXd n(N);
    for(size_t i=0; i<N; ++i)
        n(i) = rng.randn();
    y = Lmat*n;

    // Copy into the Armadillo vector
    for(size_t i=0; i<N; ++i)
        y_fft[i] = y[i];

    // Take the fft
    y_fft = arma::fft(y_fft);
}



} // namespace WhittledAway

