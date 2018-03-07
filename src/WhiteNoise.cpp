#include "WhiteNoise.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>

namespace WhittledAway
{

void WhiteNoise::generate(InfoNest::RNG& rng)
{
    sigma = exp(-10.0 + 20.0*rng.rand());

    generate_data(rng);
    calculate_logl();
}

void WhiteNoise::generate_data(InfoNest::RNG& rng)
{
    for(size_t i=0; i<N; ++i)
        y[i] = sigma*rng.randn();
}

void WhiteNoise::calculate_C()
{
    // Fill covariance matrix
    for(size_t i=0; i<N; ++i)
    {
        for(size_t j=i; j<N; ++j)
        {
            C(i, j) = 0.0;

            if(i != j)
                C(j, i) = C(i, j);
        }
    }

    for(size_t i=0; i<N; ++i)
        C(i, i) += sigma*sigma;

    L = C.llt();
    Lmat = L.matrixL();
}

double WhiteNoise::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    // Perturb noise level
    sigma = log(sigma);
    sigma += 20.0*rng.randh();
    InfoNest::wrap(sigma, -10.0, 10.0);
    sigma = exp(sigma);

//    if(!whittle)
//        calculate_C();
    calculate_logl();

    return logH;
}


void WhiteNoise::calculate_logl()
{
    if(whittle)
    {
        logl = 0.0;

        double model_psd, data_pgram;

        // Loop over first half of fft
        for(size_t j=1; j<=(size_t)floor(0.5*(N-1)); ++j)
        {
            // Model PSD
            model_psd = sigma*sigma;

            // Data periodogram
            data_pgram = (pow(y_fft[j].real(), 2) + pow(y_fft[j].imag(), 2)) / N;

            // Whittle
            logl += -log(model_psd) - data_pgram/model_psd;
        }
    }
    else
    {
        // Exact likelihood in the time domain
        logl = 0.0;

        logl += -0.5*log(2*M_PI*sigma*sigma)*N;
        double tau = 1.0/(sigma*sigma);
        for(size_t i=0; i<N; ++i)
            logl += -0.5*pow(y[i], 2)*tau;
    }

    if(std::isnan(logl) || std::isinf(logl))
        logl = -1E300;
}

void WhiteNoise::print(std::ostream& out) const
{
    out << sigma << ' ';
    for(size_t i=0; i<N; ++i)
        out << y[i] << ' ';
}

double WhiteNoise::parameter_distance(const WhiteNoise& s1, const WhiteNoise& s2)
{
    return std::abs(log(s2.sigma) - log(s1.sigma));
}

} // namespace WhittledAway

