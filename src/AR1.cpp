#include "AR1.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>

namespace WhittledAway
{

void AR1::generate(InfoNest::RNG& rng)
{
    beta = exp(-10.0 + 20.0*rng.rand());
    L = exp(log(0.1) + log(10.0*N)*rng.rand());

    generate_data(rng);
    calculate_logl();
}

void AR1::generate_data(InfoNest::RNG& rng)
{
    double alpha = exp(-1.0/L);
    double sigma = beta/sqrt(1.0 - alpha*alpha);
    y[0] = sigma*rng.randn();
    for(size_t i=1; i<N; ++i)
        y[i] = alpha*y[i-1] + beta*rng.randn();

    // Copy into the Armadillo vector
    for(size_t i=0; i<N; ++i)
        y_fft[i] = y[i];

    // Take the fft
    y_fft = arma::fft(y_fft);
}

void AR1::calculate_C()
{

}

double AR1::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(2);

    if(which == 0)
    {
        beta = log(beta);
        beta += 20.0*rng.randh();
        InfoNest::wrap(beta, -10.0, 10.0);
        beta = exp(beta);
    }
    else
    {
        L = log(L);
        L += log(10.0*N)*rng.randh();
        InfoNest::wrap(L, log(0.1), log(N));
        L = exp(L);
    }

    calculate_logl();
    return logH;
}


void AR1::calculate_logl()
{
    if(whittle)
    {
        logl = 0.0;

        double model_psd, data_pgram, w;
        double alpha = exp(-1.0/L);
        double alpha_sq = alpha*alpha;
        double sigma = beta/sqrt(1.0 - alpha*alpha);
        double sigma_sq = sigma*sigma;

        // Loop over first half of fft
        for(size_t j=1; j<=(size_t)floor(0.5*(N-1)); ++j)
        {
            // Angular frequency
            w = 2.0*M_PI*static_cast<double>(j)/N;

            // Model PSD
            model_psd = sigma_sq*(1.0 - alpha_sq)
                            / (1.0 + alpha_sq - 2.0*alpha*cos(w)); 

            // Data periodogram
            data_pgram = (pow(y_fft[j].real(), 2)
                                    + pow(y_fft[j].imag(), 2)) / N;

            // Whittle
            logl += -log(model_psd) - data_pgram/model_psd;
        }
    }
    else
    {
        // Exact likelihood in the time domain
        logl = 0.0;

        double C = -0.5*log(2*M_PI);
        double alpha = exp(-1.0/L);
        double sigma = beta/sqrt(1.0 - alpha*alpha);
        double log_beta = log(beta);
        double prec = 1.0/(beta*beta);

        logl += C - log(sigma) - 0.5*pow(y[0]/sigma, 2);
        for(size_t i=1; i<N; ++i)
        {
            logl += C - log_beta
                        - 0.5*prec*pow(y[i] - alpha*y[i-1], 2);
        }
    }

    if(std::isnan(logl) || std::isinf(logl))
        logl = -1E300;
}

void AR1::print(std::ostream& out) const
{
    out << beta << ' ' << L;
}

double AR1::parameter_distance(const AR1& s1, const AR1& s2)
{
    double dsq = 0.0;
    dsq += pow(log(s1.beta) - log(s2.beta), 2);
    dsq += pow(log(s1.L) - log(s2.L), 2);
    return sqrt(dsq);
}

} // namespace WhittledAway

