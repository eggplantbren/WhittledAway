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
    mu = 1000.0*rng.randn();
    beta = exp(-10.0 + 20.0*rng.rand());
    L = exp(-10.0 + 20.0*rng.rand());

    generate_data(rng);
    calculate_logl();
}

void AR1::generate_data(InfoNest::RNG& rng)
{
    double alpha = exp(-1.0/L);
    double sigma = beta/sqrt(1.0 - alpha*alpha);
    y[0] = mu + sigma*rng.randn();
    for(size_t i=1; i<N; ++i)
        y[i] = mu + alpha*(y[i-1] - mu) + beta*rng.randn();
}

void AR1::calculate_C()
{

}

double AR1::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(3);

    if(which == 0)
    {
        logH -= -0.5*pow(mu/1000.0, 2);
        mu += 1000.0*rng.randh();
        logH += -0.5*pow(mu/1000.0, 2);
    }
    else if(which == 1)
    {
        beta = log(beta);
        beta += 20.0*rng.randh();
        InfoNest::wrap(beta, -10.0, 10.0);
        beta = exp(beta);
    }
    else
    {
        L = log(L);
        L += 20.0*rng.randh();
        InfoNest::wrap(L, -10.0, 10.0);
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

        double model_psd, data_psd, f;
        double alpha = exp(-1.0/L);
        double sigma = beta/sqrt(1.0 - alpha*alpha);

        // Loop over first half of fft
        for(size_t j=0; j<N/2; ++j)
        {
            f = static_cast<double>(j)/N;

            // Model PSD
            model_psd = sigma*sigma/(1.0 + alpha*alpha - 2*alpha*cos(2*M_PI*f));

            // Data PSD
            data_psd = (pow(y_fft[j].real(), 2) + pow(y_fft[j].imag(), 2));

            // Whittle
            logl += -log(model_psd) - data_psd/model_psd;
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

        logl += C - log(sigma) - 0.5*pow((y[0] - mu)/sigma, 2);
        for(size_t i=1; i<N; ++i)
        {
            logl += C - log_beta
                        - 0.5*prec*pow(y[i] - (mu + alpha*(y[i-1] - mu)), 2);
        }
    }

    if(std::isnan(logl) || std::isinf(logl))
        logl = -1E300;
}

void AR1::print(std::ostream& out) const
{
    out << mu << ' ' << beta << ' ' << L;
}

double AR1::parameter_distance(const AR1& s1, const AR1& s2)
{
    return std::abs(log(s2.L) - log(s1.L));
}

} // namespace WhittledAway

