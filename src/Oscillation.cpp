#include "Oscillation.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>

namespace WhittledAway
{

void Oscillation::generate(InfoNest::RNG& rng)
{
    A = exp(rng.randn());
    log10_period = log10(N) - rng.rand();
    quality = exp(log(1.0) + log(1000.0)*rng.rand());

    generate_data(rng);
    calculate_logl();
}

void Oscillation::calculate_C()
{
    // Fill covariance matrix
    double w0 = 2*M_PI/pow(10.0, log10_period);
    double eta = sqrt(std::abs(1.0 - 1.0/(4.0*quality*quality)));
    double S0 = A/w0/quality;
    double tau;

    for(size_t i=0; i<N; ++i)
    {
        for(size_t j=i; j<N; ++j)
        {
            tau = std::abs(static_cast<double>(i) - static_cast<double>(j));
            C(i, j) = cos(eta*w0*tau) + sin(eta*w0*tau)/(2.0*eta*quality);
            C(i, j) *= S0*w0*quality*exp(-w0*tau/(2.0*quality));

            if(i != j)
                C(j, i) = C(i, j);
        }
    }

    for(size_t i=0; i<N; ++i)
        C(i, i) += sigma*sigma;

    L = C.llt();
    Lmat = L.matrixL();
}

double Oscillation::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    // Perturb one of the three parameters
    int which = rng.rand_int(3);
    if(which == 0)
    {
        A = log(A);
        logH -= -0.5*pow(A/1.0, 2);
        A += 1.0*rng.randh();
        logH += -0.5*pow(A/1.0, 2);
        A = exp(A);
    }
    else if(which == 1)
    {
        log10_period += rng.randh();
        InfoNest::wrap(log10_period, log10(N) - 1.0, log10(N));
    }
    else
    {
        quality = log(quality);
        quality += log(1000.0)*rng.randh();
        InfoNest::wrap(quality, log(1.0), log(1000.0));
        quality = exp(quality);
    }

    if(!whittle)
        calculate_C();
    calculate_logl();

    return logH;
}


void Oscillation::calculate_logl()
{
    if(whittle)
    {

        // https://www4.stat.ncsu.edu/~reich/SpatialStats/code/Whittle.html
        logl = 0.0;

        double model_psd, data_psd, f, w;
        double w0 = 2*M_PI/pow(10.0, log10_period);
        double S0 = A/w0/quality;
        double coeff = sqrt(2.0/M_PI)*S0*pow(w0, 4);

        // Loop over first half of fft
        for(size_t j=0; j<N/2; ++j)
        {
            // Model PSD (Celerite paper, Eqn 20)
            f = static_cast<double>(j) / N;
            w = 2*M_PI*f;
            model_psd = coeff/(pow(w*w - w0*w0, 2)
                                        + w0*w0*w*w/(quality*quality));
            model_psd += sigma*sigma;

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

        logl += -0.5*log(2*M_PI)*N;

        // Log determinant
        double log_det = 0.0;
        for(size_t i=0; i<N; ++i)
            log_det += 2*log(Lmat(i, i));

        // Solve
        logl += -0.5*log_det - 0.5*y.dot(L.solve(y));
    }

    if(std::isnan(logl) || std::isinf(logl))
        logl = -1E300;
}

void Oscillation::print(std::ostream& out) const
{
    out<<A<<' '<<log10_period<<' '<<quality<<' ';
    for(size_t i=0; i<N; ++i)
        out << y[i] << ' ';
}

double Oscillation::parameter_distance(const Oscillation& s1, const Oscillation& s2)
{
    return std::abs(s2.log10_period - s1.log10_period);
}

} // namespace WhittledAway

