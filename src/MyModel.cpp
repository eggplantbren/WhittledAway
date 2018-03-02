#include "MyModel.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace WhittledAway
{

MyModel::MyModel()
:y(N)
,C(N, N)
{

}

void MyModel::generate(InfoNest::RNG& rng)
{
    A = exp(0.1*rng.randn());
    log10_period = log10(N) - rng.rand();
    quality = exp(log(1.0) + log(1000.0)*rng.rand());

    calculate_C();
    generate_data(rng);
    calculate_logl();
}

void MyModel::calculate_C()
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

double MyModel::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    // Perturb one of the three parameters
    int which = rng.rand_int(2);
    if(which == 0)
    {
        A = log(A);
        logH -= -0.5*pow(A/0.1, 2);
        A += 0.1*rng.randh();
        logH += -0.5*pow(A/0.1, 2);
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
    calculate_C();
    calculate_logl();

    return logH;
}

void MyModel::generate_data(InfoNest::RNG& rng)
{
    Eigen::VectorXd n(N);
    for(size_t i=0; i<N; ++i)
        n(i) = rng.randn();
    y = Lmat*n;
}

void MyModel::calculate_logl(bool whittle)
{
    if(whittle)
    {
        // TODO: Implement Whittle likelihood
        logl = 0.0;

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
        logl += -0.5*y.dot(L.solve(y));
    }

    if(std::isnan(logl) || std::isinf(logl))
        logl = -1E300;
}

double MyModel::perturb(InfoNest::RNG& rng)
{
    double logH = -logl;

    logH += perturb_parameters(rng);
    calculate_logl();

    logH += logl;

    return logH;
}

void MyModel::print(std::ostream& out) const
{
    out<<A<<' '<<log10_period<<' '<<quality<<' ';
    for(size_t i=0; i<N; ++i)
        out << y[i] << ' ';
}

double MyModel::parameter_distance(const MyModel& s1, const MyModel& s2)
{
    return std::abs(s2.log10_period - s1.log10_period);
}

} // namespace WhittledAway

