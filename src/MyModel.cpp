#include "MyModel.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>

namespace WhittledAway
{

MyModel::MyModel()
:t(N)
,mu(N)
,y(N)
{
    for(size_t i=0; i<N; ++i)
    {
        // Evenly spaced points at integer times
        t[i] = i;
    }
}

void MyModel::generate(InfoNest::RNG& rng)
{
    A = exp(0.1*rng.randn());
    log10_period = log10(N) - rng.rand();
    phi = 2.0 * M_PI * rng.rand();

    calculate_mu();

    for(size_t i=0; i<N; ++i)
        y[i] = mu[i] + sigma*rng.randn();

    calculate_logl();
}

void MyModel::calculate_mu()
{
    double omega = 2*M_PI/pow(10.0, log10_period); // Angular frequency
    for(size_t i=0; i<N; ++i)
        mu[i] = A * sin(omega*t[i] + phi);
}

double MyModel::perturb_parameters(InfoNest::RNG& rng)
{
    double logH = 0.0;

    // Perturb one of the three parameters
    int which = rng.rand_int(3);
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
        phi += 2 * M_PI * rng.randh();
        InfoNest::wrap(phi, 0.0, 2 * M_PI);
    }

    return logH;
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
        double C = log(1.0 / sqrt(2.0 * M_PI) / sigma);
        double tau = 1.0 / (sigma * sigma);

        for(size_t i=0; i<N; ++i)
            logl += C - 0.5 * pow(y[i] - mu[i], 2) * tau;
    }
}

double MyModel::perturb(InfoNest::RNG& rng)
{
    double logH = 0.0;

    int proposal_type = 1;//rng.rand_int(4);

    if(proposal_type == 0)
    {
        // Perturb parameters, changing data along with it
        std::vector<double> mu_old = mu;

        logH += perturb_parameters(rng);
        calculate_mu();

        double n;
        for(size_t i=0; i<N; ++i)
        {
            n = (y[i] - mu_old[i]) / sigma;
            y[i] = mu[i] + sigma * n;
        }

        calculate_logl();
    }
    else if(proposal_type == 1)
    {
        // Perturb parameters, keeping data constant
        // (aka Metropolis step of the posterior!)
        logH -= logl;

        logH += perturb_parameters(rng);

        calculate_mu();
        calculate_logl();

        logH += logl;        
    }
    else if(proposal_type == 2)
    {
        // Just change one datum
        int which = rng.rand_int(N);

        logH -= -0.5*pow((y[which] - mu[which])/sigma, 2);
        y[which] += sigma * rng.randh();
        logH += -0.5*pow((y[which] - mu[which])/sigma, 2);

        calculate_logl();
    }
    else
    {
        // Potentially regenerate many of the data points
        int reps = pow(N, rng.rand());
        int which;
        for(int i=0; i<reps; ++i)
        {
            which = rng.rand_int(N);
            y[which] = mu[which] + sigma * rng.randn();
        }

        calculate_logl();
    }

    return logH;
}

void MyModel::print(std::ostream& out) const
{
    out<<A<<' '<<log10_period<<' '<<phi<<' ';
    for(double yy: y)
        out<<yy<<' ';
}

double MyModel::parameter_distance(const MyModel& s1, const MyModel& s2)
{
    return std::abs(s2.log10_period - s1.log10_period);
}

} // namespace WhittledAway

