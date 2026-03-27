#include <complex>
#include "kpm.h"
#include "random.h"

using namespace std::complex_literals;

KPM::KPM(const Hamiltonian& h)
    : h_(h)
    , v_(h.dimension(), 2)
    , index_(0)
{
    v_.setZero();
}

void KPM::KPM_iteration() {
    if (index_ == 0)
        v_.col(1) = h_.act_on(v_.col(0));
    else
        v_.col((index_ + 1) % 2) = 2 * h_.act_on(v_.col(index_ % 2)) - v_.col((index_ + 1) % 2);

    ++index_;
}

void KPM::fill_random_phase() {
    for (Index i = 0; i < h_.dimension(); ++i) {
        v_(i, 0) = std::exp(1i * Random::uniform_double(0,2 * M_PI));
    }

    index_ = 0;
}

void KPM::fill_random_cplx_gaussian() {
    for (Index i = 0; i < h_.dimension(); ++i) {
        v_(i, 0) = Random::gaussian_complex(0, 1) / std::sqrt(2.0);
    }

    index_ = 0;
}

void KPM::accumulate(long N_pol, const Eigen::ArrayXd& moments, Eigen::ArrayXcd& out) {
    for (Index n = 0; n < N_pol; ++n) {
        out += moments(n) * this->vector().array();
        this->KPM_iteration();
    }
}

Eigen::ArrayXd KPM::gaussian_chebyshev_moments(double mu, double sigma, long N) {
    Eigen::ArrayXd moments(N);

    for (Index n = 0; n < N; ++n) {
        const double one_minus_mu2 = 1.0 - mu * mu;
        const double sqrt_term = std::sqrt(one_minus_mu2);
        const double ratio = sigma / sqrt_term;   // sigma / sqrt(1 - mu^2)
        const double ratio2 = ratio * ratio;
        const double ratio4 = ratio2 * ratio2;

        const auto nn = static_cast<double>(n);
        const double n2 = nn * nn;
        const double n4 = n2 * n2;

        // a0
        const double a0 =
            1.0 - (-4.0 + 7.0 * mu * mu) * sigma * sigma * (3.0 - 6.0 * n2 * ratio2 + n4 * ratio4)
            / (24.0 * one_minus_mu2);

        // a1
        const double a1 =
            (nn * mu * sigma * sigma) / (2.0 * sqrt_term)
            * (3.0 - n2 * ratio2);

        // exponential prefactor
        const double prefactor =
            std::exp(-n2 * sigma * sigma / (2.0 * one_minus_mu2)) / sqrt_term;

        const double theta = std::acos(mu);

        moments(n) = prefactor * (a0 * std::cos(nn * theta) + a1 * std::sin(nn * theta));
    }

    return moments;
}
