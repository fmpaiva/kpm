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

void KPM::fill_random_vector() {
    for (Index i = 0; i < h_.dimension(); ++i) {
        v_(i, 0) = std::exp(1.0i * Random::get(0, 2 * M_PI));
    }

    index_ = 0;
}