#ifndef KPM_H
#define KPM_H

#include "Eigen/Dense"
#include <cstdint>
#include "hamiltonian.h"

using Eigen::Index;

class KPM {
public:
    explicit KPM(const Hamiltonian&);
    void KPM_iteration();
    static Eigen::Array<double, -1, 1> gaussian_chebyshev_moments(double mu, double sigma, long N);

    template <typename Sampler>
    void fill_random_vector(const Sampler& sampler) {
        for (Index i = 0; i < h_.dimension(); ++i) {
            v_(i, 0) = sampler();
        }

        index_ = 0;
    }

    [[nodiscard]] constexpr auto vector() const {
        return v_.col(index_ % 2);
    }
    [[nodiscard]] constexpr auto prev_vector() const {
        return v_.col((index_ + 1) % 2);
    }

private:
    const Hamiltonian& h_;
    Eigen::Matrix<std::complex<double>, -1, 2> v_;
    Index index_;
};

#endif