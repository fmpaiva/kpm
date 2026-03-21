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
    void fill_random_vector();

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