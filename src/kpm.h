#ifndef KPM_H
#define KPM_H

#include "hamiltonian.h"

using Eigen::Index;

class KPM {
public:
    explicit KPM(const Hamiltonian&);
    void KPM_iteration();
    void fill_random_phase();
    void fill_random_cplx_gaussian();
    void accumulate(long N_pol, const Eigen::ArrayXd& moments, Eigen::ArrayXcd& out);

    [[nodiscard]] auto vector() const {
        return v_.col(index_ % 2);
    }
    [[nodiscard]] auto prev_vector() const {
        return v_.col((index_ + 1) % 2);
    }

    static Eigen::ArrayXd gaussian_chebyshev_moments(double mu, double sigma);

private:
    const Hamiltonian& h_;
    Eigen::Matrix<std::complex<double>, -1, 2> v_; // TODO: Trocar para array
    Index index_;
};

#endif