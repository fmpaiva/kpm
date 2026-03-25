#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <Eigen/Dense>

class Hamiltonian {
public:
    Hamiltonian(int Lx, int Ly, int nu, bool x_periodic, bool y_periodic, bool disorder);

    [[nodiscard]] constexpr auto dimension() const { return D_;}

    [[nodiscard]] Eigen::VectorXcd act_on(const Eigen::Ref<const Eigen::VectorXcd>&) const;

private:
    long int D_;
    int Lx_;
    int Ly_;
    std::complex<double> tx_;
    std::complex<double> ty_;
    double flux_;
    bool x_periodic_;
    bool y_periodic_;
    bool disorder_;
    Eigen::VectorXcd arr_disorder_;
};

#endif