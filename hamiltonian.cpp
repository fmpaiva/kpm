#include "hamiltonian.h"
#include <complex>

using Eigen::Index;
using namespace std::complex_literals;

Hamiltonian::Hamiltonian(int Lx, int Ly, int nu, bool x_periodic, bool y_periodic)
    : D_{Lx * Ly}
    , Lx_{Lx}
    , Ly_{Ly}
    , tx_{std::complex(0.25, 0.0)}
    , ty_{std::complex(0.25, 0.0)}
    , flux_{static_cast<double>(nu) / Lx_}
    , x_periodic_{x_periodic}
    , y_periodic_{y_periodic}
{

}

[[nodiscard]] Eigen::VectorXcd Hamiltonian::act_on(const Eigen::Ref<const Eigen::VectorXcd>& psi) const {
    Eigen::VectorXcd out{psi.size()};
    out.setZero();

    for (Index y = 0; y < Ly_; ++y) {
        const std::complex exp{std::exp(2i * M_PI * flux_ * static_cast<double>(y))};

        for (Index x = 0; x < Lx_; ++x) {
            const Index i = x + Lx_ * y;

            // +x neighbor
            if (x + 1 < Lx_) {
                out(i) += -tx_ * exp * psi((x + 1) + Lx_ * y);
            } else if (x_periodic_) {
                out(i) += -tx_ * exp * psi(Lx_ * y);  // x=0
            }

            // -x neighbor
            if (x > 0) {
                out(i) += -tx_ * std::conj(exp) * psi((x - 1) + Lx_ * y);
            } else if (x_periodic_) {
                out(i) += -tx_ * std::conj(exp) * psi((Lx_ - 1) + Lx_ * y);
            }

            // +y neighbor
            if (y + 1 < Ly_) {
                out(i) += -ty_ * psi(x + Lx_ * (y + 1));
            } else if (y_periodic_) {
                out(i) += -ty_ * psi(x);  // y=0
            }

            // -y neighbor
            if (y > 0) {
                out(i) += -ty_ * psi(x + Lx_ * (y - 1));
            } else if (y_periodic_) {
                out(i) += -ty_ * psi(x + Lx_ * (Ly_ - 1));
            }
        }
    }

    return out;
}