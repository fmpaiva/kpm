#include "hamiltonian.h"
#include <complex>

#include "random.h"

using Eigen::Index;
using namespace std::complex_literals;

Hamiltonian::Hamiltonian(int Lx, int Ly, int nu, bool x_periodic, bool y_periodic, bool disorder)
    : scale_{0.23}
    , D_{Lx * Ly}
    , Lx_{Lx}
    , Ly_{Ly}
    , tx_{std::complex(1., 0.0) * scale_}
    , ty_{std::complex(1., 0.0) * scale_}
    , flux_{static_cast<double>(nu) / Lx_}
    , x_periodic_{x_periodic}
    , y_periodic_{y_periodic}
    , disorder_(disorder)
    , arr_disorder_(D_)
{
    if (disorder_) {
        for (Index i = 0; i < D_; ++i) {
            arr_disorder_(i) = Random::uniform_double(-0.3, 0.3) * scale_;
        }
    } else {
        arr_disorder_.setZero();
    }
}

[[nodiscard]] Eigen::VectorXcd Hamiltonian::act_on(const Eigen::Ref<const Eigen::VectorXcd>& psi) const {
    Eigen::VectorXcd out{psi.size()};
    out.setZero();

    for (Index y = 0; y < Ly_; ++y) {
        const std::complex phase{std::exp(2i * M_PI * flux_ * static_cast<double>(y))};

        for (Index x = 0; x < Lx_; ++x) {
            const Index i = x + Lx_ * y;

            // +x neighbor
            if (x + 1 < Lx_) {
                out(i) += -tx_ * phase * psi((x + 1) + Lx_ * y);
            } else if (x_periodic_) {
                out(i) += -tx_ * phase * psi(Lx_ * y);  // x=0
            }

            // -x neighbor
            if (x > 0) {
                out(i) += -tx_ * std::conj(phase) * psi((x - 1) + Lx_ * y);
            } else if (x_periodic_) {
                out(i) += -tx_ * std::conj(phase) * psi((Lx_ - 1) + Lx_ * y);
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

    // Meter desordem 0
    if (disorder_) {
        for (Index i = 0; i < D_; ++i) {
            out(i) += arr_disorder_(i) * psi(i);
        }
    }

    return out;
}