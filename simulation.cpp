#include <complex>
#include <iostream>

#include "hamiltonian.h"
#include "random.h"
#include "simulation.h"

#include <fstream>

#include "Constants.h"
#include "kpm.h"
#include "Eigen/Dense"

using namespace std::complex_literals;


void Simulation::dos(const Hamiltonian& h, long int N_pol, const std::filesystem::path& path) {
    if (N_pol % 2 != 0)
        throw std::invalid_argument("N_pol must be even");

    KPM kpm(h);
    Eigen::ArrayXd moments(N_pol, 1);
    moments.setZero();

    for (Index i = 0; i < Constants::N_samples_DOS; ++i) {
        kpm.fill_random_phase();

        double current_moment_0  = kpm.vector().squaredNorm();
        moments(0) += current_moment_0;

        kpm.KPM_iteration();
        double current_moment_1 = kpm.prev_vector().dot(kpm.vector()).real();
        moments(1) += current_moment_1;

        for (Index n = 1; n < N_pol / 2; ++n) {
            kpm.KPM_iteration();

            moments(2 * n) += 2 * kpm.prev_vector().squaredNorm() - current_moment_0;
            moments(2 * n + 1) += 2 * kpm.vector().dot(kpm.prev_vector()).real() - current_moment_1;
        }
    }
    moments /= Constants::N_samples_DOS * static_cast<double>(h.dimension());

    std::ofstream file{path};
    if (!file)
        throw std::runtime_error("Could not open " + path.string() + ".\n");

    file << "mu\n";
    for (Index i = 0; i < N_pol; ++i) {
        file << moments(i) << "\n";
    }
}

// TODO: Document in notes and comments
// TODO: Ver cálculo dos momentos e número de polinómios, código henrique
void Simulation::ldos(const Hamiltonian& h, double E, double sigma, long int N_pol, const std::filesystem::path& path) {
    if (N_pol % 2 != 0)
        throw std::invalid_argument("N_pol must be even");

    KPM kpm(h);
    Eigen::ArrayXd gaussian_moments =
        std::pow(2.0, 0.75) * std::pow(M_PI, 0.25) * std::sqrt(sigma * h.scale()) *
        KPM::gaussian_chebyshev_moments(E * h.scale(), std::sqrt(2) * sigma * h.scale(), N_pol);
    Eigen::ArrayXd ldos(h.dimension());
    Eigen::ArrayXd var(h.dimension());
    ldos.setZero();
    var.setZero();

    Eigen::ArrayXcd running_map(h.dimension());
    for (Index i = 0; i < Constants::N_samples_LDOS; ++i) { // TODO: Parallelize
        kpm.fill_random_cplx_gaussian();
        Eigen::ArrayXcd starting_random_vector(kpm.vector());
        running_map.setZero();
        kpm.accumulate(N_pol, gaussian_moments, running_map);

        running_map *= starting_random_vector.conjugate();
        Eigen::ArrayXd running_ldos = running_map.cwiseAbs2();

        auto N = static_cast<double>(i);
        const Eigen::ArrayXd mean = i == 0 ? ldos : ldos / N;
        var = (var + (mean - running_ldos).square() / (N + 1)) *  N / (N + 1);

        ldos += running_ldos;
    }
    ldos /= Constants::N_samples_LDOS;
    var /= Constants::N_samples_LDOS;

    std::ofstream file{path};
    if (!file)
        throw std::runtime_error("Could not open " + path.string() + ".\n");

    file << "ldos,variance,rel_err\n";
    for (Index i = 0; i < h.dimension(); ++i) {
        file << ldos(i) << "," << var(i) << "," << std::sqrt(var(i)) / ldos(i) << "\n";
    }
}