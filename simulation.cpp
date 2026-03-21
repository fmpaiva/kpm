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
    Eigen::Array<double, -1, 1> moments(N_pol, 1);
    moments.setZero();

    for (Index i = 0; i < Constants::N_samples; ++i) {
        kpm.fill_random_vector();

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
    moments /= Constants::N_samples * static_cast<double>(h.dimension());

    std::ofstream file{path};
    if (!file)
        throw std::runtime_error("Could not open " + path.string() + ".\n");

    file << "mu\n";
    for (Index i = 0; i < N_pol; ++i) {
        file << moments(i) << "\n";
    }
}