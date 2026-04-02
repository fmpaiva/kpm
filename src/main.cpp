#include <chrono> // for std::chrono functions
#include <filesystem>
#include <iostream>
#include <fstream>
#include "Constants.h"
#include "simulation.h"
#include "kpm.h"

class Timer {
private:
    // Type aliases to make accessing nested type easier
    using Clock = std::chrono::steady_clock;
    using Second = std::chrono::duration<double, std::ratio<1>>;

    std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
    void reset() {
        m_beg = Clock::now();
    }

    [[nodiscard]] double elapsed() const {
        return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
    }
};

int main() {
    int Lx = 256;
    int Ly = 256;
    int nu = 8;
    bool y_periodic = false;
    bool disorder = true;
    const Hamiltonian hamiltonian{Lx, Ly, nu, true, y_periodic, disorder};
    std::filesystem::path path{"data/ldos-" + std::to_string(Lx) + "x" + std::to_string(Ly)
        + "-" + std::to_string(nu) + "_L-y_periodic" + std::to_string(y_periodic) + "-disorder"
        + std::to_string(disorder) + "-LL.dat"};

    Timer t;
    // Simulation::dos(hamiltonian, 1024, path);
    Simulation::ldos(hamiltonian, -3.86, 0.0005, path);
    std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

    // Timer t;
    // auto moments = KPM::gaussian_chebyshev_moments(0, 0.000001);
    // std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
    //
    // std::filesystem::path path{"data/moments.dat"};
    // std::ofstream file{path};
    //
    // if (!file) {
    //     std::cout << "Error opening file\n";
    //     exit(1);
    // }
    //
    // for (long i = 0; i < moments.size(); ++i) {
    //     file << moments(i) << "\n";
    // }
}

// Meter constantes configuráveis num ficheiro para trocar facilmente

// Variação média quadrática do potencial no espaço todo da ideia da desordem.
// Função espectral: Elemento de matriz do delta no espaço dos k
