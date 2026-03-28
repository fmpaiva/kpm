#include <chrono> // for std::chrono functions
#include <filesystem>
#include <iostream>
#include <fstream>
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
    int Lx = 512;
    int Ly = 512;
    int nu = 1;
    bool y_periodic = true;
    bool disorder = false;
    const Hamiltonian hamiltonian{Lx, Ly, nu, true, y_periodic, disorder};
    std::filesystem::path path{"../data/dos-" + std::to_string(Lx) + "x" + std::to_string(Ly)
        + "-" + std::to_string(nu) + "_L-y_periodic" + std::to_string(y_periodic) + "-disorder"
        + std::to_string(disorder) + ".dat"};

    Timer t;
    Simulation::dos(hamiltonian, static_cast<long>(std::pow(2, 14)), path);
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

// Fazer estudo da convergência da aproximação à gaussiana para termos fórmula empírica fiável
// Ver output do henrique
// Meter constantes configuráveis num ficheiro para trocar facilmente

// Ver campos magnéticos mais fracos para aproximar resultados do contínuos
// Ver como isto evolui com o número de polinómios. Reduzir sigma e ir aumentando número de polinómios.
// Variação média quadrática do potencial no espaço todo da ideia da desordem.
// Outra incerteza tem a ver com a resolução
// Ver para cada sigma que aquilo converge com o número de polinómios. Depois temos de variar os sigmas
// Ver se somar LDOS dá DOS
// Podemos também calcular com o cálculo não estocástico
// Função espectral: Elemento de matriz do delta no espaço dos k
