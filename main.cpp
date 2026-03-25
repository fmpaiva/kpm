#include <chrono> // for std::chrono functions
#include <filesystem>
#include <iostream>
#include "simulation.h"

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
    int nu = 21;
    double E1 = 0.65;
    double E2 = 1.;
    const Hamiltonian hamiltonian{Lx, Ly, nu, true, false, true};
    std::filesystem::path output_path1{"../data/dos-" + std::to_string(Lx) + "x" + std::to_string(Ly)
        + "-" + std::to_string(nu) + "_L-yopen.txt"};
    std::filesystem::path output_path2{"../data/ldos-" + std::to_string(E1) + "-" + std::to_string(Lx)
        + "x" + std::to_string(Ly) + "-" + std::to_string(nu) + "_L-yopen.txt"};
    std::filesystem::path output_path3{"../data/ldos-" + std::to_string(E2) + "-" + std::to_string(Lx)
        + "x" + std::to_string(Ly) + "-" + std::to_string(nu) + "_L-yopen.txt"};

    Timer t;
    Simulation::dos(hamiltonian, static_cast<long int>(std::pow(2, 14)), output_path1);
    Simulation::ldos(hamiltonian, E1, 0.009, static_cast<long>(std::pow(2, 9)), output_path2);
    Simulation::ldos(hamiltonian, E2, 0.009, static_cast<long>(std::pow(2, 9)), output_path3);

    std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
}

// Ver campos magnéticos mais fracos para aproximar resultados do contínuos
// Ver como isto evolui com o número de polinómios. Reduzir sigma e ir aumentando número de polinómios.
// Variação média quadrática do potencial no espaço todo da ideia da desordem.
// Ver variância on the fly. Temos de controlar as incertezas numéricas!
// Outra incerteza tem a ver com a resolução
// Ver para cada sigma que aquilo converge com o número de polinómios. Depois temos de variar os sigmas
// Ver se somar LDOS dá DOS
// Podemos também calcular com o cálculo não estocástico
// Função espectral: Elemento de matriz do delta no espaço dos k