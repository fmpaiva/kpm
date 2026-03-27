#include <chrono> // for std::chrono functions
#include <filesystem>
#include <iostream>
#include <fstream>

#include "kpm.h"
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
    int nu = 1;
    const Hamiltonian hamiltonian{Lx, Ly, nu, true, true, false};
    // std::filesystem::path output_path1{"../data/dos-" + std::to_string(Lx) + "x" + std::to_string(Ly)
    // + "-" + std::to_string(nu) + "_L-yopen.txt"};

    std::filesystem::path output_path1{"../data/moments.dat"};

    Timer t;
    double sigma = 0.0001;

    Eigen::ArrayXd moments = KPM::gaussian_chebyshev_moments(0.8, sigma);
    std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

    std::ofstream file{output_path1};
    if (!file)
        throw std::runtime_error("Could not open " + output_path1.string() + ".\n");

    file << "moments\n";
    for (Index i = 0; i < N_pol; ++i) {
        file << moments(i) << "\n";
    }
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