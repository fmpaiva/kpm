#include <chrono> // for std::chrono functions
#include <filesystem>
#include <iostream>

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

#include "simulation.h"

int main() {
    int Lx = 256;
    int Ly = 256;
    int nu = 21;
    double E1 = 0.65 / 5;
    double E2 = 1. / 5;
    const Hamiltonian hamiltonian{Lx, Ly, nu, true, false, true};
    std::filesystem::path output_path1{"../data/dos-" + std::to_string(Lx) + "x" + std::to_string(Ly)
        + "-" + std::to_string(nu) + "_L-yopen.txt"};
    std::filesystem::path output_path2{"../data/ldos-" + std::to_string(E1) + "-" + std::to_string(Lx)
        + "x" + std::to_string(Ly) + "-" + std::to_string(nu) + "_L-yopen.txt"};
    std::filesystem::path output_path3{"../data/ldos-" + std::to_string(E2) + "-" + std::to_string(Lx)
        + "x" + std::to_string(Ly) + "-" + std::to_string(nu) + "_L-yopen.txt"};

    Timer t;
    Simulation::dos(hamiltonian, static_cast<long int>(std::pow(2, 14)), output_path1);
    Simulation::ldos(hamiltonian, E1, 0.003, static_cast<long>(std::pow(2, 12)), output_path2);
    Simulation::ldos(hamiltonian, E2, 0.003, static_cast<long>(std::pow(2, 12)), output_path3);

    std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
}