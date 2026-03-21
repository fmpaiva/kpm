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
    const Hamiltonian hamiltonian{2048, 2048, 37, true, true};
    std::filesystem::path output_path{"../1024x1024-34_L-xyopen.txt"};

    Timer t;
    Simulation::dos(hamiltonian, static_cast<long int>(std::pow(2, 12)), output_path);

    std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
}