#ifndef SIMULATION_H
#define SIMULATION_H

#include <filesystem>
#include "hamiltonian.h"

class Simulation {
public:
    static void dos(const Hamiltonian& h, long int N_pol, const std::filesystem::path& output_file) ;
    static void ldos(const Hamiltonian& h, double E, double sigma, const std::filesystem::path& path);
    static void time_evolve(const Hamiltonian& h, 
                            const Eigen::ArrayXcd& wavepacket, 
                            double t, 
                            const std::filesystem::path& path
    );
};

#endif //SIMULATION_H
