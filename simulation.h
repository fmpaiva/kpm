#ifndef SIMULATION_H
#define SIMULATION_H

#include <filesystem>
#include "hamiltonian.h"

class Simulation {
public:
    static void dos(const Hamiltonian& h, long int N_pol, const std::filesystem::path& output_file) ;
};

#endif //SIMULATION_H
