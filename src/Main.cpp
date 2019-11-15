#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>

#include "Assertions.h"
#include "CavityHamiltonianGenerator.h"
#include "GapRatioCalculator.h"
#include "FockBaseGenerator.h"
#include "Simulation.h"

namespace {
    class UniformGenerator {
    private:
        std::mt19937 generator;
        std::uniform_real_distribution<double> distribution;

    public:
        UniformGenerator(double min, double max, unsigned long seed) {
            this->generator.seed(seed);
            this->distribution = std::uniform_real_distribution<double>(min, max);
        };

        UniformGenerator(UniformGenerator &dummy) = delete;
        UniformGenerator operator=(UniformGenerator dummy) = delete;

        double operator()() {
            return this->distribution(this->generator);
        }
    };
}

int main(int argc, char **argv) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " [number of sites] [number of bosons] [J] [W] [U] [U1] ";
        std::cerr << "[number of simulations] [seed]" << std::endl;
        return EXIT_FAILURE;
    }

    std::size_t sites = std::stoul(argv[1]);
    std::size_t bosons = std::stoul(argv[2]);
    double J = std::stod(argv[3]);
    double W = std::stod(argv[4]);
    double U = std::stod(argv[5]);
    double U1 = std::stod(argv[6]);
    std::size_t numberOfSimulations = std::stoul(argv[7]);
    std::size_t seed = std::stoul(argv[8]);

    Validate(sites > 0);
    Validate(bosons > 0);
    Validate(J > 0);
    Validate(W > 0);
    Validate(U > 0);
    Validate(U1 > 0);
    Validate(numberOfSimulations > 0);

    FockBaseGenerator baseGenerator;
    FockBase base = baseGenerator.generate(sites, bosons);

    using TheDisorderGenerator = UniformGenerator;
    using TheHamiltonianGenerator = CavityHamiltonianGenerator<TheDisorderGenerator>;
    using TheSimulation = Simulation<TheHamiltonianGenerator, GapRatioCalculator>;

    auto disorderGenerator = std::make_unique<TheDisorderGenerator>(-W, W, seed);
    auto hamiltonianGenerator = std::make_unique<TheHamiltonianGenerator>(base, J, U, U1, std::move(disorderGenerator));
    TheSimulation simulation(std::move(hamiltonianGenerator), numberOfSimulations, 0.5, 0.1);

    simulation.perform(std::cout);
    Quantity result = simulation.getMeanGapRatio();
    result.separator = Quantity::Separator::PLUS_MINUS;
    std::cout << "Final result: " << result << std::endl;

    return EXIT_SUCCESS;
}