#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>

#include "CavityHamiltonianGenerator.h"
#include "GapRatioCalculator.h"
#include "FockBaseGenerator.h"
#include "Simulation.h"
#include "Utils.h"
#include "Parameters.h"

namespace {
    class UniformGenerator {
    private:
        std::mt19937 generator;
        std::uniform_real_distribution<double> distribution;

    public:
        UniformGenerator(double min, double max, unsigned long seed) : generator(seed), distribution(min, max) { }

        UniformGenerator(UniformGenerator &dummy) = delete;
        UniformGenerator operator=(UniformGenerator dummy) = delete;

        double operator()() {
            return this->distribution(this->generator);
        }
    };

    template<typename HamiltonianGenerator_t>
    class OnsiteDisorderAveragingModel {
    public:
        static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, std::size_t simulationIndex,
                                              std::size_t numberOfSimulations)
        {
            static_cast<void>(simulationIndex);
            static_cast<void>(numberOfSimulations);
            hamiltonianGenerator.resampleOnsiteEnergy();
        }
    };
}

int main(int argc, char **argv) {
    if (argc != 2)
        die(std::string("Usage: ") + argv[0] + " [input file]");

    std::string filename(argv[1]);
    std::ifstream input(filename);
    if (!input)
        die("Cannot open " + filename + " to read input parameters");

    Parameters params(input);
    params.print(std::cout);
    std::cout << std::endl;

    FockBaseGenerator baseGenerator;
    FockBase base = baseGenerator.generate(params.numberOfSites, params.numberOfBosons);

    using TheHamiltonianGenerator = CavityHamiltonianGenerator<UniformGenerator>;
    auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W, params.seed);
    TheHamiltonianGenerator::Parameters hamiltonianParams;
    hamiltonianParams.J = params.J;
    hamiltonianParams.U = params.U;
    hamiltonianParams.U1 = params.U1;
    hamiltonianParams.beta = params.beta;
    hamiltonianParams.phi0 = params.phi0;
    auto hamiltonianGenerator = std::make_unique<TheHamiltonianGenerator>(base, hamiltonianParams,
                                                                          std::move(disorderGenerator),
                                                                          params.usePeriodicBC);
    using TheAveragingModel = OnsiteDisorderAveragingModel<TheHamiltonianGenerator>;
    using TheSimulation = Simulation<TheHamiltonianGenerator, TheAveragingModel>;
    TheSimulation simulation(std::move(hamiltonianGenerator), params.numberOfSimulations, 0.5, 0.1);

    simulation.perform(std::cout);
    Quantity meanGapRatio = simulation.getMeanGapRatio();
    meanGapRatio.separator = Quantity::Separator::PLUS_MINUS;
    std::cout << "Mean gap ratio: " << meanGapRatio << std::endl;

    return EXIT_SUCCESS;
}