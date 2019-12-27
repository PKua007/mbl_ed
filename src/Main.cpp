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
        UniformGenerator(double min, double max, unsigned long seed) : generator(seed), distribution(min, max) { }

        UniformGenerator(UniformGenerator &dummy) = delete;
        UniformGenerator operator=(UniformGenerator dummy) = delete;

        double operator()() {
            return this->distribution(this->generator);
        }
    };

    class Parameters {
    public:
        std::size_t numberOfSites{};
        std::size_t numberOfBosons{};
        double J{};
        double W{};
        double U{};
        double U1{};
        bool usePeriodicBC{};
        std::size_t numberOfSimulations{};
        std::size_t seed{};

        Parameters(int argc, char **argv) {
            if (argc != 10) {
                std::cerr << "Usage: " << argv[0] << " [number of sites] [number of bosons] [J] [W] [U] [U1] ";
                std::cerr << "[use periodic BC true/false] [number of simulations] [seed]" << std::endl;
                exit(EXIT_FAILURE);
            }

            this->numberOfSites = std::stoul(argv[1]);
            this->numberOfBosons = std::stoul(argv[2]);
            this->J = std::stod(argv[3]);
            this->W = std::stod(argv[4]);
            this->U = std::stod(argv[5]);
            this->U1 = std::stod(argv[6]);
            this->usePeriodicBC = (std::string(argv[7]) == "true");
            this->numberOfSimulations = std::stoul(argv[8]);
            this->seed = std::stoul(argv[9]);

            Validate(numberOfSites > 0);
            Validate(numberOfBosons > 0);
            Validate(J >= 0);
            Validate(W >= 0);
            Validate(U >= 0);
            Validate(U1 >= 0);
            Validate(numberOfSimulations > 0);
        }

        void print(std::ostream &out) const {
            out << "number of sites       : " << this->numberOfSites << std::endl;
            out << "number of bosons      : " << this->numberOfBosons << std::endl;
            out << "J                     : " << this->J << std::endl;
            out << "W                     : " << this->W << std::endl;
            out << "U                     : " << this->U << std::endl;
            out << "U1                    : " << this->U1 << std::endl;
            out << "usePeriodicBC         : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
            out << "number of simulations : " << this->numberOfSimulations << std::endl;
            out << "seed                  : " << this->seed << std::endl;
        }
    };
}

int main(int argc, char **argv) {
    Parameters params(argc, argv);
    params.print(std::cout);
    std::cout << std::endl;

    FockBaseGenerator baseGenerator;
    FockBase base = baseGenerator.generate(params.numberOfSites, params.numberOfBosons);

    using TheHamiltonianGenerator = CavityHamiltonianGenerator<UniformGenerator>;
    auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W, params.seed);
    auto hamiltonianGenerator = std::make_unique<TheHamiltonianGenerator>(base, params.J, params.U, params.U1,
                                                                          std::move(disorderGenerator),
                                                                          params.usePeriodicBC);
    Simulation<TheHamiltonianGenerator> simulation(std::move(hamiltonianGenerator), params.numberOfSimulations, 0.5,
                                                   0.1);

    simulation.perform(std::cout);
    Quantity meanGapRatio = simulation.getMeanGapRatio();
    meanGapRatio.separator = Quantity::Separator::PLUS_MINUS;
    std::cout << "Mean gap ratio: " << meanGapRatio << std::endl;

    return EXIT_SUCCESS;
}