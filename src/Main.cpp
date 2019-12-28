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
#include "Config.h"
#include "Utils.h"

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
        double beta{};
        double phi0{};
        bool usePeriodicBC{};
        std::size_t numberOfSimulations{};
        std::size_t seed{};

        Parameters(int argc, char **argv) {
            if (argc != 2) {
                std::cerr << "Usage: " << argv[0] << " [input file]" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::string filename(argv[1]);
            std::ifstream input(filename);
            if (!input)
                die("Cannot open " + filename + " to read input parameters");

            auto config = Config::parse(input, '=', true);
            for (const auto &key : config.getKeys()) {
                if (key == "numberOfSites")
                    this->numberOfSites = config.getUnsignedLong("numberOfSites");
                else if (key == "numberOfBosons")
                    this->numberOfBosons = config.getUnsignedLong("numberOfBosons");
                else if (key == "J")
                    this->J = config.getDouble("J");
                else if (key == "W")
                    this->W = config.getDouble("W");
                else if (key == "U")
                    this->U = config.getDouble("U");
                else if (key == "U1")
                    this->U1 = config.getDouble("U1");
                else if (key == "beta")
                    this->beta = config.getDouble("beta");
                else if (key == "phi0")
                    this->phi0 = config.getDouble("phi0");
                else if (key == "usePeriodicBC")
                    this->usePeriodicBC = config.getBoolean("usePeriodicBC");
                else if (key == "numberOfSimulations")
                    this->numberOfSimulations = config.getUnsignedLong("numberOfSimulations");
                else if (key == "seed")
                    this->seed = config.getUnsignedLong("seed");
                else
                    std::cerr << "[Parameters::Parameters] Warning: unknown parameter " << key << std::endl;
            }

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
            out << "beta                  : " << this->beta << std::endl;
            out << "phi0                  : " << this->phi0 << std::endl;
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
    TheHamiltonianGenerator::Parameters hamiltonianParams;
    hamiltonianParams.J = params.J;
    hamiltonianParams.U = params.U;
    hamiltonianParams.U1 = params.U1;
    hamiltonianParams.beta = params.beta;
    hamiltonianParams.phi0 = params.phi0;
    auto hamiltonianGenerator = std::make_unique<TheHamiltonianGenerator>(base, hamiltonianParams,
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