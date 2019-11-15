//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_SIMULATION_H
#define MBL_ED_SIMULATION_H


#include <iosfwd>
#include <random>
#include <memory>
#include <armadillo>

#include "Quantity.h"
#include "GapRatioCalculator.h"

template<typename ConcreteHamiltonianGenerator, typename ConcreteGapRatioCalculator = GapRatioCalculator>
class Simulation {
private:
    std::unique_ptr<ConcreteHamiltonianGenerator> hamiltonianGenerator;
    std::size_t numberOfSimulations{};
    double relativeMiddleEnergy{};
    double relativeMargin{};

    Quantity meanGapRatio{};

    std::vector<double> performSingleSimulation() {
        this->hamiltonianGenerator->resampleOnsiteEnergies();
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();
        arma::vec armaEigenenergies = arma::eig_sym(hamiltonian);

        std::vector<double> eigenenergies;
        eigenenergies.reserve(armaEigenenergies.size());
        std::copy(armaEigenenergies.begin(), armaEigenenergies.end(), std::back_inserter(eigenenergies));

        return eigenenergies;
    }

public:
    Simulation(std::unique_ptr<ConcreteHamiltonianGenerator> hamiltonianGenerator, size_t numberOfSimulations,
               double relativeMiddleEnergy, double relativeMargin)
            : hamiltonianGenerator(std::move(hamiltonianGenerator)), numberOfSimulations{numberOfSimulations},
              relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
    { }

    [[nodiscard]] Quantity getMeanGapRatio() const {
        return meanGapRatio;
    }

    void perform(std::ostream &logger) {
        ConcreteGapRatioCalculator calculator(this->relativeMiddleEnergy, this->relativeMargin);

        for (std::size_t i = 0; i < this->numberOfSimulations; i++) {
            logger << "[Simulation::perform] Performing simulation " << i << "... " << std::flush;
            auto eigenenergies = this->performSingleSimulation();
            calculator.addEigenenergies(eigenenergies);
            logger << "done." << std::endl;
        }

        this->meanGapRatio = calculator.calculateMean();
    }
};


#endif //MBL_ED_SIMULATION_H
