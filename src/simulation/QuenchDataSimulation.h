//
// Created by pkua on 30.08.2020.
//

#ifndef MBL_ED_QUENCHDATASIMULATION_H
#define MBL_ED_QUENCHDATASIMULATION_H

#include <utility>

#include "evolution/TimeEvolution.h"
#include "simulation/SimulationsSpan.h"
#include "evolution/ChebyshevEvolver.h"
#include "core/HamiltonianGenerator.h"
#include "core/AveragingModel.h"
#include "core/RND.h"
#include "core/QuenchCalculator.h"
#include "RestorableSimulation.h"

/**
 * @brief A class performing quenches and gathering information about them.
 * @details Namely, it calculates average energy of quench, sample error (sqrt of variance of energy distribution)
 * and average quantum energy error. The end results are available via getResultsFields().
 * Template parameters are only for testing purposes.
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel,
         typename QuenchCalculator_t = QuenchCalculator>
class QuenchDataSimulation : public RestorableSimulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> initialHamiltonianGenerator;
    std::unique_ptr<RND> initialRnd;
    std::unique_ptr<HamiltonianGenerator_t> finalHamiltonianGenerator;
    std::unique_ptr<RND> finalRnd;

    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<QuenchCalculator_t> quenchCalculator;

public:
    /**
     * @brief Constructs the class.
     * @details Most of parameters are self-explainatory. We have initial and final hamiltonian generator and
     * corresponding random generators, which should be seeded with the same seed.
     */
    QuenchDataSimulation(std::unique_ptr<HamiltonianGenerator_t> initialHamiltonianGenerator,
                         std::unique_ptr<RND> initialRnd,
                         std::unique_ptr<HamiltonianGenerator_t> finalHamiltonianGenerator,
                         std::unique_ptr<RND> finalRnd, std::unique_ptr<AveragingModel_t> averagingModel,
                         std::unique_ptr<QuenchCalculator_t> quenchCalculator)
            : initialHamiltonianGenerator{std::move(initialHamiltonianGenerator)}, initialRnd{std::move(initialRnd)},
              finalHamiltonianGenerator{std::move(finalHamiltonianGenerator)}, finalRnd{std::move(finalRnd)},
              averagingModel{std::move(averagingModel)}, quenchCalculator{std::move(quenchCalculator)}
    { }

    [[nodiscard]] std::vector<std::string> getResultsHeader() const {
        return {"epsilon", "avgError", "quantumError"};
    }

    [[nodiscard]] std::vector<std::string> getResultsFields() const {
        return {std::to_string(this->quenchCalculator->getMeanEpsilon()),
                std::to_string(this->quenchCalculator->getEpsilonAveragingSampleError()),
                std::to_string(this->quenchCalculator->getMeanEpsilonQuantumUncertainty())};
    }

    void seedRandomGenerators(unsigned long seed) override {
        this->initialRnd->seed(seed);
        this->finalRnd->seed(seed);
    }

    /**
     * @brief Performs the quench for the simulation of index @a simulationIndex of out @a totalSimulations
     * @details Quench is performed by finding the ground state of the initial hamiltonian and calculating its energy
     * and spread for the final hamiltonian. Hamiltonian is prepared using AveragingModel passed in the constructor.
     */
    void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger &logger) override {
        Expects(totalSimulations > 0);
        Expects(simulationIndex < totalSimulations);

        arma::wall_clock timer;
        logger.verbose() << "Performing quench " << simulationIndex << " started... " << std::endl;
        timer.tic();

        this->averagingModel->setupHamiltonianGenerator(*this->initialHamiltonianGenerator, *this->initialRnd,
                                                        simulationIndex, totalSimulations);
        this->averagingModel->setupHamiltonianGenerator(*this->finalHamiltonianGenerator, *this->finalRnd,
                                                        simulationIndex, totalSimulations);
        arma::sp_mat initialHamiltonian = this->initialHamiltonianGenerator->generate();
        arma::sp_mat finalHamiltonian = this->finalHamiltonianGenerator->generate();

        this->quenchCalculator->addQuench(initialHamiltonian, finalHamiltonian);

        logger.info() << "Performing quench " << simulationIndex << " done (" << timer.toc() << " s). ";
        logger << "epsilon: " << this->quenchCalculator->getLastQuenchEpsilon();
        logger << "; quantum error: " << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty();
        logger << std::endl;
    }

    void storeState(std::ostream &binaryOut) const override {
        this->quenchCalculator->storeState(binaryOut);
    }

    void joinRestoredState(std::istream &binaryIn) override {
        this->quenchCalculator->joinRestoredState(binaryIn);
    }

    void clear() override {
        this->quenchCalculator->clear();
    }

    [[nodiscard]] std::string getTagName() const override {
        return "quench";
    }
};


#endif //MBL_ED_QUENCHDATASIMULATION_H
