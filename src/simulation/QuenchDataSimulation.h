//
// Created by pkua on 30.08.2020.
//

#ifndef MBL_ED_QUENCHDATASIMULATION_H
#define MBL_ED_QUENCHDATASIMULATION_H

#include <utility>

#include "evolution/CorrelationsTimeEvolution.h"
#include "simulation/SimulationsSpan.h"
#include "evolution/ChebyshevEvolver.h"
#include "core/HamiltonianGenerator.h"
#include "core/AveragingModel.h"
#include "core/RND.h"
#include "core/QuenchCalculator.h"

/**
 * @brief A class performing a series of quenches and gathering information about them.
 * @details Namely, it calculates average energy of quench, sample error (sqrt of variance of energy distribution)
 * and average quantum energy error. The end results are available via getResultsFields().
 * Template parameters are only for testing purposes.
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel,
         typename QuenchCalculator_t = QuenchCalculator>
class QuenchDataSimulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> initialHamiltonianGenerator;
    std::unique_ptr<RND> initialRnd;
    std::unique_ptr<HamiltonianGenerator_t> finalHamiltonianGenerator;
    std::unique_ptr<RND> finalRnd;

    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<QuenchCalculator_t> quenchCalculator;
    SimulationsSpan simulationsSpan;

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
                         std::unique_ptr<QuenchCalculator_t> quenchCalculator, const SimulationsSpan &simulationsSpan)
            : initialHamiltonianGenerator{std::move(initialHamiltonianGenerator)}, initialRnd{std::move(initialRnd)},
              finalHamiltonianGenerator{std::move(finalHamiltonianGenerator)}, finalRnd{std::move(finalRnd)},
              averagingModel{std::move(averagingModel)}, quenchCalculator{std::move(quenchCalculator)},
              simulationsSpan{simulationsSpan}
    {
        Expects(this->simulationsSpan.total > 0);
        Expects(this->simulationsSpan.from < this->simulationsSpan.to);
        Expects(this->simulationsSpan.to <= this->simulationsSpan.total);
    }

    /**
     * @brief Performs the quenches for simulation ranges given by SimulaitonsSpan from the constructor.
     * @details Quench is performed by finding the ground state of the initial hamiltonian and calculating its energy
     * and spread for the final hamiltonian. Averages are calculates according to AveragingModel passed in the
     * constructor.
     */
    void perform(std::ostream &logger) {
        this->quenchCalculator->clear();
        for (std::size_t i = this->simulationsSpan.from; i < this->simulationsSpan.to; i++) {
            arma::wall_clock timer;
            logger << "[QuenchDataSimulator::perform] Performing quench " << i << "... " << std::flush;
            timer.tic();

            this->averagingModel->setupHamiltonianGenerator(*this->initialHamiltonianGenerator, *this->initialRnd, i,
                                                            this->simulationsSpan.total);
            this->averagingModel->setupHamiltonianGenerator(*this->finalHamiltonianGenerator, *this->finalRnd, i,
                                                            this->simulationsSpan.total);
            arma::sp_mat initialHamiltonian = this->initialHamiltonianGenerator->generate();
            arma::sp_mat finalHamiltonian = this->finalHamiltonianGenerator->generate();

            this->quenchCalculator->addQuench(initialHamiltonian, finalHamiltonian);

            logger << "done (" << timer.toc() << " s). epsilon: " << this->quenchCalculator->getLastQuenchEpsilon();
            logger << "; quantum error: " << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty();
            logger << std::endl;
        }
    }

    [[nodiscard]] std::vector<std::string> getResultsHeader() const {
        return {"epsilon", "avgError", "quantumError"};
    }

    [[nodiscard]] std::vector<std::string> getResultsFields() const {
        return {std::to_string(this->quenchCalculator->getMeanEpsilon()),
                std::to_string(this->quenchCalculator->getEpsilonAveragingSampleError()),
                std::to_string(this->quenchCalculator->getMeanEpsilonQuantumUncertainty())};
    }
};


#endif //MBL_ED_QUENCHDATASIMULATION_H
