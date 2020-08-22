//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLUTION_H
#define MBL_ED_CHEBYSHEVEVOLUTION_H

#include <utility>

#include "CorrelationsTimeEvolution.h"
#include "simulation/SimulationsSpan.h"
#include "ChebyshevEvolver.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/AveragingModel.h"
#include "simulation/RND.h"
#include "QuenchCalculator.h"

/**
 * @brief A class performing a whole series of time evolutions using Chebyshev expansion technique.
 * @details See CorrelationsTimeEvolution and its "slave" classes to see what is calculated. The hamiltonian generator
 * is prepared for each simulation according to a specific averaging model and, if desired, the quench stated is
 * prepared for evolution (see constructor). The template parameters default to standard classes and exist solely for
 * mocking purposes. See the default classes description for details of what they do.
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel,
         typename CorrelationsTimeEvolution_t = CorrelationsTimeEvolution,
         typename QuenchCalculator_t = QuenchCalculator, typename ChebyshevEvolver_t = ChebyshevEvolver>
class ChebyshevEvolution {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution;

    SimulationsSpan simulationsSpan;

    std::unique_ptr<QuenchCalculator_t> quenchCalculator;
    std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator;
    std::unique_ptr<RND> quenchRnd;


    auto prepareHamiltonianAndPossiblyQuenchVector(std::size_t simulationIndex, std::ostream &logger) const {
        std::vector<arma::cx_vec> additionalVectors;

        this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, simulationIndex,
                                                        this->simulationsSpan.totalSimulations);
        arma::sp_mat hamiltonian = this->hamiltonianGenerator->generate();

        if (this->quenchCalculator != nullptr) {
            this->averagingModel->setupHamiltonianGenerator(*this->quenchHamiltonianGenerator, *this->quenchRnd,
                                                            simulationIndex, this->simulationsSpan.totalSimulations);
            arma::sp_mat initialHamiltonian = this->quenchHamiltonianGenerator->generate();

            this->quenchCalculator->addQuench(initialHamiltonian, hamiltonian);

            // Yup, Armadillo wouldn't let you initialize arma::cx_vec by arma::vec easily, we have to do this nonsense
            const arma::vec &quenchedStateDouble = this->quenchCalculator->getLastQuenchedState();
            arma::cx_vec quenchedState(quenchedStateDouble.size());
            std::copy(quenchedStateDouble.begin(), quenchedStateDouble.end(), quenchedState.begin());
            additionalVectors.push_back(quenchedState);

            logger << "quenched state epsilon: " << this->quenchCalculator->getLastQuenchEpsilon();
            logger << "; quantum error: " << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty() << ". ";
        }
        return std::make_pair(hamiltonian, additionalVectors);
    }

public:
    /**
     * @brief Constructor, where HamiltonianGenerator and RND for quench should be passed, or set to nullptr if
     * quench should not be done.
     * @details If quench is to be done, @a parameters.initialVectors has to have exatcly one external vector slot
     * provided, 0 otherwise.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                       const SimulationsSpan &simulationsSpan,
                       std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution,
                       std::unique_ptr<QuenchCalculator_t> quenchCalculator,
                       std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator,
                       std::unique_ptr<RND> quenchRnd)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, averagingModel{std::move(averagingModel)},
              rnd{std::move(rnd)}, correlationsTimeEvolution{std::move(correlationsTimeEvolution)},
              simulationsSpan{simulationsSpan}, quenchCalculator{std::move(quenchCalculator)},
              quenchHamiltonianGenerator{std::move(quenchHamiltonianGenerator)}, quenchRnd{std::move(quenchRnd)}
    {
        Expects(this->simulationsSpan.totalSimulations > 0);
        Expects(this->simulationsSpan.from < this->simulationsSpan.to);
        Expects(this->simulationsSpan.to <= this->simulationsSpan.totalSimulations);

        if (this->quenchCalculator == nullptr) {
            Expects(this->correlationsTimeEvolution->countExternalVectors() == 0);
        } else {
            Expects(this->quenchHamiltonianGenerator != nullptr);
            Expects(this->quenchRnd != nullptr);
            Expects(this->correlationsTimeEvolution->countExternalVectors() == 1);
        }
    }

    /**
     * @brief Constructor with no quench.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                       const SimulationsSpan &parameters,
                       std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution)
            : ChebyshevEvolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                 parameters, std::move(correlationsTimeEvolution), nullptr, nullptr, nullptr)
    { }

    /**
     * @brief Performs the simuations. The results (see CorrelationsTimeEvolution) are stored afterwards using
     * FileOstreamProvider from the constructor.
     * @details The name of the file is in the code, go check ;).
     */
    void perform(std::ostream &logger) {
        if (this->quenchCalculator != nullptr)
            this->quenchCalculator->clear();

        for (std::size_t i = this->simulationsSpan.from; i < this->simulationsSpan.to; i++) {
            arma::wall_clock wholeTimer;
            arma::wall_clock timer;

            wholeTimer.tic();
            logger << "[ChebyshevEvolution::perform] Performing evolution " << i << "... " << std::endl;
            logger << "[ChebyshevEvolution::perform] Preparing hamiltonian... " << std::flush;
            timer.tic();
            auto [hamiltonian, additionalVectors] = this->prepareHamiltonianAndPossiblyQuenchVector(i, logger);
            logger << "done (" << timer.toc() << " s)" << std::endl;

            logger << "[ChebyshevEvolution::perform] Preparing evolver... " << std::endl;
            timer.tic();
            ChebyshevEvolver_t evolver(hamiltonian, logger);
            logger << "[ChebyshevEvolution::perform] Preparing evolver done (" << timer.toc() << " s)." << std::endl;
            this->correlationsTimeEvolution->addEvolution(evolver, logger, additionalVectors);
            logger << "[ChebyshevEvolution::perform] Whole evolution took " << wholeTimer.toc() << " s." << std::endl;
        }

        if (this->quenchCalculator != nullptr) {
            logger << "[ChebyshevEvolution::perform] Mean quench data: epsilon: ";
            logger << this->quenchCalculator->getMeanEpsilon() << "; avg error: ";
            logger << this->quenchCalculator->getEpsilonAveragingSampleError() << "; quantum error: ";
            logger << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty() << std::endl;
        }
    }

    void storeResults(std::ostream &out) const {
        this->correlationsTimeEvolution->storeResult(out);
    }
};


#endif //MBL_ED_CHEBYSHEVEVOLUTION_H
